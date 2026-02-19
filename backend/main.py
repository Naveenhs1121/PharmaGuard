"""
PharmaGuard — FastAPI Backend
Accepts VCF uploads + drug name(s), returns structured pharmacogenomic risk JSON.
"""

from fastapi import FastAPI, UploadFile, Form, File, HTTPException
from fastapi.middleware.cors import CORSMiddleware
import shutil
import os
import logging
from datetime import datetime, timezone
from typing import List, Optional
from vcf_parser import extract_variants
from risk_engine import predict_multi_drug, predict_risk, SUPPORTED_DRUGS
from explanation_engine import generate_explanation

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("PharmaGuard.API")

app = FastAPI(
    title="PharmaGuard API",
    description="Pharmacogenomic risk prediction aligned with CPIC guidelines.",
    version="2.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

TEMP_DIR = "/tmp/pharmaguard"
os.makedirs(TEMP_DIR, exist_ok=True)
MAX_FILE_SIZE_MB = 5


@app.get("/")
def home():
    return {
        "message": "PharmaGuard API v2.0 running",
        "supported_drugs": SUPPORTED_DRUGS,
        "supported_genes": ["CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"],
    }


@app.post("/analyze")
async def analyze(
    file: UploadFile = File(...),
    drug: str = Form(...),
):
    """Single-drug analysis endpoint (backwards compatible)."""
    return await _run_analysis(file, [drug])


@app.post("/analyze/multi")
async def analyze_multi(
    file: UploadFile = File(...),
    drugs: str = Form(...),   # comma-separated list
):
    """Multi-drug analysis endpoint. Pass drugs as comma-separated string."""
    drug_list = [d.strip() for d in drugs.split(",") if d.strip()]
    if not drug_list:
        raise HTTPException(status_code=400, detail="No drugs provided.")
    return await _run_analysis(file, drug_list)


async def _run_analysis(file: UploadFile, drugs: List[str]):
    """Shared analysis logic."""
    timestamp = datetime.now(timezone.utc).isoformat()

    # ── Validate file ────────────────────────────────────────────────────────
    if not file.filename.endswith(".vcf"):
        raise HTTPException(
            status_code=400,
            detail="Invalid file type. Only .vcf files are supported."
        )

    file_location = os.path.join(TEMP_DIR, f"temp_{file.filename}")
    content = await file.read()

    if len(content) > MAX_FILE_SIZE_MB * 1024 * 1024:
        raise HTTPException(
            status_code=413,
            detail=f"File too large. Maximum size is {MAX_FILE_SIZE_MB} MB."
        )

    with open(file_location, "wb") as buffer:
        buffer.write(content)

    # ── Parse VCF ────────────────────────────────────────────────────────────
    parse_success = True
    parse_error = None
    variants = []

    try:
        variants = extract_variants(file_location)
        logger.info(f"Parsed {len(variants)} pharmacogenomic variants from VCF")
    except Exception as e:
        parse_success = False
        parse_error = str(e)
        logger.error(f"VCF parsing failed: {e}")

    # ── Risk prediction ───────────────────────────────────────────────────────
    multi_result = predict_multi_drug(variants, drugs)

    # ── Explanation ───────────────────────────────────────────────────────────
    explanations = {}
    for drug_result in multi_result["drug_results"]:
        drug_name = drug_result["drug"]
        gene_profile = multi_result["gene_profiles"].get(drug_result["gene_used"], {})
        explanation = generate_explanation(
            variants=variants,
            drug=drug_name,
            risk={
                "label":    drug_result["risk_label"],
                "severity": drug_result["severity"],
                "confidence": drug_result["confidence_score"],
            },
            gene_profile=gene_profile,
            clinical_action=drug_result["clinical_action"],
        )
        explanations[drug_name] = explanation

    # ── Build structured response ─────────────────────────────────────────────
    drug_reports = []
    for drug_result in multi_result["drug_results"]:
        drug_name = drug_result["drug"]
        gene = drug_result["gene_used"]
        gene_data = multi_result["gene_profiles"].get(gene, {})

        drug_reports.append({
            "drug": drug_name,
            "timestamp": timestamp,

            "risk_assessment": {
                "risk_label":       drug_result["risk_label"],
                "severity":         drug_result["severity"],
                "confidence_score": drug_result["confidence_score"],
                "evidence_strength": drug_result["evidence_strength"],
            },

            "pharmacogenomic_profile": {
                "gene":               gene,
                "diplotype":          gene_data.get("diplotype", "Unknown"),
                "phenotype":          drug_result["phenotype"],
                "detected_variants":  gene_data.get("detected_rsids", []),
                "activity_score":     gene_data.get("total_activity_score", None),
            },

            "clinical_recommendation": {
                "action":           drug_result["clinical_action"],
                "cpic_guideline":   drug_result["cpic_guideline"],
                "supporting_variants": drug_result["supporting_variants"],
            },

            "reasoning": drug_result["reasoning"],
            "llm_explanation": explanations.get(drug_name, ""),
        })

    response = {
        "patient_id":  "PATIENT_DEMO",
        "timestamp":   timestamp,

        "drug_reports": drug_reports,

        "genomic_summary": {
            gene: {
                "diplotype":     data["diplotype"],
                "phenotype":     data["phenotype"],
                "activity_score": data["total_activity_score"],
                "detected_rsids": data["detected_rsids"],
                "reasoning":     data["phenotype_reasoning"],
            }
            for gene, data in multi_result["gene_profiles"].items()
        },

        "quality_metrics": {
            "vcf_parsing_success":  parse_success,
            "parse_error":          parse_error,
            "total_variants_found": len(variants),
            "annotated_variants":   sum(
                len(data["detected_rsids"])
                for data in multi_result["gene_profiles"].values()
            ),
            "skipped_drugs":        multi_result.get("skipped_drugs", []),
        },
    }

    # Cleanup
    try:
        os.remove(file_location)
    except Exception:
        pass

    return response
