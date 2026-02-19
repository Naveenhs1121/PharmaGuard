"""
PharmaGuard Risk Engine
=======================
CPIC-aligned pharmacogenomic risk prediction engine.

Pipeline:
  VCF variants → rsID annotation → diplotype inference
  → phenotype classification → drug-phenotype rule lookup
  → confidence scoring → structured clinical output

Supported Genes : CYP2D6, CYP2C19, CYP2C9, SLCO1B1, TPMT, DPYD
Supported Drugs : Codeine, Warfarin, Clopidogrel, Simvastatin,
                  Azathioprine, Fluorouracil
"""

import logging
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Any

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("PharmaGuard.RiskEngine")

# ---------------------------------------------------------------------------
# Constants — Phenotype codes
# ---------------------------------------------------------------------------
PM  = "Poor Metabolizer"
IM  = "Intermediate Metabolizer"
NM  = "Normal Metabolizer"
RM  = "Rapid Metabolizer"
URM = "Ultrarapid Metabolizer"
INDETERMINATE = "Indeterminate"

# Metabolizer activity scores (used for diplotype → phenotype mapping)
ACTIVITY_SCORE = {
    "LOF": 0.0,   # Loss-of-function allele
    "DEF": 0.5,   # Decreased-function allele
    "NF":  1.0,   # Normal-function allele
    "INF": 2.0,   # Increased-function allele (gene duplication / gain)
}

# ---------------------------------------------------------------------------
# rsID → Allele Annotation Database
# Maps known pharmacogenomic rsIDs to:
#   gene, star_allele, function, activity_score, evidence_strength
# Evidence strength: "high" | "moderate" | "low"
# ---------------------------------------------------------------------------
RSID_DATABASE: Dict[str, Dict] = {
    # ── CYP2D6 ──────────────────────────────────────────────────────────────
    "rs3892097":  {"gene": "CYP2D6",  "star": "*4",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs35742686": {"gene": "CYP2D6",  "star": "*3",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs5030655":  {"gene": "CYP2D6",  "star": "*6",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs28371725": {"gene": "CYP2D6",  "star": "*41", "function": "DEF", "activity": 0.5, "evidence": "high"},
    "rs1065852":  {"gene": "CYP2D6",  "star": "*10", "function": "DEF", "activity": 0.5, "evidence": "high"},
    "rs16947":    {"gene": "CYP2D6",  "star": "*2",  "function": "NF",  "activity": 1.0, "evidence": "high"},
    "rs1135840":  {"gene": "CYP2D6",  "star": "*2",  "function": "NF",  "activity": 1.0, "evidence": "moderate"},
    "rs5030865":  {"gene": "CYP2D6",  "star": "*8",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},

    # ── CYP2C19 ─────────────────────────────────────────────────────────────
    "rs4244285":  {"gene": "CYP2C19", "star": "*2",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs4986893":  {"gene": "CYP2C19", "star": "*3",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs28399504": {"gene": "CYP2C19", "star": "*4",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},
    "rs56337013": {"gene": "CYP2C19", "star": "*5",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},
    "rs12248560": {"gene": "CYP2C19", "star": "*17", "function": "INF", "activity": 1.5, "evidence": "high"},
    "rs41291556": {"gene": "CYP2C19", "star": "*17", "function": "INF", "activity": 1.5, "evidence": "moderate"},

    # ── CYP2C9 ──────────────────────────────────────────────────────────────
    "rs1799853":  {"gene": "CYP2C9",  "star": "*2",  "function": "DEF", "activity": 0.5, "evidence": "high"},
    "rs1057910":  {"gene": "CYP2C9",  "star": "*3",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs28371686": {"gene": "CYP2C9",  "star": "*5",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},
    "rs9332131":  {"gene": "CYP2C9",  "star": "*6",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},
    "rs7900194":  {"gene": "CYP2C9",  "star": "*8",  "function": "DEF", "activity": 0.5, "evidence": "moderate"},

    # ── SLCO1B1 ─────────────────────────────────────────────────────────────
    "rs4149056":  {"gene": "SLCO1B1", "star": "*5",  "function": "DEF", "activity": 0.5, "evidence": "high"},
    "rs2306283":  {"gene": "SLCO1B1", "star": "*1B", "function": "NF",  "activity": 1.0, "evidence": "moderate"},
    "rs11045819": {"gene": "SLCO1B1", "star": "*15", "function": "DEF", "activity": 0.5, "evidence": "moderate"},
    "rs4363657":  {"gene": "SLCO1B1", "star": "*5",  "function": "DEF", "activity": 0.5, "evidence": "high"},

    # ── TPMT ────────────────────────────────────────────────────────────────
    "rs1142345":  {"gene": "TPMT",    "star": "*3C", "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs1800460":  {"gene": "TPMT",    "star": "*3B", "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs1800462":  {"gene": "TPMT",    "star": "*2",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs1800584":  {"gene": "TPMT",    "star": "*4",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},

    # ── DPYD ────────────────────────────────────────────────────────────────
    "rs3918290":  {"gene": "DPYD",    "star": "*2A", "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs55886062": {"gene": "DPYD",    "star": "*13", "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs67376798": {"gene": "DPYD",    "star": "HapB3","function": "DEF","activity": 0.5, "evidence": "high"},
    "rs75017182": {"gene": "DPYD",    "star": "HapB3","function": "DEF","activity": 0.5, "evidence": "moderate"},
}

# ---------------------------------------------------------------------------
# Gene → Phenotype rules (activity-score based, CPIC aligned)
# Each rule: (lower_bound_inclusive, upper_bound_exclusive, phenotype)
# ---------------------------------------------------------------------------
GENE_PHENOTYPE_RULES: Dict[str, List] = {
    "CYP2D6": [
        (0.0, 0.01, PM),
        (0.01, 1.25, IM),
        (1.25, 2.25, NM),
        (2.25, 3.0,  RM),
        (3.0,  99.0, URM),
    ],
    "CYP2C19": [
        (0.0, 0.01, PM),
        (0.01, 1.25, IM),
        (1.25, 2.25, NM),
        (2.25, 3.0,  RM),
        (3.0,  99.0, URM),
    ],
    "CYP2C9": [
        (0.0, 0.01, PM),
        (0.01, 1.25, IM),
        (1.25, 99.0, NM),
    ],
    "SLCO1B1": [
        # For SLCO1B1: activity score maps to transporter function
        (0.0, 0.3,  PM),    # Poor function → high myopathy risk
        (0.3, 0.9,  IM),    # Decreased function
        (0.9, 99.0, NM),    # Normal function
    ],
    "TPMT": [
        (0.0, 0.01, PM),
        (0.01, 1.25, IM),
        (1.25, 99.0, NM),
    ],
    "DPYD": [
        (0.0, 0.01, PM),
        (0.01, 0.9,  IM),
        (0.9, 99.0, NM),
    ],
}

# ---------------------------------------------------------------------------
# CPIC Drug–Phenotype Rules
# Structure:
#   drug → gene → phenotype → {label, severity, confidence_base,
#                               clinical_action, cpic_guideline}
#
# Risk labels (hackathon spec): Safe | Adjust Dosage | Toxic | Ineffective | Unknown
# ---------------------------------------------------------------------------
DRUG_RULES: Dict[str, Dict[str, Dict[str, Dict]]] = {
    "codeine": {
        "CYP2D6": {
            PM: {
                "label": "Ineffective",
                "severity": "high",
                "confidence_base": 0.95,
                "clinical_action": (
                    "Avoid codeine. CYP2D6 Poor Metabolizers cannot convert codeine "
                    "to morphine, resulting in no analgesia. Use a non-opioid analgesic "
                    "or a non-CYP2D6-metabolised opioid (e.g., morphine, oxycodone)."
                ),
                "cpic_guideline": "CPIC Codeine Guideline (2021) — PMID 22585173",
            },
            IM: {
                "label": "Adjust Dosage",
                "severity": "moderate",
                "confidence_base": 0.80,
                "clinical_action": (
                    "Use codeine with caution. Intermediate Metabolizers produce "
                    "reduced morphine levels. Consider a lower dose or an alternative analgesic."
                ),
                "cpic_guideline": "CPIC Codeine Guideline (2021) — PMID 22585173",
            },
            NM: {
                "label": "Safe",
                "severity": "none",
                "confidence_base": 0.90,
                "clinical_action": "Standard codeine dosing as per label.",
                "cpic_guideline": "CPIC Codeine Guideline (2021)",
            },
            RM: {
                "label": "Safe",
                "severity": "low",
                "confidence_base": 0.85,
                "clinical_action": (
                    "Standard dosing. Slightly higher morphine conversion — monitor "
                    "for opioid side effects."
                ),
                "cpic_guideline": "CPIC Codeine Guideline (2021)",
            },
            URM: {
                "label": "Toxic",
                "severity": "high",
                "confidence_base": 0.95,
                "clinical_action": (
                    "CONTRAINDICATED. Ultrarapid Metabolizers convert codeine to morphine "
                    "at extremely high rates, risking life-threatening opioid toxicity. "
                    "Use an alternative non-opioid or tramadol with caution."
                ),
                "cpic_guideline": "CPIC Codeine Guideline (2021) — PMID 22585173",
            },
        }
    },
    "warfarin": {
        "CYP2C9": {
            PM: {
                "label": "Adjust Dosage",
                "severity": "high",
                "confidence_base": 0.92,
                "clinical_action": (
                    "Significant warfarin dose reduction required (up to 50–80% lower). "
                    "CYP2C9 Poor Metabolizers have severely reduced warfarin clearance. "
                    "Start low, titrate slowly, and monitor INR closely."
                ),
                "cpic_guideline": "CPIC Warfarin Guideline — PMID 21900891",
            },
            IM: {
                "label": "Adjust Dosage",
                "severity": "moderate",
                "confidence_base": 0.85,
                "clinical_action": (
                    "Reduce starting warfarin dose by 20–40%. Monitor INR frequently "
                    "during initiation. CYP2C9 Intermediate Metabolizers clear warfarin "
                    "more slowly than normal."
                ),
                "cpic_guideline": "CPIC Warfarin Guideline — PMID 21900891",
            },
            NM: {
                "label": "Safe",
                "severity": "none",
                "confidence_base": 0.90,
                "clinical_action": "Standard warfarin dosing per clinical guidelines.",
                "cpic_guideline": "CPIC Warfarin Guideline",
            },
        }
    },
    "clopidogrel": {
        "CYP2C19": {
            PM: {
                "label": "Ineffective",
                "severity": "high",
                "confidence_base": 0.95,
                "clinical_action": (
                    "Avoid clopidogrel. CYP2C19 Poor Metabolizers cannot activate "
                    "clopidogrel (prodrug), resulting in minimal antiplatelet effect and "
                    "increased risk of cardiovascular events. Use prasugrel or ticagrelor."
                ),
                "cpic_guideline": "CPIC Clopidogrel Guideline — PMID 23698643",
            },
            IM: {
                "label": "Adjust Dosage",
                "severity": "moderate",
                "confidence_base": 0.80,
                "clinical_action": (
                    "Consider an alternative antiplatelet agent (prasugrel or ticagrelor). "
                    "If clopidogrel is used, higher doses may be needed — assess bleeding risk."
                ),
                "cpic_guideline": "CPIC Clopidogrel Guideline — PMID 23698643",
            },
            NM: {
                "label": "Safe",
                "severity": "none",
                "confidence_base": 0.90,
                "clinical_action": "Standard clopidogrel dosing as per label.",
                "cpic_guideline": "CPIC Clopidogrel Guideline",
            },
            RM: {
                "label": "Safe",
                "severity": "none",
                "confidence_base": 0.88,
                "clinical_action": "Standard clopidogrel dosing. Slightly enhanced activation — monitor.",
                "cpic_guideline": "CPIC Clopidogrel Guideline",
            },
            URM: {
                "label": "Safe",
                "severity": "low",
                "confidence_base": 0.82,
                "clinical_action": (
                    "Standard dosing. Ultrarapid Metabolizers may have enhanced "
                    "antiplatelet effect — monitor for bleeding."
                ),
                "cpic_guideline": "CPIC Clopidogrel Guideline",
            },
        }
    },
    "simvastatin": {
        "SLCO1B1": {
            PM: {
                "label": "Toxic",
                "severity": "high",
                "confidence_base": 0.92,
                "clinical_action": (
                    "High risk of simvastatin-induced myopathy (including rhabdomyolysis). "
                    "SLCO1B1 Poor Function severely impairs hepatic uptake, increasing plasma "
                    "simvastatin levels. Switch to pravastatin or rosuvastatin (less SLCO1B1-dependent)."
                ),
                "cpic_guideline": "CPIC Simvastatin Guideline — PMID 22617227",
            },
            IM: {
                "label": "Adjust Dosage",
                "severity": "moderate",
                "confidence_base": 0.82,
                "clinical_action": (
                    "Prescribe a lower simvastatin dose (≤20 mg/day) or switch to an "
                    "alternative statin (pravastatin, fluvastatin, or rosuvastatin). "
                    "Monitor for muscle pain or weakness."
                ),
                "cpic_guideline": "CPIC Simvastatin Guideline — PMID 22617227",
            },
            NM: {
                "label": "Safe",
                "severity": "none",
                "confidence_base": 0.90,
                "clinical_action": "Standard simvastatin dosing as per label.",
                "cpic_guideline": "CPIC Simvastatin Guideline",
            },
        }
    },
    "azathioprine": {
        "TPMT": {
            PM: {
                "label": "Toxic",
                "severity": "high",
                "confidence_base": 0.97,
                "clinical_action": (
                    "TPMT Poor Metabolizers accumulate toxic thioguanine nucleotides, "
                    "causing life-threatening myelosuppression. Reduce dose by 90% or "
                    "switch to an alternative immunosuppressant. Increase monitoring frequency."
                ),
                "cpic_guideline": "CPIC Thiopurines Guideline — PMID 21270794",
            },
            IM: {
                "label": "Adjust Dosage",
                "severity": "moderate",
                "confidence_base": 0.88,
                "clinical_action": (
                    "Reduce azathioprine starting dose by 30–70% of standard dose. "
                    "Monitor CBC frequently during therapy. Consider testing NUDT15 as well."
                ),
                "cpic_guideline": "CPIC Thiopurines Guideline — PMID 21270794",
            },
            NM: {
                "label": "Safe",
                "severity": "none",
                "confidence_base": 0.90,
                "clinical_action": "Standard azathioprine dosing as per label. Routine CBC monitoring.",
                "cpic_guideline": "CPIC Thiopurines Guideline",
            },
        }
    },
    "fluorouracil": {
        "DPYD": {
            PM: {
                "label": "Toxic",
                "severity": "high",
                "confidence_base": 0.97,
                "clinical_action": (
                    "DPYD Poor Metabolizers cannot adequately catabolize fluorouracil, "
                    "leading to severe, life-threatening toxicity (mucositis, neutropenia, "
                    "neurotoxicity). Do NOT use standard doses. Consider alternative "
                    "fluoropyrimidine therapy or reduce dose by ≥50% with close monitoring."
                ),
                "cpic_guideline": "CPIC Fluoropyrimidines Guideline — PMID 23988873",
            },
            IM: {
                "label": "Adjust Dosage",
                "severity": "moderate",
                "confidence_base": 0.88,
                "clinical_action": (
                    "Reduce fluorouracil starting dose by 50%. Increase toxicity monitoring. "
                    "DPYD Intermediate Metabolizers have reduced drug clearance — standard doses "
                    "may cause significant adverse events."
                ),
                "cpic_guideline": "CPIC Fluoropyrimidines Guideline — PMID 23988873",
            },
            NM: {
                "label": "Safe",
                "severity": "none",
                "confidence_base": 0.90,
                "clinical_action": "Standard fluorouracil dosing as per protocol.",
                "cpic_guideline": "CPIC Fluoropyrimidines Guideline",
            },
        }
    },
}

# ---------------------------------------------------------------------------
# Data classes for structured output
# ---------------------------------------------------------------------------
@dataclass
class VariantAnnotation:
    rsid: str
    gene: str
    star_allele: str
    function: str
    activity_score: float
    evidence_strength: str
    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: List[str] = field(default_factory=list)


@dataclass
class GeneResult:
    gene: str
    detected_rsids: List[str]
    annotated_variants: List[VariantAnnotation]
    total_activity_score: float
    diplotype: str
    phenotype: str
    phenotype_reasoning: str


@dataclass
class DrugRiskResult:
    drug: str
    gene_used: str
    phenotype: str
    risk_label: str
    severity: str
    confidence_score: float
    clinical_action: str
    cpic_guideline: str
    supporting_variants: List[str]
    reasoning: str
    evidence_strength: str


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def _normalise_drug(drug: str) -> str:
    """Lowercase and strip the drug name for rule lookup."""
    return drug.strip().lower().replace("-", "").replace(" ", "")


def _annotate_variant(variant: Dict) -> Optional[VariantAnnotation]:
    """
    Match a VCF variant dict to the rsID database.
    Returns VariantAnnotation or None if not in database.
    """
    rsid_raw = variant.get("rsid")
    if not rsid_raw:
        return None

    rsid_list = rsid_raw if isinstance(rsid_raw, list) else [rsid_raw]

    for rsid in rsid_list:
        if rsid and rsid in RSID_DATABASE:
            db = RSID_DATABASE[rsid]
            return VariantAnnotation(
                rsid=rsid,
                gene=db["gene"],
                star_allele=db["star"],
                function=db["function"],
                activity_score=db["activity"],
                evidence_strength=db["evidence"],
                chrom=str(variant.get("chrom", "")),
                pos=int(variant.get("pos", 0)),
                ref=str(variant.get("ref", "")),
                alt=[str(a) for a in variant.get("alt", [])],
            )
    return None


def _classify_phenotype(gene: str, activity_score: float) -> str:
    """Map a gene's total activity score to its metabolizer phenotype."""
    rules = GENE_PHENOTYPE_RULES.get(gene)
    if rules is None:
        return INDETERMINATE

    for lo, hi, phenotype in rules:
        if lo <= activity_score < hi:
            return phenotype

    # Edge case: score exactly at the max end
    return rules[-1][2]


def _build_diplotype_label(annotated: List[VariantAnnotation], gene: str) -> str:
    """
    Build a human-readable diplotype string from detected star alleles.
    Assumes biallelic diploid (two allele slots).
    """
    stars = sorted(set(a.star_allele for a in annotated))
    if not stars:
        return f"{gene}:*1/*1 (reference assumed)"
    if len(stars) == 1:
        return f"{gene}:{stars[0]}/{stars[0]}"
    return f"{gene}:{stars[0]}/{stars[1]}"


def _evidence_strength_factor(annotated: List[VariantAnnotation]) -> float:
    """
    Returns a multiplier [0.6, 1.0] based on the weakest evidence level
    in the supporting variants, penalising low-evidence calls.
    """
    if not annotated:
        return 0.6
    weights = {"high": 1.0, "moderate": 0.85, "low": 0.65}
    return min(weights.get(a.evidence_strength, 0.65) for a in annotated)


def _phenotype_reasoning(gene: str, phenotype: str, score: float,
                          annotated: List[VariantAnnotation]) -> str:
    """Generate a plain-English phenotype reasoning string."""
    variant_list = ", ".join(
        f"{a.rsid} ({a.star_allele}, {a.function})" for a in annotated
    ) or "no database-annotated variants"
    return (
        f"{gene} activity score = {score:.2f}. "
        f"Detected variants: {variant_list}. "
        f"Phenotype classified as {phenotype} per CPIC activity-score model."
    )


def _risk_reasoning(drug: str, gene: str, phenotype: str,
                    rule: Dict, score: float) -> str:
    """Generate a plain-English risk reasoning string."""
    return (
        f"{drug.capitalize()} risk assessment based on {gene} phenotype "
        f"({phenotype}, activity score {score:.2f}). "
        f"CPIC classification: '{rule['label']}' "
        f"(severity: {rule['severity']})."
    )


def _compute_confidence(
    base: float,
    activity_score: float,
    n_variants: int,
    evidence_factor: float,
    has_variants: bool,
) -> float:
    """
    Compute final confidence score on [0.0, 1.0].

    Factors:
    - base: rule-level CPIC confidence
    - variant_support: max 0.05 bonus for ≥ 2 supporting variants
    - evidence_factor: penalty for low-evidence rsIDs
    - no_variant_penalty: small deduction if phenotype inferred from absence
    """
    if not has_variants:
        # Inferred as NM from absence of known LOF variants — less certain
        return round(min(base * evidence_factor * 0.85, 1.0), 3)

    variant_bonus = min(n_variants * 0.02, 0.05)
    confidence = (base + variant_bonus) * evidence_factor
    return round(min(confidence, 1.0), 3)


# ---------------------------------------------------------------------------
# Core prediction functions
# ---------------------------------------------------------------------------

def analyse_gene(gene: str, variants: List[Dict]) -> GeneResult:
    """
    Analyse all variants for a single gene.

    Returns a GeneResult with diplotype, phenotype, and reasoning.
    """
    logger.info(f"Analysing gene: {gene} with {len(variants)} candidate variants")

    gene_variants = [v for v in variants if v.get("gene") == gene]
    annotated: List[VariantAnnotation] = []

    for v in gene_variants:
        ann = _annotate_variant(v)
        if ann and ann.gene == gene:
            annotated.append(ann)
            logger.debug(f"  {gene}: matched {ann.rsid} → {ann.star_allele} ({ann.function})")
        else:
            rsid = v.get("rsid", "unknown")
            logger.debug(f"  {gene}: rsID {rsid!r} not in database — skipped")

    # Activity score: sum of all detected variant scores
    # Biallelic assumption: if no LOF detected, assume *1/*1 (score = 2.0)
    if annotated:
        activity_score = sum(a.activity_score for a in annotated)
    else:
        activity_score = 2.0   # Reference (*1/*1) assumed

    phenotype = _classify_phenotype(gene, activity_score)
    diplotype = _build_diplotype_label(annotated, gene)
    reasoning = _phenotype_reasoning(gene, phenotype, activity_score, annotated)

    logger.info(f"  {gene}: activity={activity_score:.2f}, phenotype={phenotype}, diplotype={diplotype}")

    return GeneResult(
        gene=gene,
        detected_rsids=[a.rsid for a in annotated],
        annotated_variants=annotated,
        total_activity_score=activity_score,
        diplotype=diplotype,
        phenotype=phenotype,
        phenotype_reasoning=reasoning,
    )


def predict_drug_risk(drug_name: str, gene_results: Dict[str, GeneResult]) -> DrugRiskResult:
    """
    Predict risk for a single drug given pre-computed gene results.

    Returns a DrugRiskResult with risk label, severity, confidence, and
    clinical action.
    """
    drug_key = _normalise_drug(drug_name)
    logger.info(f"Predicting risk for drug: {drug_name!r} (key={drug_key!r})")

    drug_rules = DRUG_RULES.get(drug_key)
    if not drug_rules:
        logger.warning(f"Drug '{drug_name}' not in rule database.")
        return DrugRiskResult(
            drug=drug_name,
            gene_used="N/A",
            phenotype=INDETERMINATE,
            risk_label="Unknown",
            severity="unknown",
            confidence_score=0.0,
            clinical_action=(
                f"No pharmacogenomic data available for '{drug_name}'. "
                "Proceed per standard clinical guidelines."
            ),
            cpic_guideline="N/A",
            supporting_variants=[],
            reasoning=f"Drug '{drug_name}' is not in the current CPIC rule set.",
            evidence_strength="none",
        )

    # Try each gene defined for this drug (priority order matters)
    for gene, phenotype_rules in drug_rules.items():
        gene_result = gene_results.get(gene)
        if gene_result is None:
            logger.debug(f"  No gene result available for {gene} — skipping rule")
            continue

        phenotype = gene_result.phenotype
        if phenotype == INDETERMINATE:
            logger.debug(f"  {gene} phenotype is Indeterminate — skipping")
            continue

        rule = phenotype_rules.get(phenotype)
        if rule is None:
            logger.debug(f"  No rule for {gene}/{phenotype} combo — defaulting to Safe")
            rule = {
                "label": "Safe",
                "severity": "none",
                "confidence_base": 0.70,
                "clinical_action": (
                    f"No specific {drug_name} recommendation for {phenotype} "
                    f"{gene} phenotype. Use standard dosing with caution."
                ),
                "cpic_guideline": "No specific CPIC guidance",
            }

        has_variants = len(gene_result.annotated_variants) > 0
        ev_factor = _evidence_strength_factor(gene_result.annotated_variants)
        confidence = _compute_confidence(
            base=rule["confidence_base"],
            activity_score=gene_result.total_activity_score,
            n_variants=len(gene_result.annotated_variants),
            evidence_factor=ev_factor,
            has_variants=has_variants,
        )

        ev_strength = (
            min((a.evidence_strength for a in gene_result.annotated_variants),
                key=lambda x: {"high": 0, "moderate": 1, "low": 2}.get(x, 3))
            if has_variants else "inferred"
        )

        reasoning = _risk_reasoning(
            drug_name, gene, phenotype, rule,
            gene_result.total_activity_score,
        )

        logger.info(
            f"  {drug_name}: gene={gene}, phenotype={phenotype}, "
            f"label={rule['label']}, confidence={confidence:.3f}"
        )

        return DrugRiskResult(
            drug=drug_name,
            gene_used=gene,
            phenotype=phenotype,
            risk_label=rule["label"],
            severity=rule["severity"],
            confidence_score=confidence,
            clinical_action=rule["clinical_action"],
            cpic_guideline=rule["cpic_guideline"],
            supporting_variants=gene_result.detected_rsids,
            reasoning=reasoning,
            evidence_strength=ev_strength,
        )

    # Reached here: no applicable gene result found
    logger.warning(f"No applicable gene/phenotype found for {drug_name}")
    return DrugRiskResult(
        drug=drug_name,
        gene_used="N/A",
        phenotype=INDETERMINATE,
        risk_label="Unknown",
        severity="unknown",
        confidence_score=0.40,
        clinical_action=(
            f"Insufficient genomic data to assess {drug_name} risk. "
            "No variants detected in relevant genes. Proceed with standard care."
        ),
        cpic_guideline="N/A",
        supporting_variants=[],
        reasoning=(
            f"Relevant genes for {drug_name} were not detected in the VCF. "
            "Risk set to Unknown."
        ),
        evidence_strength="none",
    )


# ---------------------------------------------------------------------------
# Public API — main entry points
# ---------------------------------------------------------------------------

TARGET_GENES = ["CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"]
SUPPORTED_DRUGS = [
    "codeine", "warfarin", "clopidogrel",
    "simvastatin", "azathioprine", "fluorouracil",
]


def predict_risk(variants: List[Dict], drug: str) -> Dict[str, Any]:
    """
    Legacy single-drug API (backwards-compatible with main.py).

    Args:
        variants : list of variant dicts from vcf_parser.extract_variants()
        drug     : drug name string

    Returns:
        dict with keys: label, severity, confidence (+ full result nested)
    """
    drug_key = _normalise_drug(drug)
    if drug_key not in SUPPORTED_DRUGS:
        logger.warning(f"Drug '{drug}' not supported — returning Unknown")
        return {
            "label": "Unknown",
            "severity": "unknown",
            "confidence": 0.0,
            "full_result": None,
        }

    gene_results = {
        gene: analyse_gene(gene, variants) for gene in TARGET_GENES
    }
    result = predict_drug_risk(drug, gene_results)
    return {
        "label": result.risk_label,
        "severity": result.severity,
        "confidence": result.confidence_score,
        "full_result": asdict(result),
    }


def predict_multi_drug(variants: List[Dict], drugs: List[str]) -> Dict[str, Any]:
    """
    Multi-drug prediction API.

    Args:
        variants : list of variant dicts from vcf_parser.extract_variants()
        drugs    : list of drug name strings

    Returns:
        dict with:
          gene_profiles  — per-gene analysis (diplotype, phenotype, reasoning)
          drug_results   — per-drug risk prediction
    """
    if not variants and not drugs:
        logger.warning("predict_multi_drug called with empty variants and drugs")

    # Validate inputs
    validated_drugs: List[str] = []
    skipped_drugs: List[str] = []
    for d in drugs:
        key = _normalise_drug(d)
        if key in SUPPORTED_DRUGS:
            validated_drugs.append(d)
        else:
            logger.warning(f"Drug '{d}' is not supported — skipped")
            skipped_drugs.append(d)

    # Step 1: Analyse all genes once
    gene_results: Dict[str, GeneResult] = {}
    for gene in TARGET_GENES:
        gene_results[gene] = analyse_gene(gene, variants)

    # Step 2: Predict risk for each validated drug
    drug_results: List[Dict] = []
    for drug in validated_drugs:
        drug_risk = predict_drug_risk(drug, gene_results)
        drug_results.append(asdict(drug_risk))

    # Step 3: Add Unknown entries for unsupported drugs
    for drug in skipped_drugs:
        drug_results.append({
            "drug": drug,
            "gene_used": "N/A",
            "phenotype": INDETERMINATE,
            "risk_label": "Unknown",
            "severity": "unknown",
            "confidence_score": 0.0,
            "clinical_action": f"'{drug}' is not in the supported drug list.",
            "cpic_guideline": "N/A",
            "supporting_variants": [],
            "reasoning": f"Drug '{drug}' is not supported by this version of PharmaGuard.",
            "evidence_strength": "none",
        })

    return {
        "gene_profiles": {
            gene: asdict(gr) for gene, gr in gene_results.items()
        },
        "drug_results": drug_results,
        "skipped_drugs": skipped_drugs,
    }
