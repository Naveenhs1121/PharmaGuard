"""
PharmaGuard — Explanation Engine
Generates clinical explanation text from risk engine output.
Designed for LLM integration; currently uses structured rule-based templates.
"""

import logging
from typing import Dict, List, Any

logger = logging.getLogger("PharmaGuard.ExplanationEngine")

# Human-friendly phenotype descriptions
PHENOTYPE_DESC = {
    "Poor Metabolizer":         "cannot efficiently metabolize this drug due to reduced enzyme activity",
    "Intermediate Metabolizer": "has partially reduced enzyme activity affecting drug metabolism",
    "Normal Metabolizer":       "metabolizes this drug at a standard rate",
    "Rapid Metabolizer":        "metabolizes this drug faster than average",
    "Ultrarapid Metabolizer":   "metabolizes this drug at an exceptionally high rate",
    "Indeterminate":            "has an uncertain metabolizer status based on available data",
}

SEVERITY_CONTEXT = {
    "none":     "No clinically significant pharmacogenomic interaction is expected.",
    "low":      "A minor interaction is noted; routine monitoring is advised.",
    "moderate": "A clinically significant interaction exists; dose adjustment is recommended.",
    "high":     "A serious pharmacogenomic interaction is identified; immediate clinical action is required.",
    "unknown":  "The interaction risk is uncertain due to insufficient genomic data.",
}


def generate_explanation(
    variants: List[Dict],
    drug: str,
    risk: Dict[str, Any],
    gene_profile: Dict = None,
    clinical_action: str = "",
) -> str:
    """
    Generate a clinical explanation string.

    Args:
        variants        : list of VCF variant dicts
        drug            : drug name
        risk            : dict with keys label, severity, confidence
        gene_profile    : optional gene analysis result dict
        clinical_action : CPIC-recommended action string

    Returns:
        A human-readable clinical explanation paragraph.
    """
    label       = risk.get("label", "Unknown")
    severity    = risk.get("severity", "unknown")
    confidence  = risk.get("confidence", 0.0)

    # ── No variants detected ──────────────────────────────────────────────
    if not variants:
        return (
            f"No pharmacogenomic variants relevant to {drug.capitalize()} were detected "
            f"in the uploaded VCF file. In the absence of known risk variants, "
            f"standard {drug.capitalize()} dosing is generally appropriate. "
            f"Clinical judgment should guide prescribing decisions."
        )

    # ── Build gene/phenotype context ──────────────────────────────────────
    gene = "the relevant gene"
    phenotype = "Indeterminate"
    diplotype = "unknown"
    rsids = []
    activity_score = None

    if gene_profile:
        gene        = gene_profile.get("gene", gene)
        phenotype   = gene_profile.get("phenotype", phenotype)
        diplotype   = gene_profile.get("diplotype", diplotype)
        rsids       = gene_profile.get("detected_rsids", [])
        activity_score = gene_profile.get("total_activity_score")

    phenotype_desc = PHENOTYPE_DESC.get(phenotype, "has an atypical metabolizer status")
    severity_ctx   = SEVERITY_CONTEXT.get(severity, SEVERITY_CONTEXT["unknown"])

    # ── Variant summary ───────────────────────────────────────────────────
    if rsids:
        variant_text = (
            f"The following variants were identified: {', '.join(rsids)}. "
            f"These map to a {diplotype} diplotype"
        )
        if activity_score is not None:
            variant_text += f" with an activity score of {activity_score:.2f}"
        variant_text += "."
    else:
        variant_text = (
            f"No database-annotated variants were found in {gene}, "
            f"so a reference diplotype (*1/*1) was assumed."
        )

    # ── Build explanation ─────────────────────────────────────────────────
    explanation_parts = [
        f"Pharmacogenomic analysis of {drug.capitalize()} indicates a risk "
        f"classification of '{label}' (severity: {severity}, confidence: {confidence:.0%}).",

        f"The patient {phenotype_desc} ({phenotype}). {variant_text}",

        severity_ctx,
    ]

    if clinical_action:
        explanation_parts.append(f"Recommended action: {clinical_action}")

    explanation_parts.append(
        f"This assessment is based on CPIC (Clinical Pharmacogenomics Implementation "
        f"Consortium) guidelines and should be interpreted in conjunction with the "
        f"patient's complete clinical picture."
    )

    return " ".join(explanation_parts)
