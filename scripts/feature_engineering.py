"""
Módulo para engenharia de features avançadas.
"""

import pandas as pd
import numpy as np
from scipy.stats import entropy
from pathlib import Path
from typing import Dict, Optional


def calculate_shannon_index(tpm_dict: Dict[str, float]) -> float:
    """
    Calcula índice de diversidade Shannon a partir de abundâncias TPM.
    """
    if not tpm_dict:
        return 0.0

    abundances = np.array(list(tpm_dict.values()))
    abundances = abundances[abundances > 0]

    if len(abundances) == 0:
        return 0.0

    proportions = abundances / abundances.sum()
    return entropy(proportions, base=2)


def calculate_virus_bacteria_ratio(kraken_report: str) -> float:
    """
    Calcula razão vírus/bactérias a partir de relatório Kraken2.
    """
    if not Path(kraken_report).exists():
        raise FileNotFoundError(f"Relatório não encontrado: {kraken_report}")

    df = pd.read_csv(
        kraken_report,
        sep="\t",
        names=["percent", "reads", "tax_reads", "rank", "taxid", "name"]
    )

    virus_mask = df["name"].str.contains("virus", case=False, na=False) | (df["rank"] == "V")
    virus_reads = df.loc[virus_mask, "reads"].sum()

    bacteria_mask = (df["rank"] == "S") & ~virus_mask
    bacteria_reads = df.loc[bacteria_mask, "reads"].sum()

    if bacteria_reads > 0:
        return virus_reads / bacteria_reads
    return float("inf") if virus_reads > 0 else 0.0


def detect_key_pathogens(kraken_report: str) -> Dict[str, bool]:
    """
    Detecta presença de patógenos-chave.
    """
    if not Path(kraken_report).exists():
        raise FileNotFoundError(f"Relatório não encontrado: {kraken_report}")

    df = pd.read_csv(
        kraken_report,
        sep="\t",
        names=["percent", "reads", "tax_reads", "rank", "taxid", "name"]
    )

    all_names = " ".join(df["name"].astype(str))

    return {
        "SARS-CoV-2": (
            "Severe acute respiratory syndrome coronavirus 2" in all_names
            or "SARS-CoV-2" in all_names
        ),
        "Influenza": "Influenza" in all_names,
        "RSV": "Respiratory syncytial" in all_names,
        "Streptococcus pneumoniae": "Streptococcus pneumoniae" in all_names,
        "Haemophilus influenzae": "Haemophilus influenzae" in all_names,
        "Moraxella catarrhalis": "Moraxella catarrhalis" in all_names,
    }


def create_feature_vector(
    tpm_dict: Dict[str, float],
    training_columns: list,
    advanced_features: Optional[Dict] = None
) -> pd.DataFrame:
    """
    Cria vetor de features alinhado com matriz de treinamento
    (implementação otimizada, sem fragmentação do DataFrame).
    """
    feature_data = {}

    for organism in training_columns:
        # Match exato
        tpm_value = tpm_dict.get(organism, 0.0)

        # Match parcial (case-insensitive), se necessário
        if tpm_value == 0.0:
            for org_name, tpm in tpm_dict.items():
                if (
                    organism.lower() in org_name.lower()
                    or org_name.lower() in organism.lower()
                ):
                    tpm_value = tpm
                    break

        feature_data[organism] = tpm_value

    # Criar DataFrame de uma vez (1 linha)
    feature_vector = pd.DataFrame([feature_data])

    # Adicionar features avançadas
    if advanced_features:
        for feature_name, feature_value in advanced_features.items():
            feature_vector[feature_name] = feature_value

    return feature_vector


def apply_log_transform(feature_vector: pd.DataFrame) -> pd.DataFrame:
    """
    Aplica transformação logarítmica (log1p) ao vetor de features.
    """
    return np.log1p(feature_vector)
