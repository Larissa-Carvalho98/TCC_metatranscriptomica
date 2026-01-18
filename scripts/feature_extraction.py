"""
Módulo para extração de features a partir de dados genômicos.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Optional


def extract_tpm_from_kraken2(kraken_report: str) -> Dict[str, float]:
    """
    Extrai valores TPM (Transcripts Per Million) do relatório Kraken2.
    
    Args:
        kraken_report: Caminho para relatório Kraken2
    
    Returns:
        dict: Dicionário com nome da espécie -> TPM
    """
    if not Path(kraken_report).exists():
        raise FileNotFoundError(f"Relatório não encontrado: {kraken_report}")
    
    df = pd.read_csv(
        kraken_report,
        sep='\t',
        names=['percent', 'reads', 'tax_reads', 'rank', 'taxid', 'name']
    )
    
    # Filtrar apenas espécies (rank == 'S')
    species = df[df['rank'] == 'S'].copy()
    
    # Calcular TPM
    total_reads = species['reads'].sum()
    if total_reads > 0:
        species['TPM'] = (species['reads'] / total_reads) * 1_000_000
    else:
        species['TPM'] = 0.0
    
    # Criar dicionário nome -> TPM
    tpm_dict = dict(zip(species['name'], species['TPM']))
    
    return tpm_dict


def extract_tpm_from_bracken(bracken_report: str) -> Dict[str, float]:
    """
    Extrai valores TPM do relatório Bracken (já normalizado).
    
    Args:
        bracken_report: Caminho para relatório Bracken
    
    Returns:
        dict: Dicionário com nome da espécie -> TPM
    """
    if not Path(bracken_report).exists():
        raise FileNotFoundError(f"Relatório não encontrado: {bracken_report}")
    
    df = pd.read_csv(bracken_report, sep='\t')
    
    # Bracken já fornece abundância estimada
    # Usar 'new_est_reads' ou 'fraction_total_reads' dependendo do formato
    if 'new_est_reads' in df.columns:
        total_reads = df['new_est_reads'].sum()
        if total_reads > 0:
            df['TPM'] = (df['new_est_reads'] / total_reads) * 1_000_000
        else:
            df['TPM'] = 0.0
    elif 'fraction_total_reads' in df.columns:
        df['TPM'] = df['fraction_total_reads'] * 1_000_000
    else:
        raise ValueError("Formato de relatório Bracken não reconhecido")
    
    # Criar dicionário nome -> TPM
    tpm_dict = dict(zip(df['name'], df['TPM']))
    
    return tpm_dict


def extract_coverage_stats(coverage_file: str) -> Dict[str, float]:
    """
    Extrai estatísticas de cobertura a partir de arquivo samtools depth.
    
    Args:
        coverage_file: Caminho para arquivo de cobertura (formato: chr pos depth)
    
    Returns:
        dict: Dicionário com estatísticas de cobertura
    """
    if not Path(coverage_file).exists():
        raise FileNotFoundError(f"Arquivo de cobertura não encontrado: {coverage_file}")
    
    df = pd.read_csv(
        coverage_file,
        sep='\t',
        names=['chr', 'pos', 'depth']
    )
    
    stats = {
        'mean_depth': df['depth'].mean(),
        'median_depth': df['depth'].median(),
        'max_depth': df['depth'].max(),
        'min_depth': df['depth'].min(),
        'std_depth': df['depth'].std(),
        'coverage_percentage': (df['depth'] > 0).sum() / len(df) * 100
    }
    
    return stats


def extract_variant_features(variants_table: str) -> Dict[str, int]:
    """
    Extrai features relacionadas a variantes.
    
    Args:
        variants_table: Caminho para tabela CSV de variantes anotadas
    
    Returns:
        dict: Dicionário com features de variantes
    """
    if not Path(variants_table).exists():
        raise FileNotFoundError(f"Tabela de variantes não encontrada: {variants_table}")
    
    df = pd.read_csv(variants_table)
    
    features = {
        'total_variants': len(df),
        'spike_mutations': len(df[df.get('GENE', '') == 'S']) if 'GENE' in df.columns else 0,
        'non_synonymous': len(df[df.get('EFFECT', '') != 'synonymous_variant']) if 'EFFECT' in df.columns else 0,
        'high_impact': len(df[df.get('IMPACT', '') == 'HIGH']) if 'IMPACT' in df.columns else 0
    }
    
    return features
