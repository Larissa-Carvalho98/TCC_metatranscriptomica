"""
Módulo para anotação de variantes usando SnpEff.
"""

import subprocess
import pandas as pd
from pathlib import Path
from typing import Optional


def annotate_variants_snpeff(
    input_vcf: str,
    output_vcf: str,
    genome: str = "sarsCov2"
) -> str:
    """
    Anota variantes usando SnpEff.
    
    Args:
        input_vcf: Arquivo VCF de entrada (com variantes)
        output_vcf: Arquivo VCF anotado de saída
        genome: Nome do genoma de referência (padrão: sarsCov2)
    
    Returns:
        str: Caminho do arquivo VCF anotado
    """
    if not Path(input_vcf).exists():
        raise FileNotFoundError(f"Arquivo VCF não encontrado: {input_vcf}")
    
    Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "snpEff",
        "-v", genome,
        input_vcf
    ]
    
    print(f"Anotando variantes com SnpEff...")
    result = subprocess.run(cmd, stdout=open(output_vcf, 'w'), stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"SnpEff falhou: {result.stderr}")
    
    print(f"Variantes anotadas: {output_vcf}")
    return output_vcf


def parse_annotated_vcf(vcf_file: str) -> pd.DataFrame:
    """
    Parseia arquivo VCF anotado e extrai informações em DataFrame.
    
    Args:
        vcf_file: Caminho para arquivo VCF anotado
    
    Returns:
        pd.DataFrame: DataFrame com variantes e anotações
    """
    if not Path(vcf_file).exists():
        raise FileNotFoundError(f"Arquivo VCF não encontrado: {vcf_file}")
    
    variants = []
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            qual = fields[5] if fields[5] != '.' else None
            info = fields[7]
            
            # Extrair informações do campo INFO
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
                else:
                    info_dict[item] = True
            
            # Extrair anotações do SnpEff (campo ANN)
            variant = {
                'CHROM': chrom,
                'POS': pos,
                'REF': ref,
                'ALT': alt,
                'QUAL': qual,
                'DP': info_dict.get('DP', None),
                'AF': info_dict.get('AF', None)
            }
            
            if 'ANN' in info_dict:
                ann_fields = info_dict['ANN'].split('|')
                # Formato SnpEff: Allele | Annotation | Annotation_Impact | Gene_Name | ...
                variant['GENE'] = ann_fields[3] if len(ann_fields) > 3 else None
                variant['EFFECT'] = ann_fields[1] if len(ann_fields) > 1 else None
                variant['IMPACT'] = ann_fields[2] if len(ann_fields) > 2 else None
            else:
                variant['GENE'] = None
                variant['EFFECT'] = None
                variant['IMPACT'] = None
            
            variants.append(variant)
    
    df = pd.DataFrame(variants)
    return df


def extract_spike_mutations(variants_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filtra mutações no gene Spike (S).
    
    Args:
        variants_df: DataFrame com variantes anotadas
    
    Returns:
        pd.DataFrame: DataFrame com apenas mutações no gene S
    """
    spike_mutations = variants_df[variants_df['GENE'] == 'S'].copy()
    return spike_mutations
