"""
Módulo para determinação de linhagem do SARS-CoV-2 usando Pangolin.
"""

import subprocess
import pandas as pd
from pathlib import Path
from typing import Optional


def determine_lineage_pangolin(
    input_file: str,
    output_csv: str,
    input_type: str = "bam"  # "bam", "vcf", ou "fasta"
) -> pd.DataFrame:
    """
    Determina linhagem do SARS-CoV-2 usando Pangolin.
    
    Args:
        input_file: Arquivo de entrada (BAM, VCF ou FASTA consenso)
        output_csv: Arquivo CSV de saída com resultados
        input_type: Tipo de arquivo de entrada ("bam", "vcf", "fasta")
    
    Returns:
        pd.DataFrame: DataFrame com resultados do Pangolin
    """
    if not Path(input_file).exists():
        raise FileNotFoundError(f"Arquivo de entrada não encontrado: {input_file}")
    
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    
    # Pangolin tem opções diferentes dependendo do tipo de entrada
    # Para BAM ou VCF, pode ser necessário converter para FASTA primeiro
    # Por ora, assumimos que o arquivo já está no formato adequado
    
    PANGOLIN_BIN = "/opt/miniconda3/envs/tcc-metatranscriptomica/bin/pangolin"

    cmd = [PANGOLIN_BIN, input_file, "--outfile", output_csv]

    
    if input_type == "bam":
        # Para BAM, pode precisar de opções adicionais
        pass
    elif input_type == "vcf":
        # Para VCF, pode precisar de opções adicionais
        pass
    
    print(f"Determinando linhagem com Pangolin...")
    print(f"  Input: {Path(input_file).name}")
    print(f"  Type: {input_type}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Pangolin falhou: {result.stderr}")
    
    # Ler resultados
    if Path(output_csv).exists():
        df = pd.read_csv(output_csv)
        print(f"Linhagem determinada. Resultados salvos em: {output_csv}")
        return df
    else:
        raise FileNotFoundError(f"Arquivo de saída não foi criado: {output_csv}")


def generate_consensus_fasta(
    bam_file: str,
    reference_fasta: str,
    output_fasta: str,
    min_depth: int = 10
) -> str:
    """
    Gera genoma consenso a partir de BAM alinhado.
    
    Args:
        bam_file: Arquivo BAM ordenado e indexado
        reference_fasta: Genoma de referência
        output_fasta: Arquivo FASTA consenso de saída
        min_depth: Profundidade mínima para chamar base
    
    Returns:
        str: Caminho do arquivo FASTA consenso
    """
    if not Path(bam_file).exists():
        raise FileNotFoundError(f"Arquivo BAM não encontrado: {bam_file}")
    if not Path(reference_fasta).exists():
        raise FileNotFoundError(f"Genoma de referência não encontrado: {reference_fasta}")
    
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "samtools", "consensus",
        "-f", "fasta",
        "-d", str(min_depth),
        bam_file,
        ">", output_fasta
    ]
    
    print(f"Gerando genoma consenso...")
    result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"samtools consensus falhou: {result.stderr}")
    
    print(f"Genoma consenso gerado: {output_fasta}")
    return output_fasta
