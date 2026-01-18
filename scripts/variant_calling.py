"""
Módulo para alinhamento e chamada de variantes do SARS-CoV-2.
"""

import subprocess
from pathlib import Path
from typing import Optional


def align_reads_bwa(
    input_r1: str,
    input_r2: str,
    reference_fasta: str,
    output_sam: str,
    threads: int = 8
) -> str:
    """
    Alinha reads ao genoma de referência usando BWA-MEM.
    
    Args:
        input_r1: Arquivo FASTQ R1
        input_r2: Arquivo FASTQ R2
        reference_fasta: Genoma de referência FASTA (SARS-CoV-2)
        output_sam: Arquivo SAM de saída
        threads: Número de threads
    
    Returns:
        str: Caminho do arquivo SAM gerado
    """
    if not Path(input_r1).exists():
        raise FileNotFoundError(f"Arquivo R1 não encontrado: {input_r1}")
    if not Path(input_r2).exists():
        raise FileNotFoundError(f"Arquivo R2 não encontrado: {input_r2}")
    if not Path(reference_fasta).exists():
        raise FileNotFoundError(f"Genoma de referência não encontrado: {reference_fasta}")
    
    Path(output_sam).parent.mkdir(parents=True, exist_ok=True)
    
    # Verificar se índice BWA existe
    if not Path(f"{reference_fasta}.amb").exists():
        print("Índice BWA não encontrado. Criando índice...")
        index_reference_bwa(reference_fasta)
    
    cmd = [
        "bwa", "mem",
        "-t", str(threads),
        reference_fasta,
        input_r1,
        input_r2,
        ">", output_sam
    ]
    
    print(f"Alinhando reads com BWA-MEM...")
    result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"BWA-MEM falhou: {result.stderr}")
    
    print(f"Alinhamento concluído: {output_sam}")
    return output_sam


def index_reference_bwa(reference_fasta: str) -> str:
    """
    Cria índice BWA para genoma de referência.
    
    Args:
        reference_fasta: Arquivo FASTA do genoma de referência
    
    Returns:
        str: Caminho do arquivo FASTA indexado
    """
    if not Path(reference_fasta).exists():
        raise FileNotFoundError(f"Arquivo de referência não encontrado: {reference_fasta}")
    
    cmd = ["bwa", "index", reference_fasta]
    
    print(f"Indexando genoma de referência com BWA...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"bwa index falhou: {result.stderr}")
    
    print(f"Índice BWA criado para: {reference_fasta}")
    return reference_fasta


def sam_to_bam(sorted_bam: str, sam_file: str, output_bam: Optional[str] = None) -> str:
    """
    Converte SAM para BAM, ordena e indexa.
    
    Args:
        sorted_bam: Caminho do BAM ordenado de saída (sem extensão .bam)
        sam_file: Arquivo SAM de entrada
        output_bam: Caminho do BAM não ordenado (opcional, se None usa sorted_bam)
    
    Returns:
        str: Caminho do BAM ordenado e indexado
    """
    if not Path(sam_file).exists():
        raise FileNotFoundError(f"Arquivo SAM não encontrado: {sam_file}")
    
    output_path = Path(sorted_bam)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # 1. Converter SAM para BAM
    if output_bam is None:
        temp_bam = str(output_path.with_suffix('.temp.bam'))
    else:
        temp_bam = output_bam
        Path(temp_bam).parent.mkdir(parents=True, exist_ok=True)
    
    cmd1 = ["samtools", "view", "-bS", sam_file, "-o", temp_bam]
    print("Convertendo SAM para BAM...")
    result = subprocess.run(cmd1, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"samtools view falhou: {result.stderr}")
    
    # 2. Ordenar BAM
    cmd2 = ["samtools", "sort", temp_bam, "-o", sorted_bam]
    print("Ordenando BAM...")
    result = subprocess.run(cmd2, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"samtools sort falhou: {result.stderr}")
    
    # 3. Indexar BAM
    cmd3 = ["samtools", "index", sorted_bam]
    print("Indexando BAM...")
    result = subprocess.run(cmd3, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"samtools index falhou: {result.stderr}")
    
    # Remover BAM temporário
    if Path(temp_bam).exists() and temp_bam != sorted_bam:
        Path(temp_bam).unlink()
    
    print(f"BAM processado: {sorted_bam}")
    return sorted_bam


def call_variants_bcftools(
    bam_file: str,
    reference_fasta: str,
    output_vcf: str
) -> str:
    """
    Chama variantes usando bcftools.
    
    Args:
        bam_file: Arquivo BAM ordenado e indexado
        reference_fasta: Genoma de referência FASTA
        output_vcf: Arquivo VCF de saída
    
    Returns:
        str: Caminho do arquivo VCF gerado
    """
    if not Path(bam_file).exists():
        raise FileNotFoundError(f"Arquivo BAM não encontrado: {bam_file}")
    if not Path(reference_fasta).exists():
        raise FileNotFoundError(f"Genoma de referência não encontrado: {reference_fasta}")
    
    Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)
    
    # Indexar referência se necessário
    if not Path(f"{reference_fasta}.fai").exists():
        cmd_idx = ["samtools", "faidx", reference_fasta]
        subprocess.run(cmd_idx, capture_output=True, text=True)
    
    # Chamada de variantes (mpileup + call)
    cmd = [
        "bcftools", "mpileup",
        "-f", reference_fasta,
        bam_file,
        "|",
        "bcftools", "call",
        "-mv",
        "-Oz",
        "-o", output_vcf
    ]
    
    print(f"Chamando variantes com bcftools...")
    result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"bcftools falhou: {result.stderr}")
    
    # Indexar VCF
    cmd_idx = ["bcftools", "index", output_vcf]
    result = subprocess.run(cmd_idx, capture_output=True, text=True)
    
    print(f"Variantes chamadas: {output_vcf}")
    return output_vcf
