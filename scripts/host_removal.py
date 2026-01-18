"""
Módulo para remoção de reads do hospedeiro (genoma humano).
"""

import subprocess
from pathlib import Path
from typing import Optional


def remove_human_reads_bowtie2(
    input_r1: str,
    input_r2: str,
    reference_index: str,
    output_r1: str,
    output_r2: str,
    aligned_sam: Optional[str] = None,
    threads: int = 8
) -> dict:
    """
    Remove reads que alinham ao genoma humano usando Bowtie2.
    
    Args:
        input_r1: Arquivo FASTQ R1 (paired)
        input_r2: Arquivo FASTQ R2 (paired)
        reference_index: Caminho para índice Bowtie2 do genoma humano (sem extensão)
        output_r1: Arquivo R1 de saída (reads não-humanas)
        output_r2: Arquivo R2 de saída (reads não-humanas)
        aligned_sam: Caminho para arquivo SAM com reads alinhadas (opcional)
        threads: Número de threads
    
    Returns:
        dict: Dicionário com informações sobre o processamento
    """
    # Validar arquivos de entrada
    if not Path(input_r1).exists():
        raise FileNotFoundError(f"Arquivo R1 não encontrado: {input_r1}")
    if not Path(input_r2).exists():
        raise FileNotFoundError(f"Arquivo R2 não encontrado: {input_r2}")
    
    # Criar diretórios de saída
    Path(output_r1).parent.mkdir(parents=True, exist_ok=True)
    Path(output_r2).parent.mkdir(parents=True, exist_ok=True)
    
    if aligned_sam:
        Path(aligned_sam).parent.mkdir(parents=True, exist_ok=True)
    
    # Montar comando Bowtie2
    # --un-conc: Mantém reads não alinhadas em formato paired
    cmd = [
        "bowtie2",
        "-x", reference_index,
        "-1", input_r1,
        "-2", input_r2,
        "--un-conc", f"{Path(output_r1).parent}/{Path(output_r1).stem}_R%.fastq.gz",
        "--threads", str(threads),
        "--very-sensitive-local"  # Modo sensível para detectar mais reads humanas
    ]
    
    if aligned_sam:
        cmd.extend(["-S", aligned_sam])
    else:
        # Descarta SAM se não especificado (apenas mantém reads não-alinhadas)
        cmd.append("--no-unal")  # Não reportar reads não-alinhadas no SAM
    
    print(f"Removendo reads do hospedeiro com Bowtie2...")
    print(f"  Input: {Path(input_r1).name}, {Path(input_r2).name}")
    print(f"  Reference: {reference_index}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Bowtie2 falhou: {result.stderr}")
    
    # Bowtie2 gera arquivos com padrão % (R1, R2)
    # Renomear para nomes especificados se necessário
    temp_r1 = f"{Path(output_r1).parent}/{Path(output_r1).stem}_R1.fastq.gz"
    temp_r2 = f"{Path(output_r1).parent}/{Path(output_r1).stem}_R2.fastq.gz"
    
    if Path(temp_r1).exists() and temp_r1 != output_r1:
        Path(temp_r1).rename(output_r1)
    if Path(temp_r2).exists() and temp_r2 != output_r2:
        Path(temp_r2).rename(output_r2)
    
    results = {
        "output_r1": output_r1,
        "output_r2": output_r2,
        "aligned_sam": aligned_sam
    }
    
    print(f"Remoção de reads do hospedeiro concluída.")
    print(f"  Output R1: {output_r1}")
    print(f"  Output R2: {output_r2}")
    
    return results


def index_reference_bowtie2(
    reference_fasta: str,
    output_prefix: str
) -> str:
    """
    Cria índice Bowtie2 para genoma de referência.
    
    Args:
        reference_fasta: Arquivo FASTA do genoma de referência
        output_prefix: Prefixo para arquivos de índice
    
    Returns:
        str: Caminho do prefixo do índice criado
    """
    if not Path(reference_fasta).exists():
        raise FileNotFoundError(f"Arquivo de referência não encontrado: {reference_fasta}")
    
    Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "bowtie2-build",
        reference_fasta,
        output_prefix,
        "--threads", "4"
    ]
    
    print(f"Indexando genoma de referência com Bowtie2...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"bowtie2-build falhou: {result.stderr}")
    
    print(f"Índice criado: {output_prefix}")
    return output_prefix
