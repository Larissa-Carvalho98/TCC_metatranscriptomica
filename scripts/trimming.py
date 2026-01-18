"""
Módulo para trimming e pré-processamento de reads FASTQ.
"""

import os
import subprocess
from pathlib import Path
from typing import Tuple, Optional


def trim_reads_fastp(
    input_r1: str,
    input_r2: str,
    output_r1: str,
    output_r2: str,
    output_dir: Optional[str] = None,
    json_report: Optional[str] = None,
    html_report: Optional[str] = None,
    detect_adapter: bool = True,
    min_length: int = 36
) -> dict:
    """
    Executa trimming de reads paired-end usando fastp.
    
    Args:
        input_r1: Arquivo FASTQ R1 (forward)
        input_r2: Arquivo FASTQ R2 (reverse)
        output_r1: Arquivo FASTQ R1 de saída (trimmed)
        output_r2: Arquivo FASTQ R2 de saída (trimmed)
        output_dir: Diretório de saída (se None, usa diretório de output_r1)
        json_report: Caminho para relatório JSON (opcional)
        html_report: Caminho para relatório HTML (opcional)
        detect_adapter: Detectar adaptadores automaticamente
        min_length: Comprimento mínimo de reads após trimming
    
    Returns:
        dict: Dicionário com informações sobre o processamento
    """
    # Validar arquivos de entrada
    if not Path(input_r1).exists():
        raise FileNotFoundError(f"Arquivo R1 não encontrado: {input_r1}")
    if not Path(input_r2).exists():
        raise FileNotFoundError(f"Arquivo R2 não encontrado: {input_r2}")
    
    # Criar diretórios de saída se necessário
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        if not json_report:
            json_report = str(Path(output_dir) / "fastp_report.json")
        if not html_report:
            html_report = str(Path(output_dir) / "fastp_report.html")
    else:
        Path(output_r1).parent.mkdir(parents=True, exist_ok=True)
        Path(output_r2).parent.mkdir(parents=True, exist_ok=True)
    
    # Montar comando fastp
    cmd = [
        "fastp",
        "-i", input_r1,
        "-I", input_r2,
        "-o", output_r1,
        "-O", output_r2,
        "--length_required", str(min_length),
    ]
    
    if detect_adapter:
        cmd.append("--detect_adapter_for_pe")
    
    if json_report:
        cmd.extend(["--json", json_report])
    
    if html_report:
        cmd.extend(["--html", html_report])
    
    print(f"Executando fastp em {Path(input_r1).name} e {Path(input_r2).name}...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"fastp falhou: {result.stderr}")
    
    results = {
        "output_r1": output_r1,
        "output_r2": output_r2,
        "json_report": json_report,
        "html_report": html_report
    }
    
    print(f"Trimming concluído. Outputs salvos em:")
    print(f"  R1: {output_r1}")
    print(f"  R2: {output_r2}")
    
    return results


def trim_reads_trimmomatic(
    input_r1: str,
    input_r2: str,
    output_r1_paired: str,
    output_r1_unpaired: str,
    output_r2_paired: str,
    output_r2_unpaired: str,
    adapter_file: str = "TruSeq3-PE.fa",
    leading: int = 3,
    trailing: int = 3,
    window_size: int = 4,
    window_quality: int = 15,
    min_length: int = 36
) -> dict:
    """
    Executa trimming usando Trimmomatic.
    
    Args:
        input_r1: Arquivo FASTQ R1 (forward)
        input_r2: Arquivo FASTQ R2 (reverse)
        output_r1_paired: Arquivo R1 paired de saída
        output_r1_unpaired: Arquivo R1 unpaired de saída
        output_r2_paired: Arquivo R2 paired de saída
        output_r2_unpaired: Arquivo R2 unpaired de saída
        adapter_file: Arquivo de adaptadores Illumina
        leading: Remover bases no início com qualidade < leading
        trailing: Remover bases no final com qualidade < trailing
        window_size: Tamanho da janela deslizante
        window_quality: Qualidade média mínima da janela
        min_length: Comprimento mínimo de reads
    
    Returns:
        dict: Dicionário com informações sobre o processamento
    """
    # Validar arquivos de entrada
    if not Path(input_r1).exists():
        raise FileNotFoundError(f"Arquivo R1 não encontrado: {input_r1}")
    if not Path(input_r2).exists():
        raise FileNotFoundError(f"Arquivo R2 não encontrado: {input_r2}")
    
    # Criar diretórios de saída
    for out_file in [output_r1_paired, output_r1_unpaired, output_r2_paired, output_r2_unpaired]:
        Path(out_file).parent.mkdir(parents=True, exist_ok=True)
    
    # Montar comando Trimmomatic
    cmd = [
        "trimmomatic", "PE",
        input_r1,
        input_r2,
        output_r1_paired,
        output_r1_unpaired,
        output_r2_paired,
        output_r2_unpaired,
        f"ILLUMINACLIP:{adapter_file}:2:30:10",
        f"LEADING:{leading}",
        f"TRAILING:{trailing}",
        f"SLIDINGWINDOW:{window_size}:{window_quality}",
        f"MINLEN:{min_length}"
    ]
    
    print(f"Executando Trimmomatic...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Trimmomatic falhou: {result.stderr}")
    
    results = {
        "output_r1_paired": output_r1_paired,
        "output_r1_unpaired": output_r1_unpaired,
        "output_r2_paired": output_r2_paired,
        "output_r2_unpaired": output_r2_unpaired
    }
    
    print("Trimming com Trimmomatic concluído.")
    return results
