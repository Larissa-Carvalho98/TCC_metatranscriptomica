"""
Módulo para controle de qualidade (QC) de dados FASTQ.
"""

import os
import subprocess
from pathlib import Path
from typing import List, Union, Optional


def run_fastqc(
    input_files: Union[str, List[str]],
    output_dir: str = "results/qc_reports",
    threads: int = 4
) -> dict:
    """
    Executa FastQC em arquivos FASTQ.
    
    Args:
        input_files: Caminho para arquivo FASTQ ou lista de arquivos
        output_dir: Diretório de saída para relatórios
        threads: Número de threads
    
    Returns:
        dict: Dicionário com caminhos dos arquivos de saída
    """
    if isinstance(input_files, str):
        input_files = [input_files]
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    results = {
        "output_dir": str(output_path),
        "html_reports": [],
        "zip_reports": []
    }
    
    for fastq_file in input_files:
        if not Path(fastq_file).exists():
            raise FileNotFoundError(f"Arquivo não encontrado: {fastq_file}")
        
        cmd = [
            "fastqc",
            "-o", str(output_path),
            "-t", str(threads),
            fastq_file
        ]
        
        print(f"Executando FastQC em: {fastq_file}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"FastQC falhou: {result.stderr}")
        
        # Identificar arquivos gerados
        base_name = Path(fastq_file).stem.replace(".fastq", "")
        html_file = output_path / f"{base_name}_fastqc.html"
        zip_file = output_path / f"{base_name}_fastqc.zip"
        
        if html_file.exists():
            results["html_reports"].append(str(html_file))
        if zip_file.exists():
            results["zip_reports"].append(str(zip_file))
    
    return results


def generate_multiqc_report(
    qc_dir: str = "results/qc_reports",
    output_dir: str = "results/qc_reports"
) -> str:
    """
    Gera relatório MultiQC consolidando relatórios FastQC.
    
    Args:
        qc_dir: Diretório contendo relatórios FastQC
        output_dir: Diretório de saída
    
    Returns:
        str: Caminho do relatório HTML gerado
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "multiqc",
        qc_dir,
        "-o", str(output_path),
        "--force"  # Sobrescrever se já existir
    ]
    
    print("Gerando relatório MultiQC...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"MultiQC falhou: {result.stderr}")
    
    report_path = output_path / "multiqc_report.html"
    if report_path.exists():
        print(f"Relatório MultiQC gerado: {report_path}")
        return str(report_path)
    else:
        raise FileNotFoundError("Relatório MultiQC não foi gerado")
