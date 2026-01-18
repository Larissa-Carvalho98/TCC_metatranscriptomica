"""
Módulo para classificação taxonômica usando Kraken2 e Bracken.
"""

import subprocess
import pandas as pd
from pathlib import Path
from typing import Optional


def run_kraken2(
    input_r1: str,
    input_r2: str,
    db_path: str,
    output_file: str,
    report_file: str,
    threads: int = 8,
    paired: bool = True
) -> dict:
    """
    Executa classificação taxonômica usando Kraken2.
    
    Args:
        input_r1: Arquivo FASTQ R1
        input_r2: Arquivo FASTQ R2 (ou None se single-end)
        db_path: Caminho para banco de dados Kraken2
        output_file: Arquivo de saída com classificação de reads
        report_file: Arquivo de relatório taxonômico
        threads: Número de threads
        paired: Se True, trata como paired-end
    
    Returns:
        dict: Dicionário com caminhos dos arquivos gerados
    """
    if not Path(input_r1).exists():
        raise FileNotFoundError(f"Arquivo R1 não encontrado: {input_r1}")
    
    if paired and input_r2 and not Path(input_r2).exists():
        raise FileNotFoundError(f"Arquivo R2 não encontrado: {input_r2}")
    
    # Criar diretórios de saída
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    Path(report_file).parent.mkdir(parents=True, exist_ok=True)
    
    # Montar comando Kraken2
    cmd = [
        "kraken2",
        "--db", db_path,
        "--output", output_file,
        "--report", report_file,
        "--threads", str(threads)
    ]
    
    if paired and input_r2:
        cmd.extend(["--paired", input_r1, input_r2])
    else:
        cmd.append(input_r1)
    
    print(f"Executando Kraken2...")
    print(f"  Database: {db_path}")
    print(f"  Input: {Path(input_r1).name}" + (f", {Path(input_r2).name}" if input_r2 else ""))
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Kraken2 falhou: {result.stderr}")
    
    results = {
        "output_file": output_file,
        "report_file": report_file
    }
    
    print(f"Kraken2 concluído.")
    print(f"  Relatório: {report_file}")
    
    return results


def run_bracken(
    kraken_report: str,
    db_path: str,
    output_file: str,
    level: str = "S",  # S = species
    read_length: int = 100
) -> str:
    """
    Executa Bracken para estimar abundância taxonômica.
    
    Args:
        kraken_report: Relatório Kraken2 de entrada
        db_path: Caminho para banco de dados Kraken2 (mesmo usado no Kraken2)
        output_file: Arquivo de saída Bracken
        level: Nível taxonômico (S=species, G=genus, etc.)
        read_length: Comprimento médio das reads
    
    Returns:
        str: Caminho do arquivo de saída
    
    Raises:
        FileNotFoundError: Se Bracken não estiver instalado ou relatório não encontrado
        RuntimeError: Se Bracken falhar na execução
    """
    import shutil
    
    # Verificar se Bracken está instalado
    if shutil.which("bracken") is None:
        raise FileNotFoundError(
            "Bracken não está instalado ou não está no PATH.\n"
            "Para instalar: conda install -c bioconda bracken\n"
            "Nota: Bracken é opcional - você pode usar apenas os resultados do Kraken2."
        )
    
    if not Path(kraken_report).exists():
        raise FileNotFoundError(f"Relatório Kraken2 não encontrado: {kraken_report}")
    
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "bracken",
        "-d", db_path,
        "-i", kraken_report,
        "-o", output_file,
        "-r", str(read_length),
        "-l", level
    ]
    
    print(f"Executando Bracken...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Bracken falhou: {result.stderr}")
    
    print(f"Bracken concluído. Output: {output_file}")
    return output_file


def parse_kraken2_report(report_file: str) -> pd.DataFrame:
    """
    Parseia relatório Kraken2 em DataFrame pandas.
    
    Args:
        report_file: Caminho para relatório Kraken2
    
    Returns:
        pd.DataFrame: DataFrame com colunas: percent, reads, tax_reads, rank, taxid, name
    """
    if not Path(report_file).exists():
        raise FileNotFoundError(f"Relatório não encontrado: {report_file}")
    
    df = pd.read_csv(
        report_file,
        sep='\t',
        names=['percent', 'reads', 'tax_reads', 'rank', 'taxid', 'name']
    )
    
    return df


def parse_bracken_report(report_file: str) -> pd.DataFrame:
    """
    Parseia relatório Bracken em DataFrame pandas.
    
    Args:
        report_file: Caminho para relatório Bracken
    
    Returns:
        pd.DataFrame: DataFrame com abundância estimada
    """
    if not Path(report_file).exists():
        raise FileNotFoundError(f"Relatório não encontrado: {report_file}")
    
    df = pd.read_csv(report_file, sep='\t')
    
    return df
