# TCC Metatranscriptômica & Machine Learning

Pipeline completo para análise metatranscriptômica e predição de COVID-19 usando Machine Learning.

## Estrutura do Projeto

```
TCC_metatrascriptomica/
├── data/                    # Dados (raw, processed, references, training)
├── notebooks/               # Jupyter Notebooks (orquestração)
├── scripts/                 # Scripts Python modulares (lógica principal)
├── results/                 # Resultados das análises
├── environment.yml          # Configuração do ambiente conda
├── requirements.txt         # Dependências Python
└── setup_environment.sh     # Script de setup automatizado
```

## Requisitos

- Python 3.9 ou 3.10
- Conda (recomendado) ou venv
- Ferramentas bioinformáticas (instaladas via conda)

## Instalação

### Opção 1: Usando Conda (Recomendado)

```bash
# Criar ambiente conda
conda env create -f environment.yml

# Ativar ambiente
conda activate tcc-metatranscriptomica

# Ou usar script automatizado
bash setup_environment.sh conda
```

### Opção 2: Usando venv

```bash
# Criar ambiente virtual
python3.10 -m venv venv
source venv/bin/activate  # No Windows: venv\Scripts\activate

# Instalar dependências Python
pip install -r requirements.txt

# Nota: Ferramentas bioinformáticas precisam ser instaladas separadamente
# (usando conda ou instalação manual)
```

## Uso

### 1. Preparar Dados

Baixar dados FASTQ do paciente e genomas de referência:

```bash
# Criar diretórios
mkdir -p data/raw data/references data/training

# Baixar dados FASTQ (exemplo)
wget https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R1_001.fastq.gz -P data/raw/
wget https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R2_001.fastq.gz -P data/raw/

# Baixar genoma de referência SARS-CoV-2
wget https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1798174254 -O data/references/NC_045512.2.fasta
```

### 2. Executar Notebooks

Os notebooks devem ser executados sequencialmente:

1. **01_metatranscriptomica_pipeline.ipynb**: Pipeline de metatranscriptômica
2. **02_variant_calling_pipeline.ipynb**: Chamada de variantes do SARS-CoV-2
3. **03_feature_engineering.ipynb**: Extração e engenharia de features
4. **04_ml_training_prediction.ipynb**: Treinamento ML e predição

```bash
# Iniciar Jupyter
jupyter notebook notebooks/
```

## Arquitetura

### Notebooks (.ipynb)
- **Função**: Camada de orquestração
- Execução passo a passo
- Visualizações interativas
- Documentação inline

### Scripts Python (.py)
- **Função**: Lógica principal modular
- Funções reutilizáveis e testáveis
- Importados nos notebooks

**Exemplo de uso:**
```python
# No notebook
from scripts.qc_processing import run_fastqc, generate_multiqc_report
from scripts.trimming import trim_reads_fastp

# Executar funções
run_fastqc([r1_file, r2_file], output_dir="results/qc_reports")
```

## Caminhos Relativos

Todos os caminhos são relativos à raiz do projeto:
- ✅ `data/raw/file.fastq.gz`
- ❌ `/Users/larissa/Desktop/TCC_metatrascriptomica/data/raw/file.fastq.gz`

## Migração para Google Colab

Veja [COLAB_MIGRATION.md](COLAB_MIGRATION.md) para instruções detalhadas.

## Dependências Principais

### Ferramentas Bioinformáticas (via conda)
- FastQC, MultiQC
- Trimmomatic, fastp
- Bowtie2, BWA
- Kraken2, Bracken
- samtools, bcftools
- SnpEff, Pangolin

### Bibliotecas Python (via pip)
- pandas, numpy, scipy
- scikit-learn, xgboost, shap
- matplotlib, seaborn
- pysam, tqdm

## Estrutura de Resultados

```
results/
├── qc_reports/          # Relatórios de qualidade
├── kraken2_reports/     # Classificação taxonômica
├── variants/            # Variantes e linhagens
├── features/            # Vetores de features
└── ml_results/          # Modelos e predições
```

## Documentação

- **[GUIA_EXECUCAO_LOCAL.md](GUIA_EXECUCAO_LOCAL.md)** - Guia passo a passo para execução local
- **[COLAB_MIGRATION.md](COLAB_MIGRATION.md)** - Guia de migração para Google Colab
- [guia-execucao-tcc.md](guia-execucao-tcc.md) - Guia técnico detalhado do pipeline
- [descricao-tcc.txt](descricao-tcc.txt) - Descrição do projeto e requisitos

## Licença

Este projeto é parte de um trabalho de conclusão de curso (TCC).
