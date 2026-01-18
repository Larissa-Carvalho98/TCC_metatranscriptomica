# Guia de Migração para Google Colab

Este guia explica como migrar o projeto desenvolvido localmente para o Google Colab, mantendo total compatibilidade.

## Pré-requisitos

1. Conta Google (para acessar Colab)
2. Projeto desenvolvido localmente com todos os arquivos necessários
3. Acesso aos dados (FASTQ, matriz de treinamento, referências)

## Passo 1: Upload do Projeto

### Opção A: Upload Manual (Recomendado para pequenos projetos)

1. No Colab, crie um novo notebook
2. Execute para montar o Google Drive:
   ```python
   from google.colab import drive
   drive.mount('/content/drive')
   ```
3. Faça upload dos arquivos do projeto via interface do Colab:
   - Clique em "Arquivos" (ícone de pasta) no painel esquerdo
   - Clique em "Fazer upload" e selecione os arquivos/diretórios

### Opção B: Clone do GitHub/GitLab (Recomendado)

Se o projeto estiver em um repositório Git:

```python
# Clonar repositório
!git clone https://github.com/seu-usuario/TCC_metatrascriptomica.git
```

### Opção C: Upload via Terminal

```python
# Instalar gdown (se necessário)
!pip install gdown

# Fazer upload compactado ou usar wget/curl
# Exemplo: fazer upload do diretório scripts/ via interface, ou:
```

## Passo 2: Configurar Estrutura de Diretórios

No início do notebook do Colab, execute:

```python
import os
from pathlib import Path

# Definir diretório raiz do projeto
# Se fez upload manual, ajuste o caminho:
PROJECT_ROOT = Path("/content/TCC_metatrascriptomica")

# Se clonou do Git:
# PROJECT_ROOT = Path("/content/TCC_metatrascriptomica")

# Criar estrutura de diretórios
directories = [
    PROJECT_ROOT / "data" / "raw",
    PROJECT_ROOT / "data" / "processed",
    PROJECT_ROOT / "data" / "references",
    PROJECT_ROOT / "data" / "training",
    PROJECT_ROOT / "results" / "qc_reports",
    PROJECT_ROOT / "results" / "kraken2_reports",
    PROJECT_ROOT / "results" / "variants",
    PROJECT_ROOT / "results" / "features",
    PROJECT_ROOT / "results" / "ml_results"
]

for directory in directories:
    directory.mkdir(parents=True, exist_ok=True)

# Mudar para diretório do projeto
os.chdir(PROJECT_ROOT)
print(f"Diretório de trabalho: {os.getcwd()}")
```

## Passo 3: Instalar Dependências Python

```python
# Instalar dependências Python
!pip install -r requirements.txt
```

**Nota**: O `requirements.txt` já contém todas as dependências Python necessárias.

## Passo 4: Instalar Ferramentas Bioinformáticas

O Colab não vem com ferramentas bioinformáticas pré-instaladas. Instale usando conda:

```python
# Instalar miniconda no Colab
!wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
!chmod +x Miniconda3-latest-Linux-x86_64.sh
!bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /usr/local
import sys
sys.path.append('/usr/local/lib/python3.10/site-packages')

# Configurar conda
!conda install -c conda-forge -y conda
!conda config --add channels bioconda
!conda config --add channels conda-forge

# Instalar ferramentas bioinformáticas (uma por uma ou via environment.yml)
!conda install -c bioconda -y fastqc multiqc trimmomatic fastp bowtie2 bwa kraken2 bracken samtools bcftools snpeff pangolin

# Ou criar ambiente completo (mais lento mas mais completo)
# !conda env create -f environment.yml
```

**Alternativa**: Instalação seletiva (apenas ferramentas necessárias no momento):

```python
# Instalar apenas o necessário para cada etapa
# Exemplo para FastQC:
!conda install -c bioconda -y fastqc

# Exemplo para Kraken2:
!conda install -c bioconda -y kraken2
```

## Passo 5: Upload de Dados

### Upload de Arquivos FASTQ

Para arquivos grandes (FASTQ), use Google Drive ou download direto:

```python
# Opção 1: Download direto (se disponível via URL)
!wget https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R1_001.fastq.gz -P data/raw/
!wget https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R2_001.fastq.gz -P data/raw/

# Opção 2: Upload via interface do Colab (limite ~100MB por arquivo)
# Usar a interface gráfica do Colab

# Opção 3: Google Drive (para arquivos grandes)
from google.colab import drive
drive.mount('/content/drive')

# Copiar do Drive para o projeto
!cp /content/drive/MyDrive/dados_TCC/*.fastq.gz data/raw/
```

### Upload de Matriz de Treinamento

```python
# Upload manual via interface ou:
!wget [URL_DO_ARQUIVO] -O data/training/pivoted-virome-organisms-atleast10tpm-species-covid-TCC-pos.csv
```

### Upload de Genomas de Referência

```python
# Download genoma SARS-CoV-2
!wget https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1798174254 -O data/references/NC_045512.2.fasta

# Download genoma humano (se necessário)
# !wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O data/references/hg38.fa.gz
# !gunzip data/references/hg38.fa.gz
```

## Passo 6: Ajustar Caminhos nos Notebooks

Os notebooks já usam caminhos relativos, mas no Colab você precisa garantir que o diretório de trabalho está correto.

**Adicione no início de cada notebook:**

```python
import os
from pathlib import Path

# Mudar para diretório do projeto
PROJECT_ROOT = Path("/content/TCC_metatrascriptomica")  # Ajustar se necessário
os.chdir(PROJECT_ROOT)

# Verificar
print(f"Diretório atual: {os.getcwd()}")
print(f"Scripts disponíveis: {list(Path('scripts').glob('*.py'))}")
```

## Passo 7: Executar Notebooks

Os notebooks desenvolvidos localmente funcionam diretamente no Colab sem modificações estruturais, desde que:

1. ✅ Caminhos são relativos
2. ✅ Scripts Python estão disponíveis (via upload/clone)
3. ✅ Dependências foram instaladas
4. ✅ Diretório de trabalho está configurado

**Execute os notebooks na ordem:**
1. `01_metatranscriptomica_pipeline.ipynb`
2. `02_variant_calling_pipeline.ipynb`
3. `03_feature_engineering.ipynb`
4. `04_ml_training_prediction.ipynb`

## Notas Importantes

### Limitações do Colab

1. **Tempo de execução**: Notebooks desconectam após inatividade (gratuito: ~90 min)
   - **Solução**: Use runtime de GPU/TPU (mais tempo) ou salve progresso intermediário

2. **Memória**: ~12-15GB RAM (gratuito), 32GB+ (Colab Pro)
   - **Solução**: Processe dados em lotes ou use runtime de alta RAM

3. **Armazenamento**: ~77GB espaço temporário
   - **Solução**: Salve resultados importantes no Google Drive

4. **Bancos de Dados Kraken2**: Grandes (~100GB+)
   - **Solução**: Use bancos pré-construídos menores (apenas vírus) ou armazene no Drive

### Dicas para Otimização

```python
# Salvar progresso intermediário no Drive
from google.colab import drive
drive.mount('/content/drive')

# Copiar resultados importantes
!cp -r results/ /content/drive/MyDrive/TCC_results_backup/

# Baixar resultados locais
from google.colab import files
files.download('results/ml_results/best_model.pkl')
```

### Verificação de Instalação

```python
# Verificar ferramentas instaladas
import subprocess

tools = ['fastqc', 'kraken2', 'samtools', 'bcftools', 'pangolin']
for tool in tools:
    try:
        result = subprocess.run([tool, '--version'], capture_output=True, text=True)
        print(f"✅ {tool}: {result.stdout.split()[0]}")
    except:
        print(f"❌ {tool}: Não instalado")
```

## Estrutura de Células no Colab

Recomenda-se organizar o notebook Colab assim:

```markdown
# Célula 1: Configuração e Instalações
# - Montar Drive
# - Criar diretórios
# - Instalar dependências

# Célula 2: Upload de Dados
# - Dados FASTQ
# - Referências
# - Matriz de treinamento

# Célula 3+: Execução do Pipeline
# - Importar scripts
# - Executar funções
# - Visualizações
```

## Troubleshooting

### Erro: "ModuleNotFoundError: No module named 'scripts'"

**Solução:**
```python
import sys
from pathlib import Path
PROJECT_ROOT = Path("/content/TCC_metatrascriptomica")
sys.path.insert(0, str(PROJECT_ROOT))
```

### Erro: "Command not found: fastqc"

**Solução:** Instalar via conda (veja Passo 4)

### Erro: "Permission denied" ao executar ferramentas

**Solução:**
```python
# Adicionar conda ao PATH
import os
os.environ['PATH'] = '/usr/local/bin:' + os.environ['PATH']
```

### Erro: "Out of memory"

**Solução:** 
- Use runtime de alta RAM (Colab Pro)
- Processe dados em lotes menores
- Libere memória: `del variavel_grande`

## Conclusão

Com este guia, você deve conseguir executar todo o pipeline no Google Colab sem modificações estruturais nos notebooks. O código desenvolvido localmente é totalmente compatível com o Colab!
