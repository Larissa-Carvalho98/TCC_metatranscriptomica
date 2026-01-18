# Guia de Execução Passo a Passo - TCC Metatranscriptômica & Machine Learning

## Índice
1. [Visão Geral](#visão-geral)
2. [Etapa 1: Pipeline de Metatranscriptômica (Viroma)](#etapa-1-pipeline-de-metatranscriptômica-viroma)
3. [Etapa 2: Pipeline de Chamada de Variantes do SARS-CoV-2](#etapa-2-pipeline-de-chamada-de-variantes-do-sars-cov-2)
4. [Etapa 3: Extração e Engenharia de Features](#etapa-3-extração-e-engenharia-de-features)
5. [Etapa 4: Treinamento e Predição com Machine Learning](#etapa-4-treinamento-e-predição-com-machine-learning)
6. [Recursos e Referências](#recursos-e-referências)

---

## Visão Geral

Este guia detalha a execução completa do TCC, dividido em 4 entregas principais:
- **Entrega 1**: Pipeline de metatranscriptômica (viroma)
- **Entrega 2**: Pipeline de chamada de variantes do SARS-CoV-2
- **Entrega 3**: Pipeline de extração e engenharia de features
- **Entrega 4**: Treinamento e predição com Machine Learning

---

## Etapa 1: Pipeline de Metatranscriptômica (Viroma)

### Objetivo
Processar os arquivos FASTQ brutos da amostra do paciente para identificar organismos presentes (principalmente vírus), confirmar a presença de SARS-CoV-2 e avaliar coinfecções.

### Dados de Entrada

**Arquivos FASTQ do Paciente João:**
- `patient_joao_VIROMA_S21_R1_001.fastq.gz` (Reads forward)
- `patient_joao_VIROMA_S21_R2_001.fastq.gz` (Reads reverse)

**Links de Download:**
```
https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R1_001.fastq.gz
https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R2_001.fastq.gz
```

**Banco de Dados Necessários:**
- Genoma humano de referência (HG38) - para remoção de reads do hospedeiro
- Banco de dados Kraken2 (vírus/bactérias) - para classificação taxonômica

### Passo a Passo Detalhado

#### 1.1 Controle de Qualidade (QC)

**Ferramenta: FastQC**

**Instalação:**
```bash
# Via conda
conda install -c bioconda fastqc

# Ou via apt (Linux)
sudo apt-get install fastqc
```

**Como Utilizar:**
```bash
# Análise de qualidade para cada arquivo
fastqc patient_joao_VIROMA_S21_R1_001.fastq.gz -o qc_reports/
fastqc patient_joao_VIROMA_S21_R2_001.fastq.gz -o qc_reports/
```

**Saída Esperada:**
- Arquivos HTML com relatórios de qualidade
- Métricas: Phred scores, conteúdo GC, comprimento de reads, presença de adaptadores

**Ferramenta: MultiQC** (para consolidar relatórios)

**Instalação:**
```bash
conda install -c bioconda multiqc
```

**Como Utilizar:**
```bash
multiqc qc_reports/ -o multiqc_report/
```

**Saída Esperada:**
- Relatório HTML consolidado com todas as métricas de qualidade

---

#### 1.2 Pré-processamento (Trimming)

**Ferramenta: Trimmomatic**

**Instalação:**
```bash
conda install -c bioconda trimmomatic
```

**Como Utilizar:**
```bash
trimmomatic PE \
  patient_joao_VIROMA_S21_R1_001.fastq.gz \
  patient_joao_VIROMA_S21_R2_001.fastq.gz \
  patient_joao_R1_paired.fastq.gz \
  patient_joao_R1_unpaired.fastq.gz \
  patient_joao_R2_paired.fastq.gz \
  patient_joao_R2_unpaired.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 \
  TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:36
```

**Parâmetros:**
- `ILLUMINACLIP`: Remove adaptadores Illumina
- `LEADING:3`: Remove bases com qualidade < 3 no início
- `TRAILING:3`: Remove bases com qualidade < 3 no final
- `SLIDINGWINDOW:4:15`: Janela deslizante de 4 bases, média de qualidade mínima 15
- `MINLEN:36`: Remove reads com menos de 36 bases

**Alternativa: fastp** (mais rápido)

**Instalação:**
```bash
conda install -c bioconda fastp
```

**Como Utilizar:**
```bash
fastp \
  -i patient_joao_VIROMA_S21_R1_001.fastq.gz \
  -I patient_joao_VIROMA_S21_R2_001.fastq.gz \
  -o patient_joao_R1_trimmed.fastq.gz \
  -O patient_joao_R2_trimmed.fastq.gz \
  --detect_adapter_for_pe \
  --cut_front \
  --cut_tail \
  --trim_poly_g \
  --length_required 36
```

**Saída Esperada:**
- Arquivos FASTQ limpos (paired e unpaired)
- Relatório HTML com estatísticas de trimming

---

#### 1.3 Remoção de Reads do Hospedeiro (Humano)

**Ferramenta: Bowtie2**

**Instalação:**
```bash
conda install -c bioconda bowtie2
```

**Onde Obter o Genoma de Referência Humano (HG38):**
```bash
# Download do genoma humano
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Ou usar versão já indexada
# https://genome-idx.s3.amazonaws.com/bt/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.tar.gz
```

**Como Utilizar:**
```bash
# 1. Indexar o genoma humano
bowtie2-build hg38.fa hg38_index

# 2. Alinhar e manter apenas reads não-alinhadas
bowtie2 -x hg38_index \
  -1 patient_joao_R1_paired.fastq.gz \
  -2 patient_joao_R2_paired.fastq.gz \
  --un-conc patient_joao_nonhuman_R%.fastq.gz \
  -S patient_joao_aligned_human.sam \
  --threads 8

# 3. Converter SAM para BAM (opcional, para inspeção)
samtools view -bS patient_joao_aligned_human.sam > patient_joao_aligned_human.bam
```

**Saída Esperada:**
- `patient_joao_nonhuman_R1.fastq.gz` e `patient_joao_nonhuman_R2.fastq.gz` (reads não-humanas)

---

#### 1.4 Classificação Taxonômica

**Ferramenta: Kraken2**

**Instalação:**
```bash
conda install -c bioconda kraken2
```

**Onde Obter o Banco de Dados Kraken2:**

**Opção 1: Banco de Vírus (mais rápido, menor)**
```bash
# Download do banco de vírus
kraken2-build --download-library viral --db viral_db
kraken2-build --build --db viral_db
```

**Opção 2: Banco Completo (vírus + bactérias + archaea)**
```bash
# Download completo (requer ~100GB de espaço)
kraken2-build --download-library viral --db full_db
kraken2-build --download-library bacteria --db full_db
kraken2-build --download-library archaea --db full_db
kraken2-build --build --db full_db
```

**Opção 3: Usar banco pré-construído (recomendado para Colab)**
```bash
# Links para bancos pré-construídos:
# https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20221209.tar.gz (vírus + bactérias)
# https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz (apenas vírus)
```

**Como Utilizar:**
```bash
kraken2 \
  --db viral_db \
  --paired \
  --threads 8 \
  --output patient_joao_kraken2_output.txt \
  --report patient_joao_kraken2_report.txt \
  patient_joao_nonhuman_R1.fastq.gz \
  patient_joao_nonhuman_R2.fastq.gz
```

**Saída Esperada:**
- `patient_joao_kraken2_output.txt`: Classificação de cada read
- `patient_joao_kraken2_report.txt`: Relatório agregado por táxon (abundância)

**Ferramenta: Bracken** (para estimar abundância mais precisa)

**Instalação:**
```bash
conda install -c bioconda bracken
```

**Como Utilizar:**
```bash
# 1. Construir banco Bracken a partir do banco Kraken2
bracken-build -d viral_db -t 8 -k 35 -l 100

# 2. Estimar abundância
bracken \
  -d viral_db \
  -i patient_joao_kraken2_report.txt \
  -o patient_joao_bracken_output.txt \
  -r 100 \
  -l S
```

**Saída Esperada:**
- `patient_joao_bracken_output.txt`: Abundância estimada por espécie

---

#### 1.5 Visualização

**Ferramenta: Krona** (visualização interativa)

**Instalação:**
```bash
conda install -c bioconda krona
```

**Como Utilizar:**
```bash
# Converter relatório Kraken2 para formato Krona
ktImportTaxonomy patient_joao_kraken2_report.txt -o patient_joao_krona.html
```

**Saída Esperada:**
- Arquivo HTML interativo com visualização hierárquica da composição taxonômica

**Alternativa: Python (matplotlib/seaborn)**
```python
import pandas as pd
import matplotlib.pyplot as plt

# Ler relatório Kraken2
df = pd.read_csv('patient_joao_kraken2_report.txt', sep='\t', 
                 names=['percent', 'reads', 'tax_reads', 'rank', 'taxid', 'name'])

# Filtrar top organismos
top_orgs = df[df['rank'] == 'S'].nlargest(10, 'reads')

# Gráfico de barras
plt.figure(figsize=(10, 6))
plt.barh(top_orgs['name'], top_orgs['reads'])
plt.xlabel('Número de Reads')
plt.title('Top 10 Organismos Detectados')
plt.tight_layout()
plt.savefig('top_organisms.png')
```

---

### Saídas Finais da Etapa 1

1. **Relatórios de Qualidade:**
   - Relatório MultiQC consolidado
   - Estatísticas de trimming

2. **Dados Processados:**
   - Arquivos FASTQ limpos e sem reads humanas
   - Relatório Kraken2 com classificação taxonômica
   - Relatório Bracken com abundância estimada

3. **Visualizações:**
   - Gráfico Krona interativo
   - Gráficos de barras dos top organismos

4. **Interpretação:**
   - Confirmação da presença de SARS-CoV-2
   - Lista de coinfecções detectadas
   - Abundância relativa de cada organismo

---

## Etapa 2: Pipeline de Chamada de Variantes do SARS-CoV-2

### Objetivo
Reconstruir o genoma viral presente na amostra e identificar mutações em relação ao genoma de referência, determinando a linhagem.

### Dados de Entrada

**Arquivos FASTQ Processados:**
- `patient_joao_nonhuman_R1.fastq.gz`
- `patient_joao_nonhuman_R2.fastq.gz`

**Genoma de Referência SARS-CoV-2:**
- Genoma Wuhan (NC_045512.2)

**Onde Obter:**
```bash
# Download direto do NCBI
wget https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1798174254 -O NC_045512.2.fasta

# Ou via efetch
efetch -db nuccore -id NC_045512.2 -format fasta > NC_045512.2.fasta
```

---

### Passo a Passo Detalhado

#### 2.1 Alinhamento ao Genoma de Referência

**Ferramenta: BWA-MEM**

**Instalação:**
```bash
conda install -c bioconda bwa
```

**Como Utilizar:**
```bash
# 1. Indexar o genoma de referência
bwa index NC_045512.2.fasta

# 2. Alinhar reads ao genoma de referência
bwa mem -t 8 \
  NC_045512.2.fasta \
  patient_joao_nonhuman_R1.fastq.gz \
  patient_joao_nonhuman_R2.fastq.gz \
  > patient_joao_aligned.sam
```

**Alternativa: Bowtie2**
```bash
# 1. Indexar
bowtie2-build NC_045512.2.fasta sars_cov2_index

# 2. Alinhar
bowtie2 -x sars_cov2_index \
  -1 patient_joao_nonhuman_R1.fastq.gz \
  -2 patient_joao_nonhuman_R2.fastq.gz \
  -S patient_joao_aligned.sam \
  --threads 8
```

**Saída Esperada:**
- Arquivo SAM com alinhamentos

---

#### 2.2 Processamento do BAM

**Ferramenta: samtools**

**Instalação:**
```bash
conda install -c bioconda samtools
```

**Como Utilizar:**
```bash
# 1. Converter SAM para BAM
samtools view -bS patient_joao_aligned.sam > patient_joao_aligned.bam

# 2. Ordenar por coordenada
samtools sort patient_joao_aligned.bam -o patient_joao_sorted.bam

# 3. Indexar
samtools index patient_joao_sorted.bam

# 4. Estatísticas de cobertura
samtools depth patient_joao_sorted.bam > patient_joao_coverage.txt
samtools stats patient_joao_sorted.bam > patient_joao_stats.txt
```

**Ferramenta: Picard** (remoção de duplicatas)

**Instalação:**
```bash
conda install -c bioconda picard
```

**Como Utilizar:**
```bash
picard MarkDuplicates \
  I=patient_joao_sorted.bam \
  O=patient_joao_dedup.bam \
  M=patient_joao_duplicates_metrics.txt \
  REMOVE_DUPLICATES=true

# Re-indexar após remoção de duplicatas
samtools index patient_joao_dedup.bam
```

**Saída Esperada:**
- BAM ordenado e indexado
- Estatísticas de cobertura e profundidade
- BAM sem duplicatas

---

#### 2.3 Chamada de Variantes

**Opção 1: BCFtools**

**Instalação:**
```bash
conda install -c bioconda bcftools
```

**Como Utilizar:**
```bash
# 1. Indexar genoma de referência (se ainda não feito)
samtools faidx NC_045512.2.fasta

# 2. Chamada de variantes
bcftools mpileup -f NC_045512.2.fasta \
  patient_joao_dedup.bam \
  | bcftools call -mv -Oz -o patient_joao_variants.vcf.gz

# 3. Indexar VCF
bcftools index patient_joao_variants.vcf.gz

# 4. Visualizar variantes
bcftools view patient_joao_variants.vcf.gz | less
```

**Opção 2: FreeBayes**

**Instalação:**
```bash
conda install -c bioconda freebayes
```

**Como Utilizar:**
```bash
freebayes \
  -f NC_045512.2.fasta \
  -b patient_joao_dedup.bam \
  > patient_joao_variants_freebayes.vcf
```

**Opção 3: LoFreq** (especializado em baixa frequência)

**Instalação:**
```bash
conda install -c bioconda lofreq
```

**Como Utilizar:**
```bash
# 1. Calibrar alinhamento
lofreq viterbi -f NC_045512.2.fasta patient_joao_dedup.bam | \
  samtools sort - > patient_joao_calibrated.bam
samtools index patient_joao_calibrated.bam

# 2. Chamada de variantes
lofreq call-parallel \
  --pp-threads 8 \
  -f NC_045512.2.fasta \
  -o patient_joao_variants_lofreq.vcf \
  patient_joao_calibrated.bam
```

**Saída Esperada:**
- Arquivo VCF com variantes identificadas (SNPs e indels)

---

#### 2.4 Anotação de Variantes

**Ferramenta: SnpEff**

**Instalação:**
```bash
conda install -c bioconda snpeff
```

**Como Obter Banco de Dados SARS-CoV-2 para SnpEff:**
```bash
# Download do banco de dados SARS-CoV-2
# O SnpEff já vem com alguns bancos, mas pode ser necessário adicionar o SARS-CoV-2
# Verificar: snpeff databases | grep -i sars
```

**Como Utilizar:**
```bash
# Anotar variantes
snpEff \
  -v sarsCov2 \
  patient_joao_variants.vcf.gz \
  > patient_joao_variants_annotated.vcf
```

**Saída Esperada:**
- VCF anotado com impacto funcional das mutações (sinônima, não-sinônima, etc.)

**Processamento do VCF Anotado:**
```python
import pandas as pd
import vcf

# Ler VCF anotado
vcf_reader = vcf.Reader(open('patient_joao_variants_annotated.vcf', 'r'))

variants = []
for record in vcf_reader:
    variant = {
        'CHROM': record.CHROM,
        'POS': record.POS,
        'REF': record.REF,
        'ALT': str(record.ALT[0]),
        'QUAL': record.QUAL,
        'DP': record.INFO.get('DP', 0),
        'AF': record.INFO.get('AF', [0])[0] if record.INFO.get('AF') else 0
    }
    # Extrair anotações do SnpEff
    if 'ANN' in record.INFO:
        ann = record.INFO['ANN'][0].split('|')
        variant['GENE'] = ann[3]
        variant['EFFECT'] = ann[1]
        variant['IMPACT'] = ann[2]
    variants.append(variant)

df_variants = pd.DataFrame(variants)
df_variants.to_csv('patient_joao_variants_table.csv', index=False)
```

---

#### 2.5 Determinação de Linhagem

**Ferramenta: Pangolin**

**Instalação:**
```bash
# Via pip
pip install pangolin

# Ou via conda
conda install -c bioconda pangolin
```

**Como Utilizar:**
```bash
# Opção 1: Usar VCF
pangolin patient_joao_variants.vcf.gz --outfile patient_joao_lineage.csv

# Opção 2: Usar BAM (recomendado - mais preciso)
pangolin patient_joao_dedup.bam --outfile patient_joao_lineage.csv

# Opção 3: Usar FASTA (se tiver genoma consenso)
# Primeiro, gerar genoma consenso
samtools consensus -f fasta patient_joao_dedup.bam > patient_joao_consensus.fasta
pangolin patient_joao_consensus.fasta --outfile patient_joao_lineage.csv
```

**Saída Esperada:**
- CSV com linhagem identificada (ex: B.1.617.2 - Delta, BA.1 - Ômicron)
- Probabilidade e qualidade da classificação

**Web Interface (Alternativa):**
- https://pangolin.cog-uk.io/
- Upload do arquivo VCF ou FASTA consenso

---

### Saídas Finais da Etapa 2

1. **Arquivos de Alinhamento:**
   - BAM ordenado, indexado e sem duplicatas
   - Estatísticas de cobertura

2. **Variantes:**
   - Arquivo VCF com variantes identificadas
   - VCF anotado com impacto funcional

3. **Tabela de Variantes Anotadas:**
   - CSV com colunas: Posição, Gene, Mutação, Tipo, Impacto, Frequência

4. **Linhagem:**
   - Relatório Pangolin com linhagem identificada

5. **Interpretação:**
   - Mutações de interesse clínico (ex: N501Y, E484K, D614G)
   - Associação com variantes de preocupação (VOCs)

---

## Etapa 3: Extração e Engenharia de Features

### Objetivo
Transformar dados genômicos brutos em features numéricas utilizáveis por algoritmos de ML, normalizando para o formato TPM compatível com a matriz de treinamento.

### Dados de Entrada

**Dados do Paciente João:**
- Relatório Kraken2: `patient_joao_kraken2_report.txt`
- Relatório Bracken: `patient_joao_bracken_output.txt`
- Variantes identificadas: `patient_joao_variants_table.csv`
- Estatísticas de cobertura: `patient_joao_coverage.txt`

**Matriz de Treinamento:**
- `pivoted-virome-organisms-atleast10tpm-species-covid-TCC-pos.csv`

**Onde Obter:**
- Fornecido pelo professor (download do arquivo CSV)

---

### Passo a Passo Detalhado

#### 3.1 Extração de Features Básicas

**Processamento do Relatório Kraken2/Bracken:**

```python
import pandas as pd
import numpy as np

# Ler relatório Kraken2
kraken_report = pd.read_csv('patient_joao_kraken2_report.txt', sep='\t',
                           names=['percent', 'reads', 'tax_reads', 'rank', 'taxid', 'name'])

# Filtrar apenas espécies (rank == 'S')
species = kraken_report[kraken_report['rank'] == 'S'].copy()

# Calcular TPM (Transcripts Per Million)
total_reads = species['reads'].sum()
species['TPM'] = (species['reads'] / total_reads) * 1_000_000

# Criar dicionário de features (nome da espécie -> TPM)
features_dict = dict(zip(species['name'], species['TPM']))
```

**Alternativa usando Bracken (mais preciso):**
```python
# Ler relatório Bracken
bracken_report = pd.read_csv('patient_joao_bracken_output.txt', sep='\t')

# Bracken já fornece abundância estimada
features_dict = dict(zip(bracken_report['name'], bracken_report['new_est_reads']))
```

---

#### 3.2 Engenharia de Features Avançadas

**Features Adicionais que Podem ser Extraídas:**

```python
# 1. Diversidade Alfa (Shannon Index)
from scipy.stats import entropy

# Calcular Shannon Index
proportions = species['reads'] / species['reads'].sum()
shannon_index = entropy(proportions, base=2)

# 2. Razão Vírus/Bactérias
virus_reads = kraken_report[kraken_report['name'].str.contains('virus', case=False, na=False)]['reads'].sum()
bacteria_reads = kraken_report[kraken_report['rank'] == 'S']['reads'].sum() - virus_reads
virus_bacteria_ratio = virus_reads / bacteria_reads if bacteria_reads > 0 else 0

# 3. Presença/Ausência de Patógenos-Chave
key_pathogens = {
    'SARS-CoV-2': 'Severe acute respiratory syndrome coronavirus 2' in species['name'].values,
    'Influenza': 'Influenza' in ' '.join(species['name'].values),
    'RSV': 'Respiratory syncytial' in ' '.join(species['name'].values),
    'Streptococcus pneumoniae': 'Streptococcus pneumoniae' in species['name'].values
}

# 4. Número de Mutações no Gene Spike
variants_df = pd.read_csv('patient_joao_variants_table.csv')
spike_mutations = variants_df[variants_df['GENE'] == 'S'].shape[0]

# 5. Carga Viral Estimada (cobertura média do SARS-CoV-2)
coverage_df = pd.read_csv('patient_joao_coverage.txt', sep='\t', names=['chr', 'pos', 'depth'])
viral_load = coverage_df['depth'].mean()

# 6. Número Total de Organismos Detectados
total_organisms = species.shape[0]

# 7. Abundância Relativa do SARS-CoV-2
sars_cov2_tpm = features_dict.get('Severe acute respiratory syndrome coronavirus 2', 0)

# Compilar todas as features
advanced_features = {
    'shannon_index': shannon_index,
    'virus_bacteria_ratio': virus_bacteria_ratio,
    'spike_mutations': spike_mutations,
    'viral_load': viral_load,
    'total_organisms': total_organisms,
    'sars_cov2_abundance': sars_cov2_tpm,
    **key_pathogens  # Adiciona presença/ausência de patógenos
}
```

---

#### 3.3 Normalização e Transformação

**Alinhamento com a Matriz de Treinamento:**

```python
# Ler matriz de treinamento para obter colunas
train_matrix = pd.read_csv('pivoted-virome-organisms-atleast10tpm-species-covid-TCC-pos.csv')

# Obter lista de organismos (features) da matriz de treinamento
feature_columns = [col for col in train_matrix.columns if col not in ['sample_id', 'covid_status']]

# Criar vetor de features do paciente alinhado com a matriz de treinamento
patient_features = pd.DataFrame(index=[0])

for organism in feature_columns:
    # Procurar correspondência no dicionário de features
    # Pode ser necessário normalizar nomes (remover espaços, etc.)
    tpm_value = features_dict.get(organism, 0)
    patient_features[organism] = tpm_value

# Adicionar features avançadas
for feature_name, feature_value in advanced_features.items():
    patient_features[feature_name] = feature_value

# Aplicar transformação logarítmica (se necessário)
# TPM + 1 para evitar log(0)
patient_features_log = np.log1p(patient_features)

# Salvar vetor de features
patient_features.to_csv('patient_joao_features_vector.csv', index=False)
patient_features_log.to_csv('patient_joao_features_vector_log.csv', index=False)
```

**Justificativa das Transformações:**
- **TPM**: Normaliza por profundidade de sequenciamento e tamanho do genoma
- **Log-transform**: Estabiliza variância e reduz impacto de outliers
- **StandardScaler**: Pode ser aplicado antes do ML para normalizar escala

---

#### 3.4 Validação e Verificação

```python
# Verificar se todas as features da matriz de treinamento estão presentes
missing_features = set(feature_columns) - set(patient_features.columns)
if missing_features:
    print(f"Atenção: {len(missing_features)} features faltando")
    # Adicionar com valor 0
    for feat in missing_features:
        patient_features[feat] = 0

# Verificar distribuição
print("Estatísticas do vetor de features:")
print(patient_features.describe())

# Visualização
import matplotlib.pyplot as plt

# Top 20 organismos por TPM
top_orgs = patient_features[feature_columns].T.sort_values(by=0, ascending=False).head(20)
plt.figure(figsize=(10, 8))
plt.barh(range(len(top_orgs)), top_orgs[0].values)
plt.yticks(range(len(top_orgs)), top_orgs.index)
plt.xlabel('TPM')
plt.title('Top 20 Organismos - Paciente João')
plt.tight_layout()
plt.savefig('patient_features_top20.png')
```

---

### Saídas Finais da Etapa 3

1. **Vetor de Features:**
   - `patient_joao_features_vector.csv`: Features em TPM
   - `patient_joao_features_vector_log.csv`: Features com transformação log

2. **Features Avançadas:**
   - Diversidade alfa, razão vírus/bactérias, mutações Spike, etc.

3. **Validação:**
   - Verificação de alinhamento com matriz de treinamento
   - Estatísticas descritivas

4. **Visualizações:**
   - Gráficos de distribuição de features

---

## Etapa 4: Treinamento e Predição com Machine Learning

### Objetivo
Treinar um modelo supervisionado para predição de COVID-19 com base em dados metatranscriptômicos e aplicar na amostra do paciente.

### Dados de Entrada

**Matriz de Treinamento:**
- `pivoted-virome-organisms-atleast10tpm-species-covid-TCC-pos.csv`
  - Colunas: organismos (features) + coluna de rótulo (covid_status: True/False)
  - 100 amostras

**Vetor de Features do Paciente:**
- `patient_joao_features_vector.csv` ou `patient_joao_features_vector_log.csv`

---

### Passo a Passo Detalhado

#### 4.1 Exploração dos Dados (EDA)

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, roc_curve
import warnings
warnings.filterwarnings('ignore')

# Ler matriz de treinamento
df = pd.read_csv('pivoted-virome-organisms-atleast10tpm-species-covid-TCC-pos.csv')

# Verificar estrutura
print("Shape:", df.shape)
print("\nPrimeiras linhas:")
print(df.head())

# Verificar balanceamento das classes
print("\nDistribuição de classes:")
print(df['covid_status'].value_counts())
print(df['covid_status'].value_counts(normalize=True))

# Visualizar balanceamento
plt.figure(figsize=(8, 5))
df['covid_status'].value_counts().plot(kind='bar')
plt.title('Distribuição de Classes')
plt.xlabel('COVID-19 Status')
plt.ylabel('Número de Amostras')
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig('class_distribution.png')

# Verificar valores faltantes
print("\nValores faltantes:")
print(df.isnull().sum().sum())

# Análise de correlação entre features
feature_cols = [col for col in df.columns if col != 'covid_status']
correlation_matrix = df[feature_cols].corr()

plt.figure(figsize=(12, 10))
sns.heatmap(correlation_matrix, cmap='coolwarm', center=0, square=True, linewidths=0.5)
plt.title('Matriz de Correlação entre Features')
plt.tight_layout()
plt.savefig('correlation_matrix.png')
```

---

#### 4.2 Pré-processamento

```python
# Separar features e target
X = df[feature_cols].copy()
y = df['covid_status'].copy()

# Converter target para binário (se necessário)
if y.dtype == 'object':
    y = y.map({True: 1, False: 0, 'True': 1, 'False': 0})

# Aplicar transformação logarítmica (se ainda não aplicada)
X_log = np.log1p(X)

# Normalização (StandardScaler)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_log)
X_scaled = pd.DataFrame(X_scaled, columns=feature_cols)

# Divisão treino/teste
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y, test_size=0.2, random_state=42, stratify=y
)

print(f"Treino: {X_train.shape[0]} amostras")
print(f"Teste: {X_test.shape[0]} amostras")
print(f"Proporção de classes no treino: {y_train.value_counts(normalize=True)}")
print(f"Proporção de classes no teste: {y_test.value_counts(normalize=True)}")
```

---

#### 4.3 Treinamento de Modelos

**Opção 1: Random Forest**

```python
# Random Forest
rf_model = RandomForestClassifier(
    n_estimators=100,
    max_depth=10,
    min_samples_split=5,
    min_samples_leaf=2,
    random_state=42,
    n_jobs=-1
)

rf_model.fit(X_train, y_train)

# Predições
y_pred_rf = rf_model.predict(X_test)
y_pred_proba_rf = rf_model.predict_proba(X_test)[:, 1]

# Métricas
print("=== Random Forest ===")
print(classification_report(y_test, y_pred_rf))
print(f"ROC-AUC: {roc_auc_score(y_test, y_pred_proba_rf):.4f}")

# Importância das features
feature_importance = pd.DataFrame({
    'feature': feature_cols,
    'importance': rf_model.feature_importances_
}).sort_values('importance', ascending=False)

print("\nTop 10 Features Mais Importantes:")
print(feature_importance.head(10))

# Visualizar importância
plt.figure(figsize=(10, 8))
top_features = feature_importance.head(20)
plt.barh(range(len(top_features)), top_features['importance'].values)
plt.yticks(range(len(top_features)), top_features['feature'].values)
plt.xlabel('Importância')
plt.title('Top 20 Features - Random Forest')
plt.tight_layout()
plt.savefig('feature_importance_rf.png')
```

**Opção 2: XGBoost**

```python
# XGBoost
xgb_model = XGBClassifier(
    n_estimators=100,
    max_depth=5,
    learning_rate=0.1,
    random_state=42,
    eval_metric='logloss'
)

xgb_model.fit(X_train, y_train)

# Predições
y_pred_xgb = xgb_model.predict(X_test)
y_pred_proba_xgb = xgb_model.predict_proba(X_test)[:, 1]

# Métricas
print("=== XGBoost ===")
print(classification_report(y_test, y_pred_xgb))
print(f"ROC-AUC: {roc_auc_score(y_test, y_pred_proba_xgb):.4f}")

# Importância das features
feature_importance_xgb = pd.DataFrame({
    'feature': feature_cols,
    'importance': xgb_model.feature_importances_
}).sort_values('importance', ascending=False)

print("\nTop 10 Features Mais Importantes:")
print(feature_importance_xgb.head(10))
```

**Opção 3: Logistic Regression**

```python
# Logistic Regression
lr_model = LogisticRegression(
    max_iter=1000,
    random_state=42,
    C=1.0
)

lr_model.fit(X_train, y_train)

# Predições
y_pred_lr = lr_model.predict(X_test)
y_pred_proba_lr = lr_model.predict_proba(X_test)[:, 1]

# Métricas
print("=== Logistic Regression ===")
print(classification_report(y_test, y_pred_lr))
print(f"ROC-AUC: {roc_auc_score(y_test, y_pred_proba_lr):.4f}")

# Coeficientes
coef_df = pd.DataFrame({
    'feature': feature_cols,
    'coefficient': lr_model.coef_[0]
}).sort_values('coefficient', key=abs, ascending=False)

print("\nTop 10 Features por Coeficiente:")
print(coef_df.head(10))
```

---

#### 4.4 Avaliação e Seleção do Melhor Modelo

```python
# Comparar modelos
models = {
    'Random Forest': (y_pred_rf, y_pred_proba_rf),
    'XGBoost': (y_pred_xgb, y_pred_proba_xgb),
    'Logistic Regression': (y_pred_lr, y_pred_proba_lr)
}

results = []
for name, (y_pred, y_pred_proba) in models.items():
    results.append({
        'Model': name,
        'Accuracy': (y_pred == y_test).mean(),
        'ROC-AUC': roc_auc_score(y_test, y_pred_proba),
        'Precision': classification_report(y_test, y_pred, output_dict=True)['1']['precision'],
        'Recall': classification_report(y_test, y_pred, output_dict=True)['1']['recall'],
        'F1-Score': classification_report(y_test, y_pred, output_dict=True)['1']['f1-score']
    })

results_df = pd.DataFrame(results)
print("\nComparação de Modelos:")
print(results_df)

# Visualizar comparação
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Métricas de classificação
metrics = ['Accuracy', 'Precision', 'Recall', 'F1-Score']
x = np.arange(len(metrics))
width = 0.25

for i, model in enumerate(results_df['Model']):
    values = [results_df.loc[i, m] for m in metrics]
    axes[0].bar(x + i*width, values, width, label=model)

axes[0].set_xlabel('Métricas')
axes[0].set_ylabel('Score')
axes[0].set_title('Comparação de Métricas de Classificação')
axes[0].set_xticks(x + width)
axes[0].set_xticklabels(metrics)
axes[0].legend()
axes[0].grid(axis='y', alpha=0.3)

# ROC-AUC
axes[1].bar(results_df['Model'], results_df['ROC-AUC'])
axes[1].set_ylabel('ROC-AUC Score')
axes[1].set_title('Comparação de ROC-AUC')
axes[1].grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('model_comparison.png')

# Selecionar melhor modelo (baseado em ROC-AUC)
best_model_name = results_df.loc[results_df['ROC-AUC'].idxmax(), 'Model']
print(f"\nMelhor modelo: {best_model_name}")

# Salvar melhor modelo
if best_model_name == 'Random Forest':
    best_model = rf_model
elif best_model_name == 'XGBoost':
    best_model = xgb_model
else:
    best_model = lr_model

import joblib
joblib.dump(best_model, 'best_model.pkl')
joblib.dump(scaler, 'scaler.pkl')
```

**Visualização de Métricas:**

```python
# Matriz de Confusão
from sklearn.metrics import ConfusionMatrixDisplay

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for idx, (name, (y_pred, _)) in enumerate(models.items()):
    ConfusionMatrixDisplay.from_predictions(y_test, y_pred, ax=axes[idx])
    axes[idx].set_title(f'{name} - Matriz de Confusão')

plt.tight_layout()
plt.savefig('confusion_matrices.png')

# Curvas ROC
fig, ax = plt.subplots(figsize=(8, 6))

for name, (_, y_pred_proba) in models.items():
    fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
    auc = roc_auc_score(y_test, y_pred_proba)
    ax.plot(fpr, tpr, label=f'{name} (AUC = {auc:.3f})')

ax.plot([0, 1], [0, 1], 'k--', label='Random')
ax.set_xlabel('Taxa de Falsos Positivos')
ax.set_ylabel('Taxa de Verdadeiros Positivos')
ax.set_title('Curvas ROC')
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('roc_curves.png')
```

---

#### 4.5 Predição da Amostra do Paciente

```python
# Ler vetor de features do paciente
patient_features = pd.read_csv('patient_joao_features_vector.csv')

# Garantir mesma ordem de colunas
patient_X = patient_features[feature_cols].copy()

# Aplicar mesma transformação
patient_X_log = np.log1p(patient_X)
patient_X_scaled = scaler.transform(patient_X_log)
patient_X_scaled = pd.DataFrame(patient_X_scaled, columns=feature_cols)

# Predição
patient_pred = best_model.predict(patient_X_scaled)[0]
patient_pred_proba = best_model.predict_proba(patient_X_scaled)[0]

print("=== Predição para Paciente João ===")
print(f"Predição: {'POSITIVO para COVID-19' if patient_pred == 1 else 'NEGATIVO para COVID-19'}")
print(f"Probabilidade de COVID-19: {patient_pred_proba[1]:.4f} ({patient_pred_proba[1]*100:.2f}%)")
print(f"Probabilidade de NÃO COVID-19: {patient_pred_proba[0]:.4f} ({patient_pred_proba[0]*100:.2f}%)")

# Salvar resultado
result = pd.DataFrame({
    'patient': ['João'],
    'prediction': [patient_pred],
    'probability_covid': [patient_pred_proba[1]],
    'probability_no_covid': [patient_pred_proba[0]]
})
result.to_csv('patient_joao_prediction.csv', index=False)
```

**Interpretação com SHAP (opcional, mas recomendado):**

```python
# Instalar: pip install shap
import shap

# Criar explainer
explainer = shap.TreeExplainer(best_model)
shap_values = explainer.shap_values(patient_X_scaled)

# Visualizar importância para esta predição específica
shap.summary_plot(shap_values, patient_X_scaled, plot_type="bar", show=False)
plt.savefig('shap_summary.png', bbox_inches='tight')

# Waterfall plot para a predição do paciente
shap.waterfall_plot(explainer(patient_X_scaled)[0], show=False)
plt.savefig('shap_waterfall_patient.png', bbox_inches='tight')
```

---

### Saídas Finais da Etapa 4

1. **Modelo Treinado:**
   - `best_model.pkl`: Modelo salvo
   - `scaler.pkl`: Scaler salvo

2. **Avaliação:**
   - Relatório de classificação (precision, recall, F1-score)
   - Matriz de confusão
   - Curva ROC
   - Comparação de modelos

3. **Predição do Paciente:**
   - `patient_joao_prediction.csv`: Resultado da predição
   - Probabilidade de COVID-19

4. **Interpretabilidade:**
   - Importância de features
   - Gráficos SHAP (se aplicado)

5. **Visualizações:**
   - Distribuição de classes
   - Matriz de correlação
   - Importância de features
   - Matriz de confusão
   - Curvas ROC

---

## Recursos e Referências

### Bancos de Dados e Referências

1. **Genoma de Referência SARS-CoV-2:**
   - https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2
   - Genoma Wuhan (padrão)

2. **Genotipagem SARS-CoV-2:**
   - https://pangolin.cog-uk.io/
   - Ferramenta web e linha de comando

3. **Análise de Resistência a Antibióticos:**
   - https://card.mcmaster.ca/analyze/rgi
   - Comprehensive Antibiotic Resistance Database (CARD)

4. **Genoma Humano (HG38):**
   - https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
   - Para remoção de reads do hospedeiro

5. **Bancos Kraken2 Pré-construídos:**
   - https://genome-idx.s3.amazonaws.com/kraken/
   - Bancos de vírus e bactérias

### Ferramentas Principais

| Ferramenta | Uso | Instalação |
|-----------|-----|------------|
| FastQC | Controle de qualidade | `conda install -c bioconda fastqc` |
| MultiQC | Consolidação de relatórios | `conda install -c bioconda multiqc` |
| Trimmomatic | Trimming de reads | `conda install -c bioconda trimmomatic` |
| fastp | Trimming rápido | `conda install -c bioconda fastp` |
| Bowtie2 | Alinhamento | `conda install -c bioconda bowtie2` |
| BWA | Alinhamento | `conda install -c bioconda bwa` |
| Kraken2 | Classificação taxonômica | `conda install -c bioconda kraken2` |
| Bracken | Estimativa de abundância | `conda install -c bioconda bracken` |
| samtools | Manipulação de BAM | `conda install -c bioconda samtools` |
| bcftools | Chamada de variantes | `conda install -c bioconda bcftools` |
| FreeBayes | Chamada de variantes | `conda install -c bioconda freebayes` |
| SnpEff | Anotação de variantes | `conda install -c bioconda snpeff` |
| Pangolin | Determinação de linhagem | `pip install pangolin` |
| Krona | Visualização taxonômica | `conda install -c bioconda krona` |

### Bibliotecas Python

```bash
pip install pandas numpy matplotlib seaborn scikit-learn xgboost shap scipy
```

### Estrutura de Arquivos Recomendada

```
tcc-metatranscriptomica/
├── data/
│   ├── raw/
│   │   ├── patient_joao_VIROMA_S21_R1_001.fastq.gz
│   │   └── patient_joao_VIROMA_S21_R2_001.fastq.gz
│   ├── processed/
│   │   ├── patient_joao_R1_trimmed.fastq.gz
│   │   └── patient_joao_R2_trimmed.fastq.gz
│   ├── references/
│   │   ├── NC_045512.2.fasta
│   │   └── hg38.fa
│   └── training/
│       └── pivoted-virome-organisms-atleast10tpm-species-covid-TCC-pos.csv
├── notebooks/
│   ├── 01_metatranscriptomica_pipeline.ipynb
│   ├── 02_variant_calling_pipeline.ipynb
│   ├── 03_feature_engineering.ipynb
│   └── 04_ml_training_prediction.ipynb
├── results/
│   ├── qc_reports/
│   ├── kraken2_reports/
│   ├── variants/
│   ├── features/
│   └── ml_results/
└── guia-execucao-tcc.md
```

---

## Checklist Final

Antes da entrega, verificar:

- [ ] Todos os notebooks executam sem erros
- [ ] Relatórios de qualidade gerados e interpretados
- [ ] Classificação taxonômica completa (Kraken2/Bracken)
- [ ] Variantes identificadas e anotadas
- [ ] Linhagem determinada (Pangolin)
- [ ] Vetor de features alinhado com matriz de treinamento
- [ ] Modelo treinado e avaliado com métricas adequadas
- [ ] Predição do paciente realizada
- [ ] Todas as questões da apresentação respondidas
- [ ] Visualizações claras e informativas
- [ ] Código comentado e documentado

---

**Boa sorte com o TCC! 🧬🔬🤖**

