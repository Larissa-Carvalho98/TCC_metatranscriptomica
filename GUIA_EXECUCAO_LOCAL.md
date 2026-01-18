# Guia de Execução Local - TCC Metatranscriptômica

Este guia detalha o processo completo de execução do pipeline localmente, desde a configuração do ambiente até a execução de todos os notebooks.

## Índice

1. [Pré-requisitos](#pré-requisitos)
2. [Configuração do Ambiente](#configuração-do-ambiente)
3. [Download de Dados e Referências](#download-de-dados-e-referências)
4. [Configuração de Bancos de Dados](#configuração-de-bancos-de-dados)
5. [Execução dos Notebooks](#execução-dos-notebooks)
6. [Troubleshooting](#troubleshooting)

---

## Pré-requisitos

### Software Necessário

1. **Python 3.9 ou 3.10** instalado
2. **Conda** (Miniconda ou Anaconda) - Recomendado para ferramentas bioinformáticas
   - Download: https://docs.conda.io/en/latest/miniconda.html
3. **Git** (opcional, se usar controle de versão)
4. **Jupyter Notebook/Lab** (será instalado via conda ou pip)

### Espaço em Disco

Recomenda-se pelo menos **200GB** de espaço livre para:
- Dados FASTQ (~2-5GB)
- Genomas de referência (~3-4GB)
- Bancos Kraken2 (~20-100GB dependendo do banco)
- Resultados intermediários e finais (~10-50GB)

---

## Configuração do Ambiente

### Passo 1: Navegar para o Diretório do Projeto

```bash
cd /Users/larissa/Desktop/TCC_metatrascriptomica
```

### Passo 2: Criar Ambiente Conda

**Opção A: Script Automatizado (Recomendado)**

```bash
# Dar permissão de execução
chmod +x setup_environment.sh

# Executar script
bash setup_environment.sh conda
```

**Opção B: Comandos Manuais**

```bash
# Criar ambiente a partir do environment.yml
conda env create -f environment.yml

# Ativar ambiente
conda activate tcc-metatranscriptomica

# Verificar instalação
conda list | grep -E "(fastqc|kraken2|samtools|pandas)"
```

### Passo 3: Verificar Instalação

```bash
# Verificar Python
python --version  # Deve ser 3.9 ou 3.10

# Verificar ferramentas bioinformáticas
fastqc --version
kraken2 --version
samtools --version
```

### Passo 4: Configurar Jupyter

```bash
# Instalar kernel do ambiente no Jupyter
python -m ipykernel install --user --name=tcc-metatranscriptomica --display-name "TCC Metatranscriptomica"

# Iniciar Jupyter
jupyter notebook
# ou
jupyter lab
```

**No Jupyter**: Ao criar um novo notebook ou abrir um existente, certifique-se de selecionar o kernel "TCC Metatranscriptomica" no menu Kernel → Change Kernel.

---

## Download de Dados e Referências

### Passo 1: Criar Estrutura de Diretórios

Os diretórios já foram criados automaticamente, mas você pode verificar:

```bash
ls -la data/
# Deve mostrar: raw/, processed/, references/, training/
```

### Passo 2: Download dos Dados FASTQ do Paciente

```bash
# Navegar para diretório de dados brutos
cd data/raw/

# Download dos arquivos FASTQ
wget https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R1_001.fastq.gz

wget https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R2_001.fastq.gz

# Verificar download
ls -lh *.fastq.gz
```

**Alternativa com curl** (se wget não estiver disponível):

```bash
curl -O https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R1_001.fastq.gz
curl -O https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R2_001.fastq.gz
```

### Passo 3: Download do Genoma de Referência SARS-CoV-2

```bash
# Navegar para diretório de referências
cd ../references/

# Download genoma Wuhan (NC_045512.2)
wget -O NC_045512.2.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1798174254"

# Verificar arquivo
head -5 NC_045512.2.fasta
```

**Alternativa usando efetch (se instalado):**

```bash
efetch -db nuccore -id NC_045512.2 -format fasta > NC_045512.2.fasta
```

### Passo 4: Download do Genoma Humano (HG38) - Opcional

Necessário apenas para remoção de reads do hospedeiro:

```bash
# Download genoma humano (pode ser grande: ~3GB)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Descompactar
gunzip hg38.fa.gz

# Verificar
ls -lh hg38.fa
```

**Nota**: O download do genoma humano pode demorar. Se preferir, pode usar uma versão já indexada ou pular esta etapa se não for fazer remoção de hospedeiro.

### Passo 5: Download da Matriz de Treinamento

```bash
# Navegar para diretório de treinamento
cd ../training/

# Download da matriz (ajustar URL se necessário)
# Se o arquivo estiver em um repositório ou servidor, fazer download manual via browser
# e mover para: data/training/pivoted-virome-organisms-atleast10tpm-species-covid-TCC-pos.csv
```

---

## Configuração de Bancos de Dados

### Banco Kraken2

O Kraken2 precisa de um banco de dados para classificação taxonômica. Existem três opções:

#### Opção 1: Banco de Vírus (Recomendado para TCC - Menor e mais rápido)

```bash
# Criar diretório para banco
mkdir -p ~/kraken2_db/viral_db

# Download do banco de vírus (pré-construído, ~4-8GB)
# Link: https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz
cd ~/kraken2_db/
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz
tar -xzf k2_viral_20221209.tar.gz

# Verificar como os arquivos foram extraídos e organizá-los
if [ -d "k2_viral_20221209" ]; then
    # Se foi extraído em subdiretório, mover conteúdo
    mv k2_viral_20221209/* viral_db/
    rmdir k2_viral_20221209
else
    # Se foi extraído diretamente, mover arquivos do banco para viral_db
    find . -maxdepth 1 -type f \( -name "*.k2d" -o -name "*.map" -o -name "*.kmer_distrib" -o -name "*.txt" -o -name "opts.k2d" \) ! -name "*.tar.gz" -exec mv {} viral_db/ \;
fi

# Verificar que arquivos essenciais estão presentes
cd viral_db
ls -lh *.k2d *.map
```

#### Opção 2: Construir Banco Manualmente

```bash
# Criar banco de vírus
kraken2-build --download-library viral --db ~/kraken2_db/viral_db
kraken2-build --build --db ~/kraken2_db/viral_db
```

#### Opção 3: Banco Completo (Vírus + Bactérias) - Maior e mais lento

```bash
# Download do banco completo (~100GB)
# Link: https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20221209.tar.gz
# Este banco leva muito tempo para download e processamento
```

**Nota**: Lembre-se do caminho do banco para usar nos notebooks: `~/kraken2_db/viral_db` ou caminho absoluto.

### Indexação do Genoma Humano (para Bowtie2)

Se você baixou o genoma humano, precisa criar o índice Bowtie2:

```bash
# Navegar para diretório de referências
cd /Users/larissa/Desktop/TCC_metatrascriptomica/data/references/

# Criar índice Bowtie2 (pode demorar 30-60 minutos)
bowtie2-build hg38.fa hg38_index

# Verificar arquivos de índice criados
ls -lh hg38_index*
```

---

## Execução dos Notebooks

### Ordem de Execução

Os notebooks devem ser executados **sequencialmente**, pois cada um depende dos resultados do anterior:

1. `01_metatranscriptomica_pipeline.ipynb`
2. `02_variant_calling_pipeline.ipynb`
3. `03_feature_engineering.ipynb`
4. `04_ml_training_prediction.ipynb`

### Passo 1: Iniciar Jupyter

```bash
# Ativar ambiente conda (se ainda não estiver ativo)
conda activate tcc-metatranscriptomica

# Navegar para diretório do projeto
cd /Users/larissa/Desktop/TCC_metatrascriptomica

# Iniciar Jupyter
jupyter notebook
# ou
jupyter lab
```

Isso abrirá o navegador automaticamente. Navegue até a pasta `notebooks/`.

### Passo 2: Verificar Ambiente (Opcional mas Recomendado)

1. Abrir `00_VERIFICACAO_AMBIENTE.ipynb` (se disponível)
2. Executar a célula para verificar se o ambiente está configurado corretamente
3. Isso ajuda a diagnosticar problemas de importação antes de executar os pipelines

### Passo 3: Executar Notebook 1 - Pipeline Metatranscriptômica

1. Abrir `01_metatranscriptomica_pipeline.ipynb`
2. **IMPORTANTE**: Executar a **primeira célula de código** (configuração inicial) antes de qualquer importação
3. **Verificar configurações:**
   - Verificar que o kernel está selecionado: "TCC Metatranscriptomica"
   - Na célula de configuração, verificar se o diretório `scripts` foi encontrado (deve mostrar ✅)
   - Verificar caminhos dos arquivos FASTQ

3. **Ajustar variáveis:**
   - Descomentar e ajustar o caminho do banco Kraken2:
     ```python
     kraken_db = "/Users/larissa/kraken2_db/viral_db"  # Ajustar seu caminho
     ```

4. **Executar células sequencialmente:**
   - Executar todas as células (Cell → Run All) ou uma por uma (Shift+Enter)
   - Verificar erros e avisos em cada etapa

5. **Verificar saídas:**
   - Relatórios FastQC em `results/qc_reports/`
   - Relatórios Kraken2/Bracken em `results/kraken2_reports/`

### Passo 4: Executar Notebook 2 - Chamada de Variantes

1. Abrir `02_variant_calling_pipeline.ipynb`

2. **Verificar dependências:**
   - Certificar-se que o Notebook 1 foi executado completamente
   - Arquivos necessários:
     - `data/processed/patient_joao_nonhuman_R1.fastq.gz`
     - `data/processed/patient_joao_nonhuman_R2.fastq.gz`
     - `data/references/NC_045512.2.fasta`

3. **Executar células:**
   - BWA criará índice automaticamente se necessário
   - Processamento de BAM pode demorar (30min-2h dependendo do tamanho)
   - Chamada de variantes também pode demorar

4. **Verificar saídas:**
   - BAM ordenado: `results/variants/patient_joao_sorted.bam`
   - VCF com variantes: `results/variants/patient_joao_variants.vcf.gz`
   - VCF anotado: `results/variants/patient_joao_variants_annotated.vcf`
   - Linhagem: `results/variants/patient_joao_lineage.csv`

### Passo 5: Executar Notebook 3 - Engenharia de Features

1. Abrir `03_feature_engineering.ipynb`

2. **Verificar dependências:**
   - Relatório Kraken2/Bracken do Notebook 1
   - Tabela de variantes do Notebook 2

3. **Ajustar caminhos se necessário:**
   - Verificar caminhos dos relatórios
   - Verificar caminho da matriz de treinamento

4. **Executar células:**
   - Extração de TPM
   - Cálculo de features avançadas
   - Criação do vetor de features alinhado

5. **Verificar saídas:**
   - `results/features/patient_joao_features_vector.csv`
   - `results/features/patient_joao_features_vector_log.csv`

### Passo 6: Executar Notebook 4 - Machine Learning

1. Abrir `04_ml_training_prediction.ipynb`

2. **Verificar dependências:**
   - Matriz de treinamento em `data/training/`
   - Vetor de features do paciente do Notebook 3

3. **Executar células:**
   - EDA (Exploração de Dados)
   - Treinamento de modelos (pode demorar alguns minutos)
   - Avaliação e seleção do melhor modelo
   - Predição do paciente

4. **Verificar saídas:**
   - Modelo salvo: `results/ml_results/best_model.pkl`
   - Scaler salvo: `results/ml_results/scaler.pkl`
   - Resultado da predição: console/notebook

---

## Comandos Úteis para Verificação

### Verificar Tamanho dos Arquivos

```bash
# Verificar dados brutos
du -sh data/raw/*.fastq.gz

# Verificar resultados
du -sh results/*/

# Verificar tamanho total do projeto
du -sh .
```

### Verificar Processos em Execução

```bash
# Ver processos Python/Jupyter
ps aux | grep -E "(python|jupyter)"

# Ver processos bioinformáticos
ps aux | grep -E "(fastqc|kraken2|samtools|bwa)"
```

### Monitorar Uso de Memória/Disco

```bash
# Uso de memória
free -h  # Linux
vm_stat  # macOS

# Espaço em disco
df -h
```

---

## Troubleshooting

### Problema: "ModuleNotFoundError: No module named 'scripts'"

**Solução:**

1. **Verificar que a célula de configuração inicial foi executada primeiro:**
   - A primeira célula de código no notebook deve adicionar o projeto ao `sys.path`
   - Execute essa célula ANTES de tentar importar qualquer coisa dos scripts

2. **Verificar diretório de trabalho:**
   ```python
   # No notebook, executar:
   import os
   from pathlib import Path
   
   # Verificar diretório atual
   print(f"Diretório atual: {os.getcwd()}")
   print(f"Diretório atual (Path): {Path().resolve()}")
   
   # Verificar se scripts existe
   project_root = Path('/Users/larissa/Desktop/TCC_metatrascriptomica')
   print(f"Scripts existe: {(project_root / 'scripts').exists()}")
   ```

3. **Forçar caminho correto manualmente:**
   ```python
   # Adicionar no início do notebook (ANTES de qualquer import)
   import sys
   import os
   from pathlib import Path
   
   # Definir raiz do projeto explicitamente
   project_root = Path('/Users/larissa/Desktop/TCC_metatrascriptomica')
   
   # Adicionar ao path
   if str(project_root) not in sys.path:
       sys.path.insert(0, str(project_root))
   
   # Mudar para diretório do notebook se necessário
   os.chdir(project_root / 'notebooks')
   
   # Verificar
   import scripts  # Deve funcionar agora
   print(f"✅ Scripts importado com sucesso!")
   ```

4. **Usar notebook de verificação:**
   - Execute o notebook `00_VERIFICACAO_AMBIENTE.ipynb` primeiro
   - Ele faz todos os diagnósticos automaticamente

### Problema: "Command not found: fastqc" (ou outras ferramentas)

**Solução:**
- Verificar que o ambiente conda está ativo: `conda activate tcc-metatranscriptomica`
- Reinstalar ferramenta: `conda install -c bioconda fastqc`
- Verificar PATH: `echo $PATH`

### Problema: "Out of memory" durante processamento

**Soluções:**
- Fechar outros programas
- Reduzir número de threads nas funções
- Processar dados em lotes menores
- Aumentar swap do sistema (Linux)

### Problema: Arquivos FASTQ não encontrados

**Solução:**
- Verificar caminhos dos arquivos nas células do notebook
- Verificar que os arquivos foram baixados corretamente
- Usar caminhos absolutos temporariamente para debug

### Problema: Banco Kraken2 não encontrado

**Solução:**
- Verificar caminho do banco (pode ser relativo ou absoluto)
- Verificar que o banco foi baixado/construído corretamente
- Testar acesso: `ls ~/kraken2_db/viral_db/`

### Problema: Jupyter não abre

**Soluções:**
```bash
# Reinstalar Jupyter
conda install -c conda-forge jupyter

# Ou via pip
pip install jupyter

# Iniciar em porta específica
jupyter notebook --port 8888
```

### Problema: Permissão negada em scripts

**Solução:**
```bash
chmod +x setup_environment.sh
chmod -R 755 scripts/
```

---

## Dicas de Otimização

### Para Processamento Mais Rápido

1. **Aumentar threads:**
   - Ajustar parâmetro `threads=8` (ou maior) nas funções
   - Depende do número de CPUs disponíveis

2. **Processamento Paralelo:**
   - Algumas ferramentas (Kraken2, BWA) já usam múltiplas threads
   - Verificar disponibilidade com `nproc` (Linux) ou `sysctl -n hw.ncpu` (macOS)

3. **Usar SSDs:**
   - Processamento de BAM e grandes arquivos é mais rápido em SSDs

### Para Economizar Espaço

1. **Remover arquivos intermediários:**
   ```bash
   # Após completar pipeline, remover SAM temporários
   rm results/variants/*.sam
   ```

2. **Comprimir resultados:**
   ```bash
   # Comprimir relatórios grandes
   gzip results/kraken2_reports/*.txt
   ```

3. **Usar .gitignore:**
   - Arquivos grandes já estão no .gitignore para não versionar

---

## Checklist de Execução

Antes de começar, verificar:

- [ ] Ambiente conda criado e ativado
- [ ] Dados FASTQ baixados (`data/raw/`)
- [ ] Genoma de referência SARS-CoV-2 baixado (`data/references/`)
- [ ] Genoma humano baixado e indexado (se necessário)
- [ ] Banco Kraken2 configurado
- [ ] Matriz de treinamento disponível (`data/training/`)
- [ ] Jupyter instalado e kernel configurado
- [ ] Espaço em disco suficiente (>50GB livre)

Durante execução:

- [ ] Notebook 1: FastQC, Trimming, Kraken2 executados
- [ ] Notebook 2: Alinhamento, Variantes, Linhagem executados
- [ ] Notebook 3: Features extraídas e vetor criado
- [ ] Notebook 4: Modelo treinado e predição realizada

---

## Próximos Passos

Após completar a execução local:

1. **Verificar resultados:**
   - Revisar todas as saídas em `results/`
   - Verificar logs de erros se houver

2. **Preparar para Colab:**
   - Ver [COLAB_MIGRATION.md](COLAB_MIGRATION.md)
   - Notebooks já estão prontos para migração

3. **Análise e Interpretação:**
   - Usar resultados para responder questões do TCC
   - Gerar visualizações adicionais se necessário

---

**Boa execução! 🧬🔬**

Para dúvidas ou problemas, consulte:
- [README.md](README.md) - Visão geral do projeto
- [COLAB_MIGRATION.md](COLAB_MIGRATION.md) - Guia para Colab
- [guia-execucao-tcc.md](guia-execucao-tcc.md) - Guia técnico detalhado
