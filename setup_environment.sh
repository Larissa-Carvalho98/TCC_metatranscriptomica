#!/bin/bash

# Script de setup automatizado para TCC Metatranscriptômica
# Uso: bash setup_environment.sh [conda|venv]

set -e  # Parar em caso de erro

ENV_TYPE=${1:-conda}  # Padrão: conda

echo "=========================================="
echo "Setup do Ambiente TCC Metatranscriptômica"
echo "=========================================="

if [ "$ENV_TYPE" == "conda" ]; then
    echo "Criando ambiente conda..."
    
    # Verificar se conda está instalado
    if ! command -v conda &> /dev/null; then
        echo "ERRO: Conda não está instalado. Por favor, instale o Miniconda ou Anaconda."
        exit 1
    fi
    
    # Criar ambiente
    conda env create -f environment.yml
    
    echo ""
    echo "Ambiente conda criado com sucesso!"
    echo "Para ativar o ambiente, execute:"
    echo "  conda activate tcc-metatranscriptomica"
    echo ""
    echo "Para instalar Jupyter kernel no ambiente:"
    echo "  conda activate tcc-metatranscriptomica"
    echo "  python -m ipykernel install --user --name=tcc-metatranscriptomica"
    
elif [ "$ENV_TYPE" == "venv" ]; then
    echo "Criando ambiente virtual Python..."
    
    # Verificar Python
    if ! command -v python3.10 &> /dev/null && ! command -v python3.9 &> /dev/null; then
        echo "ERRO: Python 3.9 ou 3.10 não encontrado."
        exit 1
    fi
    
    # Usar python3.10 se disponível, senão python3.9
    PYTHON_CMD=$(command -v python3.10 || command -v python3.9)
    echo "Usando: $PYTHON_CMD"
    
    # Criar venv
    $PYTHON_CMD -m venv venv
    
    # Ativar e instalar dependências
    source venv/bin/activate
    pip install --upgrade pip
    pip install -r requirements.txt
    
    echo ""
    echo "Ambiente virtual criado com sucesso!"
    echo "Para ativar o ambiente, execute:"
    echo "  source venv/bin/activate"
    echo ""
    echo "Para instalar Jupyter kernel:"
    echo "  source venv/bin/activate"
    echo "  pip install jupyter ipykernel"
    echo "  python -m ipykernel install --user --name=tcc-metatranscriptomica"
    
else
    echo "ERRO: Tipo de ambiente inválido. Use 'conda' ou 'venv'"
    exit 1
fi

echo ""
echo "=========================================="
echo "Setup concluído!"
echo "=========================================="
echo ""
echo "NOTA: Ferramentas bioinformáticas (fastqc, kraken2, etc.)"
echo "      podem precisar ser instaladas separadamente se usar venv."
echo "      Recomenda-se usar conda para instalar essas ferramentas."
