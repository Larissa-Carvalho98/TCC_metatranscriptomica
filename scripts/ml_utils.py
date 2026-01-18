"""
Módulo utilitário para Machine Learning: treinamento, avaliação e predição.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Tuple, Optional, Dict
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score
)

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except ImportError:
    XGBOOST_AVAILABLE = False
    print("Aviso: XGBoost não disponível. Instale com: pip install xgboost")


def load_training_data(csv_file: str, target_column: str = 'covid_status') -> Tuple[pd.DataFrame, pd.Series]:
    """
    Carrega matriz de treinamento.
    
    Args:
        csv_file: Caminho para CSV com dados de treinamento
        target_column: Nome da coluna target
    
    Returns:
        Tuple: (X, y) features e target
    """
    if not Path(csv_file).exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {csv_file}")
    
    df = pd.read_csv(csv_file)
    
    if target_column not in df.columns:
        raise ValueError(f"Coluna target '{target_column}' não encontrada")
    
    X = df.drop(columns=[target_column])
    y = df[target_column]
    
    # Converter target para binário se necessário
    if y.dtype == 'object' or y.dtype == 'bool':
        y = y.map({True: 1, False: 0, 'True': 1, 'False': 0, 'true': 1, 'false': 0})
    
    return X, y

import numpy as np
import pandas as pd
from typing import Tuple
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer

def prepare_data(
    X: pd.DataFrame,
    y: pd.Series,
    test_size: float = 0.2,
    random_state: int = 42,
    apply_log: bool = True,
    apply_scaling: bool = True,
    variance_threshold: float = 0.0
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series, StandardScaler]:

    # 1️⃣ Garantir tipo numérico
    X = X.apply(pd.to_numeric, errors="coerce")

    # 2️⃣ Remover colunas 100% NaN
    X = X.dropna(axis=1, how="all")
    print(f"🔹 Após drop all-NaN: {X.shape[1]} features")

    # 3️⃣ Log transform
    if apply_log:
        X = np.log1p(X)

    # 4️⃣ Variance Threshold
    vt = VarianceThreshold(threshold=variance_threshold)
    X_vt = vt.fit_transform(X)

    kept_features = X.columns[vt.get_support()]
    X = pd.DataFrame(X_vt, columns=kept_features, index=X.index)

    print(f"🔹 Após VarianceThreshold: {X.shape[1]} features")

    # 5️⃣ Split
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=test_size,
        random_state=random_state,
        stratify=y
    )

    # 6️⃣ Remover colunas all-NaN no treino
    all_nan_cols = X_train.columns[X_train.isna().all()]
    if len(all_nan_cols) > 0:
        print(f"⚠️ Removendo {len(all_nan_cols)} colunas all-NaN no treino")
        X_train = X_train.drop(columns=all_nan_cols)
        X_test = X_test.drop(columns=all_nan_cols)

    # 7️⃣ Imputação
    imputer = SimpleImputer(strategy="median")
    X_train = pd.DataFrame(
        imputer.fit_transform(X_train),
        columns=X_train.columns,
        index=X_train.index
    )
    X_test = pd.DataFrame(
        imputer.transform(X_test),
        columns=X_test.columns,
        index=X_test.index
    )

    # 8️⃣ Scaling
    scaler = None
    if apply_scaling:
        scaler = StandardScaler()
        X_train = pd.DataFrame(
            scaler.fit_transform(X_train),
            columns=X_train.columns,
            index=X_train.index
        )
        X_test = pd.DataFrame(
            scaler.transform(X_test),
            columns=X_test.columns,
            index=X_test.index
        )

        # 🔑 salvar features usadas no scaler
        scaler.feature_names_ = X_train.columns.tolist()

    # 9️⃣ Checks finais
    assert not np.isnan(X_train.values).any()
    assert not np.isnan(X_test.values).any()

    print(f"✅ Dataset final: {X_train.shape[1]} features")

    return X_train, X_test, y_train, y_test, scaler

def train_random_forest(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    n_estimators: int = 100,
    max_depth: int = 10,
    min_samples_split: int = 5,
    min_samples_leaf: int = 2,
    random_state: int = 42
) -> RandomForestClassifier:
    """
    Treina modelo Random Forest.
    
    Returns:
        RandomForestClassifier: Modelo treinado
    """
    model = RandomForestClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_samples_split=min_samples_split,
        min_samples_leaf=min_samples_leaf,
        random_state=random_state,
        n_jobs=-1
    )
    
    model.fit(X_train, y_train)
    return model


def train_xgboost(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    n_estimators: int = 100,
    max_depth: int = 5,
    learning_rate: float = 0.1,
    random_state: int = 42
) -> XGBClassifier:
    """
    Treina modelo XGBoost.
    
    Returns:
        XGBClassifier: Modelo treinado
    """
    if not XGBOOST_AVAILABLE:
        raise ImportError("XGBoost não está instalado. Instale com: pip install xgboost")
    
    model = XGBClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=learning_rate,
        random_state=random_state,
        eval_metric='logloss'
    )
    
    model.fit(X_train, y_train)
    return model


def train_logistic_regression(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    max_iter: int = 1000,
    C: float = 1.0,
    random_state: int = 42
) -> LogisticRegression:
    """
    Treina modelo Logistic Regression.
    
    Returns:
        LogisticRegression: Modelo treinado
    """
    model = LogisticRegression(
        max_iter=max_iter,
        C=C,
        random_state=random_state
    )
    
    model.fit(X_train, y_train)
    return model


def evaluate_model(
    model,
    X_test: pd.DataFrame,
    y_test: pd.Series,
    model_name: str = "Model"
) -> Dict:
    """
    Avalia modelo e retorna métricas.
    
    Returns:
        dict: Dicionário com métricas
    """
    y_pred = model.predict(X_test)
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    
    metrics = {
        'accuracy': accuracy_score(y_test, y_pred),
        'precision': precision_score(y_test, y_pred, zero_division=0),
        'recall': recall_score(y_test, y_pred, zero_division=0),
        'f1_score': f1_score(y_test, y_pred, zero_division=0),
        'roc_auc': roc_auc_score(y_test, y_pred_proba)
    }
    
    print(f"\n=== {model_name} ===")
    print(classification_report(y_test, y_pred))
    print(f"ROC-AUC: {metrics['roc_auc']:.4f}")
    
    return metrics


def predict_patient(
    model,
    patient_features: pd.DataFrame,
    scaler: Optional[StandardScaler] = None,
    apply_log: bool = True
) -> Dict:
    """
    Realiza predição para amostra do paciente com alinhamento de features.
    """

    patient_X = patient_features.copy()

    # 1️⃣ Garantir numérico
    patient_X = patient_X.apply(pd.to_numeric, errors="coerce")

    # 2️⃣ Log transform
    if apply_log:
        patient_X = np.log1p(patient_X)

    # 3️⃣ Alinhamento de features
    if scaler and hasattr(scaler, "feature_names_"):
        model_features = scaler.feature_names_

        aligned_X = pd.DataFrame(
            0,
            index=patient_X.index,
            columns=model_features
        )

        common_features = patient_X.columns.intersection(model_features)
        aligned_X[common_features] = patient_X[common_features]

        patient_X = aligned_X

    # 4️⃣ Scaling
    if scaler:
        patient_X = pd.DataFrame(
            scaler.transform(patient_X),
            columns=patient_X.columns,
            index=patient_X.index
        )

    # 5️⃣ Predição
    prediction = model.predict(patient_X)[0]
    probabilities = model.predict_proba(patient_X)[0]

    return {
        'prediction': int(prediction),
        'probability_covid': float(probabilities[1]),
        'probability_no_covid': float(probabilities[0])
    }


def get_feature_importance(model, feature_names: list) -> pd.DataFrame:
    """
    Extrai importância de features do modelo.
    
    Returns:
        pd.DataFrame: DataFrame com features e importâncias
    """
    if hasattr(model, 'feature_importances_'):
        importances = model.feature_importances_
    elif hasattr(model, 'coef_'):
        # Para Logistic Regression, usar valor absoluto dos coeficientes
        importances = np.abs(model.coef_[0])
    else:
        raise ValueError("Modelo não tem atributo de importância de features")
    
    df = pd.DataFrame({
        'feature': feature_names,
        'importance': importances
    }).sort_values('importance', ascending=False)
    
    return df
