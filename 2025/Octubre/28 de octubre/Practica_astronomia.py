# -*- coding: utf-8 -*-
"""
Práctica 11 - Distancia a cúmulos por ajuste de Secuencia Principal (Main Sequence Fitting)
Entrada: TSV de VizieR (V/19) con comentarios (#) y separador ';' incluyendo filas de unidades.
Salidas:
 - Limpieza y filtrado de datos (B-V, Vmag)
 - Gráfica Color–Magnitud (invertida en V)
 - Estimación robusta de m-M (mediana tras sigma-clipping) y distancia D (pc)
 - Errores por dispersión (1σ) y N efectivo
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============= CONFIGURACIÓN =============
FILEPATH = "asu (2).tsv"
OUTPUT_PREFIX = "cluster25_hyades"
BV_MIN, BV_MAX = 0.4, 1.0  # rango F–K típico de Hyades


# Sigma-clipping para rechazar outliers en m-M
SIGMA_CLIP = 2.0
MAX_ITERS = 5

# ================= UTILIDAD =================
def load_vizier_tsv(path):
    """
    Lee TSV de VizieR con comentarios y filas extra (unidades, separadores).
    Devuelve DataFrame con columnas 'B-V' y 'Vmag' numéricas, sin NaN ni duplicados exactos.
    """
    # Leemos todo ignorando filas que empiezan con '#'
    df = pd.read_csv(path, sep=';', comment='#', header=None, engine='python')

    # Esperamos que las dos primeras filas remanentes sean etiquetas y unidades; forzamos nombres fijos.
    # Luego coaccionamos a numérico para barrer filas como 'mag' y '------'.
    df.columns = ['B-V', 'Vmag'] if df.shape[1] >= 2 else ['B-V']
    for col in ['B-V', 'Vmag']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Mantén solo filas válidas
    df = df[['B-V','Vmag']].dropna().copy()

    # Limpieza dura de rangos plausibles
    df = df[(df['B-V'] > -0.3) & (df['B-V'] < 2.2)]
    df = df[(df['Vmag'] > -2.0) & (df['Vmag'] < 22.0)]

    # Duplicados exactos fuera
    df = df.drop_duplicates().reset_index(drop=True)
    return df

def main_sequence_reference_Mv():
    """
    Locus MS alineado con las anclas del documento de clase:
      B-V=0.2 -> Mv=2.0
      B-V=0.5 -> Mv=3.5
    Curva suavizada alrededor de F–K:
    """
    bv = np.array([0.00, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 1.00, 1.20, 1.40], dtype=float)
    mv = np.array([0.9 , 1.4 , 2.0 , 2.9 , 3.5 , 4.5 , 5.6 , 6.7 , 7.9 , 9.0 ], dtype=float)
    return bv, mv


def interpolate_Mv_from_BV(bv_values, ref_bv, ref_mv):
    """
    Interpolación lineal de Mv en función de B-V usando los puntos de referencia de la SP.
    Extrapola linealmente en extremos si hace falta.
    """
    return np.interp(bv_values, ref_bv, ref_mv, left=np.nan, right=np.nan)

def sigma_clip(values, sigma=2.0, max_iters=5):
    """
    Sigma-clipping iterativo simple. Devuelve máscara booleana de valores retenidos.
    """
    mask = np.isfinite(values)
    vals = values[mask]
    for _ in range(max_iters):
        if vals.size < 3:
            break
        mu = np.nanmedian(vals)
        sd = np.nanstd(vals, ddof=1)
        if sd == 0 or not np.isfinite(sd):
            break
        keep = np.abs(vals - mu) <= sigma * sd
        if keep.all():
            break
        vals = vals[keep]
        # Actualiza mask global manteniendo los índices que siguen dentro
        temp = np.zeros_like(mask, dtype=bool)
        temp_indices = np.where(mask)[0][keep]
        temp[temp_indices] = True
        mask = temp
    return mask

def estimate_distance_modulus(df, bv_min=0.2, bv_max=1.2, sigma=2.0, max_iters=5):
    """
    Calcula m-M por ajuste vertical (V - Mv_ref(B-V)) en un rango de color.
    Aplica sigma-clipping y devuelve:
      mM_med, mM_std, N_eff, indices_usados (máscara)
    """
    # Filtra a un rango de color donde la SP es más limpia
    sel = (df['B-V'] >= bv_min) & (df['B-V'] <= bv_max)
    sub = df.loc[sel].copy()
    if sub.empty:
        raise ValueError("No hay datos en el rango de color especificado. Ajusta BV_MIN/BV_MAX.")

    ref_bv, ref_mv = main_sequence_reference_Mv()
    Mv_ref = interpolate_Mv_from_BV(sub['B-V'].values, ref_bv, ref_mv)

    # Algunas extrapolaciones pueden dar NaN si B-V cae fuera del rango de referencia
    good = np.isfinite(Mv_ref)
    sub = sub.loc[good].copy()
    Mv_ref = Mv_ref[good]

    # Distancia módulo por estrella: m - M
    mM = sub['Vmag'].values - Mv_ref

    # Sigma-clipping sobre mM
    mask = sigma_clip(mM, sigma=sigma, max_iters=max_iters)
    mM_used = mM[mask]
    if mM_used.size == 0:
        raise ValueError("Sigma-clipping eliminó todos los puntos. Relaja parámetros o revisa datos.")

    mM_med = np.median(mM_used)
    mM_std = np.std(mM_used, ddof=1) if mM_used.size > 1 else np.nan
    return mM_med, mM_std, mM_used.size, sub.iloc[mask].index.values

def modulus_to_distance_pc(m_minus_M):
    """Convierte módulo de distancia a parsecs."""
    return 10 ** ((m_minus_M + 5.0) / 5.0)

def distance_uncertainty_pc(mM_med, mM_std):
    """
    Propaga incertidumbre: dD/d(m-M) = ln(10)/5 * D
    Devuelve sigma_D (pc). Si mM_std es NaN, devuelve NaN.
    """
    if not np.isfinite(mM_std):
        return np.nan
    D = modulus_to_distance_pc(mM_med)
    factor = np.log(10.0) / 5.0
    return factor * D * mM_std

def plot_color_magnitude(df, used_idx=None, savepath=None, title="Color–Magnitud (V vs B–V)"):
    plt.figure(figsize=(6, 6))
    # Todos los puntos
    plt.scatter(df['B-V'], df['Vmag'], s=14, alpha=0.6, c="#6666cc", label="Datos (crudos)")
    # Puntos usados en el ajuste
    if used_idx is not None and len(used_idx) > 0:
        sel = np.zeros(len(df), dtype=bool)
        sel[used_idx] = True
        plt.scatter(df.loc[sel, 'B-V'], df.loc[sel, 'Vmag'], s=28, edgecolor='k', facecolor="#ffcc00",
                    label="Usados (ajuste)")

    plt.gca().invert_yaxis()  # magnitud: menor es más brillante → arriba
    plt.xlabel("B − V (mag)")
    plt.ylabel("V (mag)")
    plt.title(title)
    plt.legend()
    plt.grid(alpha=0.3)
    if savepath:
        plt.tight_layout()
        plt.savefig(savepath, dpi=200)
    plt.show()

def main():
    # 1) Cargar y limpiar
    if not os.path.exists(FILEPATH):
        raise FileNotFoundError(f"No encuentro el archivo: {FILEPATH}")
    df = load_vizier_tsv(FILEPATH)

    # 2) Estimar módulo de distancia y distancia
    mM_med, mM_std, n_eff, used_idx = estimate_distance_modulus(
        df,
        bv_min=BV_MIN, bv_max=BV_MAX,
        sigma=SIGMA_CLIP, max_iters=MAX_ITERS
    )
    D_pc = modulus_to_distance_pc(mM_med)
    D_sigma_pc = distance_uncertainty_pc(mM_med, mM_std)

    # 3) Reporte en consola
    print("\n=== RESULTADOS (Cluster 20: Pléyades) ===")
    print(f"m - M (mediana, sigma-clipped): {mM_med:.3f} mag")
    if np.isfinite(mM_std):
        print(f"Dispersión 1σ en m - M     : {mM_std:.3f} mag (N={n_eff})")
    else:
        print(f"Dispersión 1σ en m - M     : N/A (N={n_eff})")
    print(f"Distancia D                : {D_pc:.1f} pc")
    if np.isfinite(D_sigma_pc):
        print(f"Incertidumbre σ_D         : ±{D_sigma_pc:.1f} pc")
        print(f"Rango aproximado          : [{D_pc - D_sigma_pc:.1f}, {D_pc + D_sigma_pc:.1f}] pc")

    # 4) Guardar CSV limpio y gráfico
    cleaned_csv = f"{OUTPUT_PREFIX}_clean.csv"
    df.to_csv(cleaned_csv, index=False)
    print(f"\nCSV limpio guardado en: {cleaned_csv}")

    plot_path = f"{OUTPUT_PREFIX}_CM.png"
    plot_color_magnitude(df, used_idx=used_idx, savepath=plot_path,
                         title="Pléyades (Cl=20) — Diagrama Color–Magnitud")
    print(f"Gráfica guardada en: {plot_path}")

if __name__ == "__main__":
    main()
