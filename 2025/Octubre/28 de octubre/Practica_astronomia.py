
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
FILEPATH = "asu (2).tsv" #Modificar segun el nombre de donde tega los datos guardados.
OUTPUT_PREFIX = "cluster25_hyades"
BV_MIN, BV_MAX = 0.4, 1.0

SIGMA_CLIP = 2.0
MAX_ITERS = 5
def load_vizier_tsv(path):
    df = pd.read_csv(path, sep=';', comment='#', header=None, engine='python')
    df.columns = ['B-V', 'Vmag'] if df.shape[1] >= 2 else ['B-V']
    for col in ['B-V', 'Vmag']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df[['B-V','Vmag']].dropna().copy()
    df = df[(df['B-V'] > -0.3) & (df['B-V'] < 2.2)]
    df = df[(df['Vmag'] > -2.0) & (df['Vmag'] < 22.0)]
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
    sel = (df['B-V'] >= bv_min) & (df['B-V'] <= bv_max)
    sub = df.loc[sel].copy()
    if sub.empty:
        raise ValueError("No hay datos en el rango de color especificado. Ajusta BV_MIN/BV_MAX.")

    ref_bv, ref_mv = main_sequence_reference_Mv()
    Mv_ref = interpolate_Mv_from_BV(sub['B-V'].values, ref_bv, ref_mv)
    good = np.isfinite(Mv_ref)
    sub = sub.loc[good].copy()
    Mv_ref = Mv_ref[good]
    mM = sub['Vmag'].values - Mv_ref

    mask = sigma_clip(mM, sigma=sigma, max_iters=max_iters)
    mM_used = mM[mask]
    if mM_used.size == 0:
        raise ValueError("Sigma-clipping eliminó todos los puntos. Relaja parámetros o revisa datos.")

    mM_med = np.median(mM_used)
    mM_std = np.std(mM_used, ddof=1) if mM_used.size > 1 else np.nan
    return mM_med, mM_std, mM_used.size, sub.iloc[mask].index.values

def modulus_to_distance_pc(m_minus_M):
    return 10 ** ((m_minus_M + 5.0) / 5.0)

def distance_uncertainty_pc(mM_med, mM_std):
    if not np.isfinite(mM_std):
        return np.nan
    D = modulus_to_distance_pc(mM_med)
    factor = np.log(10.0) / 5.0
    return factor * D * mM_std

def plot_color_magnitude(df, used_idx=None, savepath=None, title="Color–Magnitud (V vs B–V)"):
    plt.figure(figsize=(6, 6))
    plt.scatter(df['B-V'], df['Vmag'], s=14, alpha=0.6, c="#6666cc", label="Datos (crudos)")
    if used_idx is not None and len(used_idx) > 0:
        sel = np.zeros(len(df), dtype=bool)
        sel[used_idx] = True
        plt.scatter(df.loc[sel, 'B-V'], df.loc[sel, 'Vmag'], s=28, edgecolor='k', facecolor="#ffcc00",
                    label="Usados (ajuste)")

    plt.gca().invert_yaxis()
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
    if not os.path.exists(FILEPATH):
        raise FileNotFoundError(f"No encuentro el archivo: {FILEPATH}")
    df = load_vizier_tsv(FILEPATH)
    mM_med, mM_std, n_eff, used_idx = estimate_distance_modulus(
        df,
        bv_min=BV_MIN, bv_max=BV_MAX,
        sigma=SIGMA_CLIP, max_iters=MAX_ITERS
    )
    D_pc = modulus_to_distance_pc(mM_med)
    D_sigma_pc = distance_uncertainty_pc(mM_med, mM_std)
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
    cleaned_csv = f"{OUTPUT_PREFIX}_clean.csv"
    df.to_csv(cleaned_csv, index=False)
    print(f"\nCSV limpio guardado en: {cleaned_csv}")

    plot_path = f"{OUTPUT_PREFIX}_CM.png"
    plot_color_magnitude(df, used_idx=used_idx, savepath=plot_path,
                         title="Pléyades (Cl=20) — Diagrama Color–Magnitud")
    print(f"Gráfica guardada en: {plot_path}")

if __name__ == "__main__":
    main()
