"""
Zonasi Curah Hujan Sumatera Utara: MSSA + Ward HCA (versi Python)
=================================================================
Data  : CHIRPS harian 1995-2024, 33 kabupaten/kota (rata-rata areal, GEE)
Alur  : 1) standardisasi  2) MSSA (setara Rssa kind="mssa")  3) W-correlation
        4) rekonstruksi sinyal musiman  5) Ward HCA + Silhouette & CH
        6) pembanding (mentah harian, klimatologi bulanan) + ARI
        7) profil musiman per zona  8) validasi BPS (opsional)

Implementasi MSSA mengikuti Golyandina et al. (Rssa): trajectory matrix
disusun horizontal L x (C*K), dekomposisi eigen matriks L x L. Hasilnya
telah dicocokkan dengan Rssa (angka identik).

Pakai:  python zonasi_mssa_sumut.py
Butuh:  numpy pandas scipy scikit-learn matplotlib openpyxl
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.lib.stride_tricks import sliding_window_view
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.stats import kruskal
from sklearn.metrics import silhouette_score, calinski_harabasz_score, adjusted_rand_score

# ---------------------------------------------------------------- 0. KONFIGURASI
PATH_DATA = "Precip_CHIRPS_Sumut_33kab_1995_2024.xlsx"
PATH_BPS = "bps_produktivitas_padi_sumut.csv"   # kolom: Kabupaten_Kota, Produktivitas
L = 365                  # window length = satu siklus tahunan
N_KOMP = 12              # komponen yang diperiksa
GRUP_MUSIMAN = [0, 1, 2, 3]   # RC1-4 (indeks mulai 0): tahunan + semi-tahunan
K_RANGE = range(2, 9)
K_FINAL = 3
OUT = Path("hasil_zonasi_py"); OUT.mkdir(exist_ok=True)

ZONE_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
MARKERS = ["o", "s", "^", "D", "v", "P"]
INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e4e3de"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})

# ---------------------------------------------------------------- 1. DATA
df = pd.read_excel(PATH_DATA)
df["date"] = pd.to_datetime(df["date"])
df = df.set_index("date").sort_index()
names = df.columns.tolist()
P = df.to_numpy(float)
N, C = P.shape
print(f"Periode {df.index.min():%Y-%m-%d} s.d. {df.index.max():%Y-%m-%d} | N={N} hari | C={C} | NA={np.isnan(P).sum()}")
print(f"Proporsi hari tanpa hujan: {(P == 0).mean():.3f}")
Z = (P - P.mean(0)) / P.std(0, ddof=1)

# ---------------------------------------------------------------- 2. MSSA
K = N - L + 1
traj = [sliding_window_view(Z[:, c], L).T for c in range(C)]      # tiap L x K
S = sum(X @ X.T for X in traj)                                    # L x L
ev, U = np.linalg.eigh(S)
order = ev.argsort()[::-1]; ev, U = ev[order], U[:, order]
share = ev / ev.sum()

def hankel_avg(M):
    out = np.zeros(N); cnt = np.zeros(N)
    for i in range(L):
        out[i:i + K] += M[i]; cnt[i:i + K] += 1
    return out / cnt

def elementary_rc(j):
    u = U[:, [j]]
    return np.column_stack([hankel_avg(u @ (u.T @ X)) for X in traj])   # N x C

print("Rekonstruksi komponen elementer...")
RC = [elementary_rc(j) for j in range(N_KOMP)]

def dominant_period(x):
    x = x - x.mean(); spec = np.abs(np.fft.rfft(x))[1:]
    return N / (spec.argmax() + 1)

tab_komp = pd.DataFrame({
    "Komponen": np.arange(1, N_KOMP + 1),
    "Kontribusi_%": (100 * share[:N_KOMP]).round(2),
    "Kumulatif_%": (100 * np.cumsum(share[:N_KOMP])).round(2),
    "Periode_hari": [round(dominant_period(r.mean(1)), 1) for r in RC],
})
print(tab_komp.to_string(index=False))

# ---------------------------------------------------------------- 3. W-CORRELATION
w_t = np.array([min(t + 1, L, K, N - t) for t in range(N)], float)
w_all = np.tile(w_t, C)
flat = [r.ravel("F") for r in RC]
def wcor(a, b):
    return (a * b * w_all).sum() / np.sqrt((a * a * w_all).sum() * (b * b * w_all).sum())
W = np.array([[wcor(flat[i], flat[j]) for j in range(N_KOMP)] for i in range(N_KOMP)])
tab_wcor = pd.DataFrame(np.abs(W).round(3), index=[f"RC{i+1}" for i in range(N_KOMP)],
                        columns=[f"RC{i+1}" for i in range(N_KOMP)])

# ---------------------------------------------------------------- 4-5. CLUSTERING
seasonal = sum(RC[j] for j in GRUP_MUSIMAN)                       # N x C
print(f"Varians dijelaskan grup musiman: {100 * share[GRUP_MUSIMAN].sum():.2f}%")

def row_z(F):
    return (F - F.mean(1, keepdims=True)) / F.std(1, ddof=1, keepdims=True)

def evaluate(F, label):
    F = row_z(F); Lk = linkage(F, "ward"); rows, labels = [], {}
    for k in K_RANGE:
        lab = fcluster(Lk, k, "maxclust"); labels[k] = lab
        sizes = "-".join(map(str, sorted(np.bincount(lab)[1:], reverse=True)))
        rows.append(dict(Data=label, k=k, Silhouette=round(silhouette_score(F, lab), 3),
                         CH=round(calinski_harabasz_score(F, lab), 1), Ukuran=sizes))
    return Lk, labels, pd.DataFrame(rows)

monthly = df.resample("MS").sum()
clim = monthly.groupby(monthly.index.month).mean()               # 12 x C (mm/bulan)

Lk_m, lab_m, t_m = evaluate(seasonal.T, "MSSA (RC1-4)")
_, lab_r, t_r = evaluate(Z.T, "Mentah harian")
_, lab_c, t_c = evaluate(clim.T.to_numpy(), "Klimatologi bulanan")
tab_valid = pd.concat([t_m, t_r, t_c], ignore_index=True)
tab_ari = pd.DataFrame([dict(k=k,
    ARI_MSSA_vs_Mentah=round(adjusted_rand_score(lab_m[k], lab_r[k]), 2),
    ARI_MSSA_vs_Klim=round(adjusted_rand_score(lab_m[k], lab_c[k]), 2)) for k in K_RANGE])
print(tab_valid.to_string(index=False)); print(tab_ari.to_string(index=False))

# ---------------------------------------------------------------- 6. ZONA FINAL
# Beri nomor zona berurutan dari wilayah barat-tengah ke timur agar stabil:
raw_lab = lab_m[K_FINAL]
mean_annual = clim.sum()
zone_order = sorted(set(raw_lab), key=lambda z: -np.mean(mean_annual.values[raw_lab == z]))
remap = {z: i + 1 for i, z in enumerate(zone_order)}           # zona 1 = paling basah
zona = pd.DataFrame({"Kabupaten_Kota": names, "Zona": [remap[z] for z in raw_lab]})
zona["Curah_hujan_tahunan_mm"] = mean_annual.values.round(0)
zona = zona.sort_values(["Zona", "Kabupaten_Kota"])
print(zona.to_string(index=False))

prof = (clim.T.assign(Zona=[remap[z] for z in raw_lab]).groupby("Zona").agg(["mean", "std"]))
tab_prof = pd.DataFrame({
    "Zona": np.repeat(prof.index, 12), "Bulan": np.tile(np.arange(1, 13), len(prof)),
    "Rata2_mm": np.concatenate([prof.xs("mean", axis=1, level=1).loc[z].values for z in prof.index]).round(1),
    "SD_mm": np.concatenate([prof.xs("std", axis=1, level=1).loc[z].values for z in prof.index]).round(1),
})

# ---------------------------------------------------------------- 7. GAMBAR
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]

# Fig 1: spektrum kontribusi + label periode
fig, ax = plt.subplots(figsize=(7.5, 3.6))
cols = ["#2a78d6" if j in GRUP_MUSIMAN else "#b9b8b0" for j in range(N_KOMP)]
ax.bar(tab_komp.Komponen, tab_komp["Kontribusi_%"], color=cols, width=0.7, edgecolor="white", linewidth=2)
for _, r in tab_komp.iterrows():
    p = r.Periode_hari
    lbl = f"{p/365.25:.1f} y" if p > 400 else f"{p:.0f} d"
    ax.text(r.Komponen, r["Kontribusi_%"] + 0.05, lbl, ha="center", va="bottom", fontsize=7.5, color=MUTED)
ax.set_xticks(tab_komp.Komponen); ax.set_xlabel("Component"); ax.set_ylabel("Variance explained (%)")
ax.grid(axis="x", visible=False); ax.set_axisbelow(True)
ax.set_title("MSSA eigenvalue spectrum (dark bars = seasonal group RC1–4; labels = dominant period)",
             fontsize=10, loc="left", color=INK)
fig.savefig(OUT / "Fig1_eigen_spectrum.png"); plt.close(fig)

# Fig 2: W-correlation heatmap (satu hue, terang -> gelap)
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list("blue", ["#fcfcfb", "#cde2fb", "#6da7ec", "#2a78d6", "#104281"])
fig, ax = plt.subplots(figsize=(5.2, 4.6))
im = ax.imshow(np.abs(W), cmap=cmap, vmin=0, vmax=1)
ax.set_xticks(range(N_KOMP)); ax.set_yticks(range(N_KOMP))
ax.set_xticklabels([f"{i+1}" for i in range(N_KOMP)]); ax.set_yticklabels([f"{i+1}" for i in range(N_KOMP)])
ax.set_xlabel("Reconstructed component"); ax.set_ylabel("Reconstructed component"); ax.grid(False)
cb = fig.colorbar(im, ax=ax, shrink=0.8); cb.set_label("|W-correlation|"); cb.outline.set_visible(False)
ax.set_title("W-correlation matrix", fontsize=10, loc="left")
fig.savefig(OUT / "Fig2_wcorrelation.png"); plt.close(fig)

# Fig 3: validasi klaster (dua panel, satu sumbu masing-masing)
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
series = ["MSSA (RC1-4)", "Mentah harian", "Klimatologi bulanan"]
legend_en = {"MSSA (RC1-4)": "MSSA seasonal (RC1–4)", "Mentah harian": "Raw daily", "Klimatologi bulanan": "Monthly climatology"}
for metric, ax in zip(["Silhouette", "CH"], axes):
    for i, s in enumerate(series):
        d = tab_valid[tab_valid.Data == s]
        ax.plot(d.k, d[metric], color=ZONE_COLORS[i], lw=2, marker=MARKERS[i], ms=6,
                markeredgecolor="white", markeredgewidth=1.5, label=legend_en[s])
    ax.axvline(K_FINAL, color=MUTED, lw=1, ls=":")
    ax.set_xlabel("Number of clusters (k)"); ax.set_ylabel(metric if metric == "Silhouette" else "Calinski–Harabasz")
    ax.set_xticks(list(K_RANGE))
axes[0].legend(frameon=False, fontsize=8.5, loc="center right")
axes[0].set_title("Average silhouette", fontsize=10, loc="left"); axes[1].set_title("Calinski–Harabasz index", fontsize=10, loc="left")
fig.tight_layout(); fig.savefig(OUT / "Fig3_cluster_validation.png"); plt.close(fig)

# Fig 4: dendrogram
fig, ax = plt.subplots(figsize=(9, 4.2))
cut = Lk_m[-(K_FINAL - 1), 2] - 1e-9
zone_of_leaf = {i: remap[z] for i, z in enumerate(raw_lab)}
node_zone = dict(zone_of_leaf)
for i, (a, b, _, _) in enumerate(Lk_m):
    za, zb = node_zone[int(a)], node_zone[int(b)]
    node_zone[C + i] = za if za == zb else None
def link_color(node):
    z = node_zone.get(node)
    return ZONE_COLORS[z - 1] if z else "#9a998f"
dendrogram(Lk_m, labels=names, leaf_rotation=90, leaf_font_size=7.5, link_color_func=link_color, ax=ax)
for tick in ax.get_xticklabels():
    tick.set_color(INK)
ax.axhline(cut, color=MUTED, lw=1, ls=":"); ax.grid(False); ax.set_ylabel("Ward linkage distance")
ax.set_title(f"Ward dendrogram of MSSA seasonal signals (cut at k = {K_FINAL})", fontsize=10, loc="left")
fig.savefig(OUT / "Fig4_dendrogram.png"); plt.close(fig)

# Fig 5: profil musiman per zona (label langsung di ujung garis)
fig, ax = plt.subplots(figsize=(7.5, 4))
for z in sorted(prof.index):
    m = prof.xs("mean", axis=1, level=1).loc[z].values
    n_mem = (zona.Zona == z).sum()
    ax.plot(range(1, 13), m, color=ZONE_COLORS[z - 1], lw=2, marker=MARKERS[z - 1], ms=6,
            markeredgecolor="white", markeredgewidth=1.5, label=f"Region {z} (n = {n_mem})")
ax.set_xticks(range(1, 13)); ax.set_xticklabels(months); ax.set_xlim(0.6, 12.4); ax.set_axisbelow(True)
ax.set_ylabel("Mean monthly precipitation (mm)"); ax.legend(frameon=False, fontsize=8.5, loc="upper left")
ax.set_title("Seasonal rainfall regime by region (CHIRPS 1995–2024)", fontsize=10, loc="left")
fig.savefig(OUT / "Fig5_seasonal_profile.png"); plt.close(fig)

# ---------------------------------------------------------------- 8. VALIDASI BPS (opsional)
if Path(PATH_BPS).exists():
    bps = pd.read_csv(PATH_BPS); bps["Kabupaten_Kota"] = bps["Kabupaten_Kota"].str.strip()
    rows = []
    for nama, lab in [("MSSA", lab_m), ("Mentah harian", lab_r), ("Klimatologi bulanan", lab_c)]:
        d = pd.DataFrame({"Kabupaten_Kota": names, "Zona": lab[K_FINAL]}).merge(bps, on="Kabupaten_Kota")
        groups = [g.Produktivitas.values for _, g in d.groupby("Zona")]
        H, p = kruskal(*groups); kk = len(groups)
        rows.append(dict(Partisi=nama, n=len(d), H=round(H, 2), p=round(p, 4),
                         epsilon2=round((H - kk + 1) / (len(d) - kk), 3)))
    tab_bps = pd.DataFrame(rows); print(tab_bps.to_string(index=False))
    tab_bps.to_csv(OUT / "Tabel_validasi_BPS.csv", index=False)
else:
    print("File BPS belum ada; validasi BPS dilewati.")

# ---------------------------------------------------------------- 9. EKSPOR TABEL
with pd.ExcelWriter(OUT / "Tabel_hasil_zonasi.xlsx") as xw:
    tab_komp.to_excel(xw, sheet_name="Komponen", index=False)
    tab_wcor.to_excel(xw, sheet_name="W_correlation")
    tab_valid.to_excel(xw, sheet_name="Validasi_klaster", index=False)
    tab_ari.to_excel(xw, sheet_name="ARI", index=False)
    zona.to_excel(xw, sheet_name="Zona_kabkota", index=False)
    tab_prof.to_excel(xw, sheet_name="Profil_musiman", index=False)
print(f"Selesai. Output di folder: {OUT.resolve()}")
