# =====================================================================
# Zonasi Curah Hujan Sumatera Utara: MSSA + Ward HCA
# Alur lengkap & reproducible (versi revisi untuk ICDSOS 2026)
#
# Struktur:
#   0. Konfigurasi
#   1. Import & standardisasi data
#   2. MSSA: kontribusi varians, periode komponen, W-correlation
#   3. Rekonstruksi sinyal musiman
#   4. Clustering Ward + validasi (Silhouette & Calinski-Harabasz)
#   5. Pembanding: data mentah harian & klimatologi bulanan (+ ARI)
#   6. Profil musiman tiap zona
#   7. Ekspor hasil
#   8. Validasi dengan statistik resmi (BPS) -- isi setelah data tersedia
# =====================================================================

library(Rssa)
library(readxl)
library(cluster)     # silhouette()
library(tidyverse)

# ---------------------------------------------------------------------
# 0. KONFIGURASI  (hanya bagian ini yang perlu diubah saat ganti data)
# ---------------------------------------------------------------------
path_data  <- "Precip_CHIRPS_Sumut_33kab_1995_2024.xlsx"   # CHIRPS harian 1995-2024, 33 kab/kota
L          <- 365        # window length: siklus tahunan
n_komp     <- 20         # jumlah komponen yang diperiksa
k_range    <- 2:8        # rentang jumlah klaster yang dievaluasi
out_dir    <- "hasil_zonasi"
dir.create(out_dir, showWarnings = FALSE)
set.seed(2026)

# ---------------------------------------------------------------------
# 1. IMPORT & STANDARDISASI
# ---------------------------------------------------------------------
raw <- read_excel(path_data)
names(raw) <- trimws(names(raw))                 # buang spasi di nama kolom
raw$date <- as.Date(raw$date)
raw <- raw %>% mutate(across(-date, as.numeric))  # kolom Excel terbaca sebagai teks

P <- raw %>% select(-date) %>% as.matrix()
cat("Periode:", format(min(raw$date)), "s.d.", format(max(raw$date)),
    "| N =", nrow(P), "hari | C =", ncol(P), "kab/kota\n")
cat("Missing value:", sum(is.na(P)), "\n")
cat("Proporsi hari tanpa hujan:", round(mean(P == 0), 3), "\n")

Z <- scale(P)                                    # z-score per kab/kota

# ---------------------------------------------------------------------
# 2. MSSA
# ---------------------------------------------------------------------
s <- ssa(Z, L = L, kind = "mssa", neig = n_komp)

# 2a. Kontribusi varians tiap komponen (laporkan apa adanya di paper)
kontrib <- contributions(s)[1:n_komp]
tabel_kontrib <- tibble(
  Komponen     = 1:n_komp,
  Kontribusi   = round(100 * kontrib, 2),
  Kumulatif    = round(100 * cumsum(kontrib), 2)
)
print(tabel_kontrib)

png(file.path(out_dir, "Fig_eigen_spectrum.png"), 900, 600)
plot(s, type = "values", numvalues = n_komp, main = "MSSA eigenvalue spectrum")
dev.off()

# 2b. Periode dominan tiap komponen (untuk memberi label: tahunan / semi-tahunan / tren)
rc_elem <- reconstruct(s, groups = as.list(1:n_komp))
periode_dominan <- function(x) {
  x  <- x - mean(x)
  sp <- Mod(fft(x))[2:floor(length(x) / 2)]
  length(x) / which.max(sp)                      # dalam hari
}
tabel_periode <- tibble(
  Komponen = 1:n_komp,
  Periode_hari = sapply(1:n_komp, function(j)
    round(periode_dominan(rowMeans(as.matrix(rc_elem[[j]]))), 1))
)
print(tabel_periode)
# Catatan (CHIRPS 1995-2024): RC1-2 = tahunan (365 hari), RC3-4 = semi-tahunan
# (183 hari), RC5 = komponen antartahunan (~7,5 tahun). Cek ulang jika data diganti.

# 2c. W-correlation
w <- wcor(s, groups = 1:n_komp)
png(file.path(out_dir, "Fig_wcorrelation.png"), 800, 800)
plot(w, main = "W-correlation matrix (RC1-20)")
dev.off()
print(round(abs(as.matrix(w))[1:8, 1:8], 2))

# ---------------------------------------------------------------------
# 3. REKONSTRUKSI SINYAL MUSIMAN
#    Tentukan grup BERDASARKAN tabel_periode & W-correlation di atas.
# ---------------------------------------------------------------------
grup_musiman <- 1:4              # tahunan (RC1-2) + semi-tahunan (RC3-4)
grup_tren    <- 5                # komponen antartahunan, tidak ikut clustering

rec <- reconstruct(s, groups = list(musiman = grup_musiman, tren = grup_tren))
S_musiman <- as.matrix(rec$musiman)              # N x C
cat("Varians yang dijelaskan grup musiman:",
    round(100 * sum(kontrib[grup_musiman]), 2), "%\n")

# ---------------------------------------------------------------------
# 4. CLUSTERING WARD + VALIDASI
# ---------------------------------------------------------------------
# Row-wise z-score: klaster berdasarkan pola/fase, bukan volume
row_z <- function(M) t(scale(t(M)))

# Calinski-Harabasz untuk label hasil cutree
ch_index <- function(X, lab) {
  k <- length(unique(lab)); n <- nrow(X); m <- colMeans(X)
  B <- 0; W <- 0
  for (g in unique(lab)) {
    Xg <- X[lab == g, , drop = FALSE]; mg <- colMeans(Xg)
    B <- B + nrow(Xg) * sum((mg - m)^2)
    W <- W + sum(sweep(Xg, 2, mg)^2)
  }
  (B / (k - 1)) / (W / (n - k))
}

evaluasi_klaster <- function(F, label_data) {
  F  <- row_z(F)
  d  <- dist(F)
  hc <- hclust(d, method = "ward.D2")
  hasil <- map_dfr(k_range, function(k) {
    lab <- cutree(hc, k)
    tibble(Data = label_data, k = k,
           Silhouette = round(mean(silhouette(lab, d)[, 3]), 3),
           CH = round(ch_index(F, lab), 1),
           Ukuran = paste(sort(table(lab), decreasing = TRUE), collapse = "-"))
  })
  list(hc = hc, tabel = hasil)
}

F_mssa <- t(S_musiman)                           # baris = kab/kota
ev_mssa <- evaluasi_klaster(F_mssa, "MSSA (RC musiman)")
print(ev_mssa$tabel)

# ---------------------------------------------------------------------
# 5. PEMBANDING
# ---------------------------------------------------------------------
# 5a. Data mentah harian (standar)
ev_raw <- evaluasi_klaster(t(Z), "Mentah harian")

# 5b. Klimatologi bulanan (rata-rata curah hujan per bulan kalender)
klim <- raw %>%
  mutate(th = year(date), bln = month(date)) %>%
  pivot_longer(-c(date, th, bln), names_to = "kab", values_to = "pr") %>%
  group_by(kab, th, bln) %>% summarise(pr = sum(pr), .groups = "drop") %>%
  group_by(kab, bln) %>% summarise(pr = mean(pr), .groups = "drop") %>%
  pivot_wider(names_from = bln, values_from = pr)
F_klim <- as.matrix(klim[, -1]); rownames(F_klim) <- klim$kab
F_klim <- F_klim[colnames(P), ]                  # samakan urutan kab/kota
ev_klim <- evaluasi_klaster(F_klim, "Klimatologi bulanan")

tabel_banding <- bind_rows(ev_mssa$tabel, ev_raw$tabel, ev_klim$tabel)
print(tabel_banding, n = Inf)

# 5c. Kesesuaian partisi (Adjusted Rand Index)
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); c2 <- function(x) sum(choose(x, 2))
  idx <- c2(tab); ea <- c2(rowSums(tab)) * c2(colSums(tab)) / choose(n, 2)
  mx <- (c2(rowSums(tab)) + c2(colSums(tab))) / 2
  (idx - ea) / (mx - ea)
}
tabel_ari <- map_dfr(k_range, function(k) tibble(
  k = k,
  ARI_MSSA_vs_Mentah = round(ari(cutree(ev_mssa$hc, k), cutree(ev_raw$hc, k)), 2),
  ARI_MSSA_vs_Klim   = round(ari(cutree(ev_mssa$hc, k), cutree(ev_klim$hc, k)), 2)
))
print(tabel_ari)

# ---------------------------------------------------------------------
# 6. PILIH k & PROFIL MUSIMAN TIAP ZONA
# ---------------------------------------------------------------------
k_final <- 3      # silhouette maksimum pada data CHIRPS; tetapkan SETELAH melihat tabel_banding (aturan ditulis di Metode)

zona <- tibble(Kabupaten_Kota = colnames(P), Zona = cutree(ev_mssa$hc, k_final))
print(zona %>% arrange(Zona), n = Inf)
print(table(zona$Zona))

png(file.path(out_dir, "Fig_dendrogram.png"), 1200, 700)
plot(ev_mssa$hc, main = "Ward HCA on MSSA seasonal reconstruction", xlab = "", sub = "")
rect.hclust(ev_mssa$hc, k = k_final, border = "blue")
dev.off()

profil <- as_tibble(F_klim, rownames = "Kabupaten_Kota") %>%
  pivot_longer(-Kabupaten_Kota, names_to = "Bulan", values_to = "mm") %>%
  mutate(Bulan = as.integer(Bulan)) %>%
  left_join(zona, by = "Kabupaten_Kota") %>%
  group_by(Zona, Bulan) %>%
  summarise(mean = mean(mm), sd = sd(mm), .groups = "drop")

g <- ggplot(profil, aes(Bulan, mean, colour = factor(Zona))) +
  geom_line(linewidth = 1) + geom_point() +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(x = NULL, y = "Mean monthly precipitation (mm)", colour = "Zone") +
  theme_minimal(base_size = 13)
ggsave(file.path(out_dir, "Fig_profil_musiman_zona.png"), g, width = 9, height = 5, dpi = 300)

# ---------------------------------------------------------------------
# 7. EKSPOR
# ---------------------------------------------------------------------
write_csv(tabel_kontrib, file.path(out_dir, "Tabel_kontribusi_varians.csv"))
write_csv(tabel_periode, file.path(out_dir, "Tabel_periode_komponen.csv"))
write_csv(tabel_banding, file.path(out_dir, "Tabel_validasi_klaster.csv"))
write_csv(tabel_ari,     file.path(out_dir, "Tabel_ARI.csv"))
write_csv(zona,          file.path(out_dir, "Tabel_zona_kabkota.csv"))
write_csv(profil,        file.path(out_dir, "Tabel_profil_musiman_zona.csv"))

# ---------------------------------------------------------------------
# 8. VALIDASI DENGAN STATISTIK RESMI (BPS)
#    Siapkan CSV: Kabupaten_Kota, Produktivitas (ku/ha) -- nama kab/kota
#    harus sama persis dengan kolom data curah hujan.
# ---------------------------------------------------------------------
path_bps <- "bps_produktivitas_padi_sumut.csv"
if (file.exists(path_bps)) {
  bps <- read_csv(path_bps) %>% mutate(Kabupaten_Kota = trimws(Kabupaten_Kota))

  uji_partisi <- function(lab, nama) {
    d <- tibble(Kabupaten_Kota = colnames(P), Zona = factor(lab)) %>%
      inner_join(bps, by = "Kabupaten_Kota")
    kw <- kruskal.test(Produktivitas ~ Zona, data = d)
    eps2 <- (kw$statistic - nlevels(d$Zona) + 1) / (nrow(d) - nlevels(d$Zona))
    tibble(Partisi = nama, n = nrow(d), H = round(kw$statistic, 2),
           p = signif(kw$p.value, 3), epsilon2 = round(eps2, 3))
  }
  tabel_bps <- bind_rows(
    uji_partisi(cutree(ev_mssa$hc, k_final), "MSSA"),
    uji_partisi(cutree(ev_raw$hc,  k_final), "Mentah harian"),
    uji_partisi(cutree(ev_klim$hc, k_final), "Klimatologi bulanan")
  )
  print(tabel_bps)
  write_csv(tabel_bps, file.path(out_dir, "Tabel_validasi_BPS.csv"))
} else {
  cat("File BPS belum ada; bagian 8 dilewati.\n")
}

cat("Selesai. Semua output ada di folder:", out_dir, "\n")
