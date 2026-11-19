// =====================================================================
// Ekstraksi curah hujan harian CHIRPS per kabupaten/kota Sumatera Utara
// Jalankan di Google Earth Engine Code Editor: https://code.earthengine.google.com
// Output: CSV format panjang (date, kabkota, precip_mm) ke Google Drive
// =====================================================================

// ---------- 1. Batas wilayah kabupaten/kota ----------
// GAUL 2015 hanya punya 25 kab/kota Sumut (pemekaran 2007-2008 belum masuk),
// jadi dipakai geoBoundaries (lebih baru). GAUL hanya untuk batas provinsi.
var prov = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM1_NAME', 'Sumatera Utara')).geometry();

var kab = ee.FeatureCollection('WM/geoLab/geoBoundaries/600/ADM2')
  .filter(ee.Filter.eq('shapeGroup', 'IDN'))
  .filterBounds(prov)
  .map(function (f) {               // ambil hanya yang titik tengahnya di Sumut
    return f.set('in_prov', prov.contains(f.geometry().centroid(100), 100));
  })
  .filter(ee.Filter.eq('in_prov', true))
  .map(function (f) { return f.set('ADM2_NAME', f.get('shapeName')); });

// LANGKAH PERTAMA: cek jumlah & nama kab/kota (harus 33)
print('Jumlah kab/kota:', kab.size());
print('Nama kab/kota:', kab.aggregate_array('ADM2_NAME').sort());
Map.centerObject(kab, 7);
Map.addLayer(kab, {color: 'blue'}, 'Kab/Kota Sumut');

// ---------- 2. Data CHIRPS harian ----------
var mulai   = '1995-01-01';   // 30 tahun -> periode normal klimatologis
var selesai = '2025-01-01';   // tanggal akhir tidak ikut dihitung

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterDate(mulai, selesai)
  .select('precipitation');

// ---------- 3. Rata-rata areal per kab/kota per hari ----------
var hasil = chirps.map(function (img) {
  var tgl = img.date().format('YYYY-MM-dd');
  return img.reduceRegions({
    collection: kab.select(['ADM2_NAME']),
    reducer: ee.Reducer.mean(),
    scale: 5566            // resolusi asli CHIRPS (~0,05 derajat)
  }).map(function (f) {
    return ee.Feature(null, {
      date: tgl,
      kabkota: f.get('ADM2_NAME'),
      precip_mm: f.get('mean')
    });
  });
}).flatten();

// ---------- 4. Ekspor ke Google Drive ----------
Export.table.toDrive({
  collection: hasil,
  description: 'CHIRPS_daily_kabkota33_Sumut_1995_2024',
  folder: 'GEE_export',
  fileFormat: 'CSV',
  selectors: ['date', 'kabkota', 'precip_mm']
});
// Setelah Run: buka tab "Tasks" (kanan atas) -> klik RUN pada tugas ekspor.
// Proses bisa memakan 15-60 menit. File muncul di Google Drive/GEE_export.
