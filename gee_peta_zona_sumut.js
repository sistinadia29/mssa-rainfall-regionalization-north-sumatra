// =====================================================================
// Peta zona curah hujan Sumatera Utara (hasil MSSA + Ward, k = 3)
// Jalankan di https://code.earthengine.google.com
// Hasil: (1) pratinjau peta di layar, (2) GeoJSON batas + zona ke Google Drive
// =====================================================================

// ---------- 1. Batas wilayah (sama dengan saat ekstraksi CHIRPS) ----------
var prov = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM1_NAME', 'Sumatera Utara')).geometry();

var gb = ee.FeatureCollection('WM/geoLab/geoBoundaries/600/ADM2')
  .filter(ee.Filter.eq('shapeGroup', 'IDN'))
  .filterBounds(prov)
  .map(function (f) {
    return f.set('in_prov', prov.contains(f.geometry().centroid(100), 100));
  })
  .filter(ee.Filter.eq('in_prov', true))
  .filter(ee.Filter.neq('shapeName', 'Danau Toba'))          // poligon danau
  .map(function (f) { return ee.Feature(f.geometry(), {nama: f.get('shapeName')}); });

var niasSel = ee.FeatureCollection('FAO/GAUL/2015/level2')
  .filter(ee.Filter.eq('ADM1_NAME', 'Sumatera Utara'))
  .filter(ee.Filter.eq('ADM2_NAME', 'Nias Selatan'))
  .map(function (f) { return ee.Feature(f.geometry(), {nama: 'Nias Selatan'}); });

var kab = gb.merge(niasSel);

// ---------- 2. Zona hasil analisis (nama sesuai ejaan batas wilayah) ----------
var zona = ee.Dictionary({
  // Zona 1 - Kepulauan Nias
  'Nias': 1, 'Nias Selatan': 1, 'Nias Barat': 1, 'Nias Utara': 1, 'Kota Gunungsitoli': 1,
  // Zona 2 - Lereng barat & dataran tinggi
  'Dairi': 2, 'Pakpak Bharat': 2, 'Samosir': 2, 'Toba Samosir': 2, 'Humbang Hasundutan': 2,
  'Tapanuli Utara': 2, 'Tapanuli Tengah': 2, 'Tapanuli Selatan': 2, 'Mandailing Natal': 2,
  'Padang Lawas': 2, 'Padang Lawas Utara': 2, 'Kota Sibolga': 2, 'Kota Padangsidimpuan': 2,
  // Zona 3 - Dataran rendah timur
  'Langkat': 3, 'Kota Binjai': 3, 'Kota Medan': 3, 'Deli Serdang': 3, 'Karo': 3,
  'Serdang Bedagai': 3, 'Kota Tebing Tinggi': 3, 'Batu Bara': 3, 'Simalungun': 3,
  'Kota Pematang Siantar': 3, 'Asahan': 3, 'Kota Tanjung Balai': 3, 'Labuhan Batu': 3,
  'Labuhan Batu Utara': 3, 'Labuhan Batu Selatan': 3
});

var warna = ee.Dictionary({'1': '2a78d6', '2': 'eb6834', '3': '1baf7a', '0': '999999'});

var kabZona = kab.map(function (f) {
  var z = ee.Number(zona.get(f.get('nama'), 0));        // 0 = nama tidak cocok
  return f.set({
    zona: z,
    style: {color: 'ffffff', width: 1,
            fillColor: ee.String(warna.get(z.format('%d'))).cat('dd')}
  });
});

// CEK: harus 33 wilayah dan 0 yang tidak cocok
print('Jumlah wilayah:', kabZona.size());
print('Tidak cocok (zona 0):', kabZona.filter(ee.Filter.eq('zona', 0)).aggregate_array('nama'));
print('Jumlah per zona:', kabZona.aggregate_histogram('zona'));

// ---------- 3. Pratinjau peta ----------
Map.setOptions('SATELLITE');
Map.centerObject(kabZona, 7);
Map.addLayer(kabZona.style({styleProperty: 'style'}), {}, 'Zona curah hujan');

var legenda = ui.Panel({style: {position: 'bottom-left', padding: '8px 12px'}});
legenda.add(ui.Label('Precipitation zones', {fontWeight: 'bold'}));
[['2a78d6', 'Zone 1 - Nias archipelago'],
 ['eb6834', 'Zone 2 - Western slopes & highlands'],
 ['1baf7a', 'Zone 3 - Eastern lowlands']].forEach(function (it) {
  legenda.add(ui.Panel([
    ui.Label('', {backgroundColor: '#' + it[0], padding: '8px', margin: '2px 6px 2px 0'}),
    ui.Label(it[1], {margin: '2px 0'})
  ], ui.Panel.Layout.Flow('horizontal')));
});
Map.add(legenda);

// ---------- 4. Ekspor batas + zona (untuk peta publikasi) ----------
var ekspor = kabZona.map(function (f) {
  return ee.Feature(f.geometry().simplify(200), {nama: f.get('nama'), zona: f.get('zona')});
});
Export.table.toDrive({
  collection: ekspor,
  description: 'zona_curah_hujan_sumut',
  folder: 'GEE_export',
  fileFormat: 'GeoJSON'
});
// Tab "Tasks" -> RUN. File zona_curah_hujan_sumut.geojson muncul di Drive/GEE_export.
