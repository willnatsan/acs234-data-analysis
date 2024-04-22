data = readmatrix("LondonAirQuality(2010_2018)_v2.xlsx");
X = data(:, 9:13);
y = data(:, 15);

mdl = fitlm(X, y)
