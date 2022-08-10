# Resampling and correcting bias of a netCDF 
Using this mini script you can readily convert specific CSV data to netCDF file with higher resolution and then the bias of model has been corrected by Qmap library. Prerequisites of this script are **CDO**, **gdal_translate**, and **gdal_rasterize**.

For running the script use the following pattern:

```
Rscript --vanilla test_.R <coordination csv data> <main csv file> <directory for output fiiles> <directory of gdal> <directory of cdo> <directory of observation file>
```

The follwoing is an example of using pattern:

```
Rscript --vanilla test_.R /Users/MosHad91/ukcp09_grid_coords.txt /Users/MosHad91/PrLnd_Abs_Pmean_Med_SCP_Scen1.txt /Users/MosHad91/ /Applications/QGIS.app/Contents/MacOS/bin/ /opt/homebrew/bin/ /Users/MosHad91/rainfall_hadukgrid_uk_1km_mon-30y_198101-201012.nc
```
