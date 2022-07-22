# CSV2netCDF
Using this mini script you can readily convert specific CSV data to netCDF.

For running the script use the following pattern:

```
Rscript --vanilla test_.R <coordination csv data> <main csv file> <directory for output fiiles> <directory of gdal> <directory of cdo>
```

For example for a random case the follwoing code will be used:

```
Rscript --vanilla test_.R /Users/MosHad91/Downloads/ukcp09_grid_coords.txt /Users/MosHad91/Downloads/PrLnd_Abs_Pmean_Med_SCP_Scen1.txt /Users/MosHad91/Anna_s_proj/ /Applications/QGIS.app/Contents/MacOS/bin/ /opt/homebrew/bin/
```
