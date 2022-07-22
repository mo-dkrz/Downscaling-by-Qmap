#Rscript --vanilla test_.R /Users/hadizadeh-m/Downloads/ukcp09_grid_coords.txt /Users/hadizadeh-m/Downloads/PrLnd_Abs_Pmean_Med_SCP_Scen1.txt /Users/hadizadeh-m/Anna_s_proj/ /Applications/QGIS.app/Contents/MacOS/bin/ /opt/homebrew/bin/

packages = c("dplyr", "stringr",
             "rgdal", "sp")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
args = commandArgs(trailingOnly=TRUE)



ukcp09_grid_coords <- read.csv(args[1], header=FALSE, sep=",")
PrLnd_Abs_Pmean_Med_SCP_Scen1 <- read.csv(args[2], header=TRUE, sep="\t")
df <- data.frame(ukcp09_grid_coords['V2'],ukcp09_grid_coords['V3'],PrLnd_Abs_Pmean_Med_SCP_Scen1[,2:ncol(PrLnd_Abs_Pmean_Med_SCP_Scen1)],check.names = TRUE)
#print(df)

for (num_col in 3:ncol(PrLnd_Abs_Pmean_Med_SCP_Scen1)){
    data_mod1 <- cbind(df[0:2], stack(df[num_col:num_col]))
    sx = str_split(str_remove(colnames(df[num_col]),"X"), "_", simplify=TRUE)
    #print(ca)
    ca = as.POSIXct(paste(sx[,1],sx[,2],01,sep="-"), format = "%Y-%B-%d")
    #print(ca)
    names(data_mod1)[2] <- "lat"
    names(data_mod1)[1] <- "long"
    data_mod1$ind <- NULL
    write.csv(data_mod1,gsub(" ", "", paste(args[3],ca,".csv")), row.names = FALSE)
    MyData <- read.csv(file=gsub(" ", "", paste(args[3],ca,".csv")), header=TRUE, sep=",")
    class(MyData)
    coordinates(MyData)<-~long+lat
    class(MyData)
    writeOGR(MyData, gsub(" ", "", paste(args[3],ca,".shp")), "values", driver = "ESRI Shapefile")
    print(paste0(args[4],'gdal_translate -of NetCDF', gsub(" ", "", paste0(args[3],ca,".tiff")), gsub(" ", "", paste0(args[3],ca,".nc"))))
    system(paste0(args[4],'gdal_rasterize -l ',ca,' -a values -tr 0.5 0.5 -a_nodata 0.0 -te -15.072848689 48.880267594 9.579930555 61.557098453 -ot Float32 -of GTiff ', gsub(" ", "", paste(args[3],ca,".shp"))," ",gsub(" ", "", paste(args[3],ca,".tiff"))))
    system(paste0(args[4],'gdal_translate -of NetCDF ', gsub(" ", "", paste0(args[3],ca,".tiff"))," ",gsub(" ", "", paste0(args[3],ca,".nc"))))
    system(paste0(args[5],'cdo -setreftime,1900-01-01,00:00:00,1day -setdate,',ca, ' -setcalendar,standard'," ", gsub(" ", "", paste(args[3],ca,".nc")) ," ", gsub(" ", "", paste0(args[3],"final_",ca,".nc"))))


}
system(paste0(args[5],"cdo -f nc -b F32 mergetime ",args[3],"final_*.nc ",args[3],"final.nc"))
