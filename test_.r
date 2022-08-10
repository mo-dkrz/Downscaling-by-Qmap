#Rscript --vanilla test_.R /Users/hadizadeh-m/Downloads/All_in_one/anna/ukcp09_grid_coords.txt /Users/hadizadeh-m/Downloads/All_in_one/anna/PrLnd_Abs_Pmean_Med_SCP_Scen1.txt /Users/hadizadeh-m/Downloads/All_in_one/anna/ /Applications/QGIS.app/Contents/MacOS/bin/ /opt/homebrew/bin/ /Users/hadizadeh-m/Downloads/All_in_one/anna/rainfall_hadukgrid_uk_1km_mon-30y_198101-201012.nc

packages = c("dplyr", "stringr",
             "rgdal", "sp", "ncdf4","qmap","reshape2","splitstackshape","janitor","tidyr","tibble")

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

time <- Sys.time()
print(paste("Start time:",time))
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
system(paste0(args[5],"cdo remapbil,",args[6]," final.nc ",args[3],"final_ress.nc"))

ncpath_gcm <- "final_ress.nc"
ncpath_obs <- args[6]

ncin_gcm <- nc_open(ncpath_gcm)
ncin_obs <- nc_open(ncpath_obs)
#print(ncin_gcm)
# get longitude and latitude


dname_obs <- "rainfall"  
dname_gcm <- "Band1"  

rain_arra_gcm <- ncvar_get(ncin_gcm,dname_gcm)
rain_arra_obs <- ncvar_get(ncin_obs,dname_obs)

t_gcm <- ncvar_get(ncin_gcm, "time")
t_obs <- ncvar_get(ncin_obs, "time")
gcmdatadates <- as.Date(ncin_gcm$dim$time$vals, origin = '1900-01-01')
#print(gcmdatadates)
obsdatadates <- as.Date(ncin_obs$dim$time$vals, origin = '1950-01-01')
# print(obsdatadates)

lon_gcm <- ncvar_get(ncin_gcm, "longitude")
lat_gcm <- ncvar_get(ncin_gcm, "latitude", verbose = F)

lon_obs <- ncvar_get(ncin_obs, "longitude")
lat_obs <- ncvar_get(ncin_obs, "latitude", verbose = F)

projection_x_obs <- ncvar_get(ncin_obs, "projection_x_coordinate")
projection_y_obs <- ncvar_get(ncin_obs, "projection_y_coordinate", verbose = F)

projection_x_gcm <- ncvar_get(ncin_gcm, "projection_x_coordinate")
projection_y_gcm <- ncvar_get(ncin_gcm, "projection_y_coordinate", verbose = F)


# print(projection_x_obs[397])
# print(projection_y_obs[433])
# print(t_gcm[1])
#obsoutput <- ncvar_get(ncin_gcm, varid = 'Band1',start= c(projection_y_obs[433],projection_x_obs[397],t_gcm[1]),

projxytimegcm <- as.matrix(expand.grid(projection_x_gcm,projection_y_gcm,t_gcm))
#lonlattime <- as.matrix(expand.grid(lon_obs,lat_obs,t_gcm))
projxytimeobs <- as.matrix(expand.grid(projection_x_obs,projection_y_obs,t_obs))

lswt_vec_long_gcm <- as.vector(rain_arra_gcm)
lswt_vec_long_obs <- as.vector(rain_arra_obs)


lswt_obs <- data.frame(cbind(projxytimeobs, lswt_vec_long_obs))
#lswt_obsx <- data.frame(cbind(lonlattime, lswt_vec_long))
lswt_obs <- na.omit(lswt_obs)

lswt_gcm <- data.frame(cbind(projxytimegcm, lswt_vec_long_gcm))
#lswt_obsx <- data.frame(cbind(lonlattime, lswt_vec_long))
lswt_gcm <- na.omit(lswt_gcm)


# print(nrow(lswt_obs))
# print(nrow(lswt_gcm))

time <- Sys.time()
print(paste("Start time for processing extracted data:",time))
lswt_obs$latlon <- paste0(lswt_obs$Var1,"-",lswt_obs$Var2)
lswt_gcm$latlon <- paste0(lswt_gcm$Var1,"-",lswt_gcm$Var2)

#head(lswt_obsx)
df_obs = subset(lswt_obs, select = -c(Var1,Var2))
df_gcm = subset(lswt_gcm, select = -c(Var1,Var2))

#head(df_obs)

df_obs <- df_obs[, c("Var3", "latlon", "lswt_vec_long_obs")]
df_gcm <- df_gcm[, c("Var3", "latlon", "lswt_vec_long_gcm")]

#head(df_obs)
time <- Sys.time()
print(paste("Start time for processing tables:",time))
df_obs.long <- dcast(df_obs, Var3~latlon, value.var = "lswt_vec_long_obs")
#df_gcm.long <- dcast(df_gcm, Var3~latlon, value.var = "lswt_vec_long")
df_gcm.long <- dcast(df_gcm, Var3~latlon, value.var = "lswt_vec_long_gcm")

#head(df_obs)

time <- Sys.time()
print(paste("Start time for processing bias correcting:",time))

final_data <- data.frame(date = gcmdatadates)

sqrtquant <- function(x,qstep=0.01){
  qq <- quantile(x, prob= seq(0,1,by=qstep))
  sqrt(qq)
}

for (id in 2:ncol(df_gcm.long)){
#for (id in 2:100){  
  qm.fit <- fitQmapQUANT(df_obs.long[,id],df_gcm.long[,id],qstep=0.1,nboot=1,wet.day=TRUE)
  qm.a <- doQmapQUANT(df_gcm.long[,id],qm.fit,type="linear")
  final_data[,id] <- qm.a
  colnames(final_data)[id] <- colnames(df_gcm.long)[id]
  print(id)
}
head(final_data)


final_data <- as.data.frame(t(final_data))

names(final_data) <- NULL
final_data <- final_data %>% row_to_names(row_number = 1)
final_data$date <- rownames(final_data)

head(final_data)



final_data <- final_data %>% separate(date, c('x', 'y'))
head(final_data)

time <- Sys.time()
print(paste("End time for processing bias correcting:",time))
write.csv(final_data, "final_bias_corrected_data.csv")



nc_close(ncin_gcm)
nc_close(ncin_obs)




