# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# https://www.tec.ac.cr

# Eng. MSc. Maikel Mendez Morales
# email: maikel.mendez@gmail.com mamendez@itcr.ac.cr
# Github: https://github.com/maikelonu
# Google Scholar: https://scholar.google.com/citations?hl=en&user=JnmSVFYAAAAJ
# YouTube: https://www.youtube.com/c/MaikelMendez
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#----------------------------------------------------------------------------------------------------
# INFO: This script is intended for the generation of the final IDW climatology for the
# entire Costa Rican territory at 1x1 km and 25x25 km, both as continuos records and as individual months
#----------------------------------------------------------------------------------------------------
# INPUT FILES:

# TIFF format: 
# multiple *.tiff IDW3 monthly precipitation maps from 1960 to 1990

# ASC Raster format:
# mod_region.asc: Costa Rica territory mask 
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# OUTPUT FILES:

#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# MANUSCRIPT FIGURES:
# NONE
#----------------------------------------------------------------------------------------------------

# Workspace is cleared
#rm(list = ls())

# CRAN libraries are loaded
require(raster)
require(automap)
require(car)
require(chron)
require(compare)
require(DescTools)
require(dplyr)
require(fitdistrplus)
require(foreach)
require(forecast)
require(ggmap)
require(ggplot2)
require(gstat)
require(hyfo)
require(intamap)
require(KScorrect)
require(lattice)
require(lubridate)
require(maps)
require(maptools)
require(MASS)
require(ncdf4)
require(ncdf4.helpers)
require(pastecs)
require(PCICt)
require(plyr)
require(psych)
require(RColorBrewer)
require(rcompanion)
require(readr)
require(readxl)
require(reshape)
require(reshape2)
require(rgdal)
require(rgeos)
require(RNetCDF)
require(scales)
require(sp)
require(spatialEco)
require(tibble)
require(tidyr)
require(tseries)
require(viridis)
require(weathermetrics)

# Working directory is defined
setwd("/mnt/BACKUP/R_ITC/XLSX")  # UBUNTU-LINUX
setwd("C:/DATOS/R_ITC/XLSX")     # MS-WINDOWS

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: PRECIS georeference is loaded
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

# The blank map is transformed into SpatialGridDataFrame
raster.precis <- read.asciigrid("mod_region.asc", as.image = FALSE)

# The map is transformed into RasterLayer
raster.precis <- raster(raster.precis)

# CRTM05 projection for CR is created
CRTM05 <- CRS("+proj=tmerc +lat_0=0 +lon_0=-84 +k=0.9999 +x_0=500000 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# CRTM05 projection is assigned to rasterbrick object
crs(raster.precis) <- CRTM05

# Object attributes are requested
# It can be seen that spa.res is different !!!!
image(raster.precis, col = topo.colors(64))

# A simple plot is requested
plot(raster.precis)

# A summery is requested
summary(raster.precis)

# A dummy raster is created
raster.precis.dummy <- raster.precis

# Dummy NA values are introduced
values(raster.precis.dummy) <- NA

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: IMN Cubist Regression Temperature Climatology at 1x1 km
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Cubist Regression Temperature Climatology at 1x1 km is loaded 
IMN_temp_1x1_km <- stack("TEMP_1960_1990_MONTHLY.tif")

# A simple plot is created on raster level [1]
plot(IMN_temp_1x1_km[[1]])

# Dummy NAs values are assigned
NAvalue(IMN_temp_1x1_km) <- (-9999)

# CRTM05 projection is assigned to raster layer
crs(IMN_temp_1x1_km) <- CRTM05

#-----------------------------------------------------------------------------------
# A 25X25 km temperature climatology is compiled for the period 1961-1990
#-----------------------------------------------------------------------------------

# Bilinear resampling is executed
C_brk25km <- raster::resample(IMN_temp_1x1_km,
                              raster.precis.dummy,
                              resample = "bilinear")

# A simple plot is created on raster level [1]
plot(C_brk25km[[1]])

# # RasterLayers are masked based on PRECIS Climatic-Region domain
C_brk25km <- mask(C_brk25km, raster.precis)

# A simple plot is requested
plot(C_brk25km[[1]])

# Rasterbrick object is exported to GeoTIFF format
writeRaster(C_brk25km, "IMN_temp_25x25_km.tif", format="GTiff", overwrite=TRUE)

# Rasterbrick object is replicated N times to complete a Quasi-Extended Time-Series
C_brk25km_30 <- stack(replicate(30, C_brk25km))

# Rasterbrick object is exported to GeoTIFF format
writeRaster(C_brk25km_30, "IMN_temp_25x25_km_30N.tif", format="GTiff", overwrite=TRUE)

# GeoTIFF object is re-imported to avoid conflicts with netCDF
Hist_IMN_Temp_C <- stack("IMN_temp_25x25_km_30N.tif")

# Historical date vectors are created
dates_hist <- seq(as.Date("1961-01-01"), as.Date("1990-12-31"), by = "1 month")
dates_1951_1980 <- seq(as.Date("1951-01-01"), as.Date("1980-12-31"), by = "1 month")
dates_1961_1990 <- seq(as.Date("1961-01-01"), as.Date("1990-12-31"), by = "1 month")
dates_1981_1995 <- seq(as.Date("1981-01-01"), as.Date("1995-12-31"), by = "1 month")
dates_1991_2005 <- seq(as.Date("1991-01-01"), as.Date("2005-12-31"), by = "1 month")
dates_2041_2070 <- seq(as.Date("2041-01-01"), as.Date("2070-12-31"), by = "1 month")
dates_2071_2100 <- seq(as.Date("2071-01-01"), as.Date("2100-12-31"), by = "1 month")
dates_2011_2040 <- seq(as.Date("2011-01-01"), as.Date("2040-12-31"), by = "1 month")

# Dates are printed into Rasterbrick object (1961_1990)
Hist_IMN_Temp_C <- setZ(Hist_IMN_Temp_C, dates_1961_1990, name="time")

# Rasterbrick object is exported to netCDF format
writeRaster(Hist_IMN_Temp_C,"IMN_temp_25x25_km_30N.nc", format = "CDF", overwrite=TRUE, 
            varname = "Precipitation", varunit = "mm", longname = "totals", xname = "lon",
            yname = "lat", zname = "time", zunit = "numeric")

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: IMN Cubist Regression ET0 Climatology at 1x1 km
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Cubist Regression Temperature Climatology at 1x1 km is loaded 
IMN_ET0_1x1_km <- stack("ET0_1960_1990_MONTHLY.tif")

# A simple plot is created on raster level [1]
plot(IMN_ET0_1x1_km[[1]])

# Dummy NAs values are assigned
NAvalue(IMN_ET0_1x1_km) <- (-9999)

# CRTM05 projection is assigned to raster layer
crs(IMN_ET0_1x1_km) <- CRTM05

#-----------------------------------------------------------------------------------
# A 25X25 km ET0erature climatology is compiled for the period 1961-1990
#-----------------------------------------------------------------------------------

# Bilinear resampling is executed
C_brk25km <- raster::resample(IMN_ET0_1x1_km,
                              raster.precis.dummy,
                              resample = "bilinear")

# A simple plot is created on raster level [1]
plot(C_brk25km[[1]])

# # RasterLayers are masked based on PRECIS Climatic-Region domain
C_brk25km <- mask(C_brk25km, raster.precis)

# A simple plot is requested
plot(C_brk25km[[12]])

# Rasterbrick object is exported to GeoTIFF format
writeRaster(C_brk25km, "IMN_ET0_25x25_km.tif", format="GTiff", overwrite=TRUE)

# Rasterbrick object is replicated N times to complete a Quasi-Extended Time-Series
C_brk25km_30 <- stack(replicate(30, C_brk25km))

# Rasterbrick object is exported to GeoTIFF format
writeRaster(C_brk25km_30, "IMN_ET0_25x25_km_30N.tif", format="GTiff", overwrite=TRUE)

# GeoTIFF object is re-imported to avoid conflicts with netCDF
Hist_IMN_Temp_C <- stack("IMN_ET0_25x25_km_30N.tif")

# Historical date vectors are created
dates_hist <- seq(as.Date("1961-01-01"), as.Date("1990-12-31"), by = "1 month")
dates_1951_1980 <- seq(as.Date("1951-01-01"), as.Date("1980-12-31"), by = "1 month")
dates_1961_1990 <- seq(as.Date("1961-01-01"), as.Date("1990-12-31"), by = "1 month")
dates_1981_1995 <- seq(as.Date("1981-01-01"), as.Date("1995-12-31"), by = "1 month")
dates_1991_2005 <- seq(as.Date("1991-01-01"), as.Date("2005-12-31"), by = "1 month")
dates_2041_2070 <- seq(as.Date("2041-01-01"), as.Date("2070-12-31"), by = "1 month")
dates_2071_2100 <- seq(as.Date("2071-01-01"), as.Date("2100-12-31"), by = "1 month")
dates_2011_2040 <- seq(as.Date("2011-01-01"), as.Date("2040-12-31"), by = "1 month")

# Dates are printed into Rasterbrick object (1961_1990)
Hist_IMN_Temp_C <- setZ(Hist_IMN_Temp_C, dates_1961_1990, name="time")

# Rasterbrick object is exported to netCDF format
writeRaster(Hist_IMN_Temp_C,"IMN_ET0_25x25_km_30N.nc", format = "CDF", overwrite=TRUE, 
            varname = "Precipitation", varunit = "mm", longname = "totals", xname = "lon",
            yname = "lat", zname = "time", zunit = "numeric")

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: PRECIS Historial Temperature at 25x25 km
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

# A temporal modelled raster layer is created
tempX <- stack("temp_C_1960_1990_hadgem.nc")

# Rasterbrick object is reprojected from WGS84 to CRTM05
tempX <- projectRaster(tempX, raster.precis.dummy, method="bilinear")

# Rasterbrick object is masked
tempX <- mask(tempX, raster.precis)

# Year 1960 is deleted
tempX <- tempX[[13:372]]

# A simple plot is requested
plot(tempX[[1]])

# RasterLayer is exported to GeoTIFF format
writeRaster(tempX, "PRECIS_temp_25x25_km_30N.tif", format="GTiff", overwrite = TRUE)

# Dates are printed into Rasterbrick object
tempX <- setZ(tempX, dates_1961_1990, name="time")

# Rasterbrick object is exported to netCDF format
writeRaster(tempX,"PRECIS_temp_25x25_km_30N.nc", format = "CDF", overwrite=TRUE, 
            varname = "Precipitation", varunit = "mm", longname = "totals", xname = "lon",
            yname = "lat", zname = "time", zunit = "numeric")

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: PRECIS RPCs Bias Correction
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

# A temporal modelled raster layer is created from original netCDF from OpenSUSE
tempX <- stack("temp_C_RPC26_sssaa_2011_2040.nc")

# Rasterbrick object is reprojected from WGS84 to CRTM05
tempX <- projectRaster(tempX, raster.precis.dummy, method="bilinear")

# Rasterbrick object is masked
tempX <- mask(tempX, raster.precis)

# A simple plot is requested
plot(tempX[[1]])

# Variable names are removed
names(tempX) <- NULL

# Dates are printed into Rasterbrick object
# dates_2071_2100 <- seq(as.Date("2071-01-01"), as.Date("2099-11-01"), by = "1 month") # for 2071-2100 ONLY !!!
tempX <- setZ(tempX, dates_2011_2040, name="time")

# Rasterbrick object is exported to netCDF format
writeRaster(tempX,"Temporal_Temperature.nc", format = "CDF", overwrite=TRUE,
            varname = "Precipitation", varunit = "mm", longname = "totals",
            xname = "lon", yname = "lat", zname = "time", zunit = "numeric")

#------------------------------------------------------------------------------------------------------
# SUB-BLOCK: {hyfo} Bias Correction
#------------------------------------------------------------------------------------------------------

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 1)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 1)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 1)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_1.nc")

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 2)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 2)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 2)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_2.nc")

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 3)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 3)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 3)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_3.nc")

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 4)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 4)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 4)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_4.nc")

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 5)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 5)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 5)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_5.nc")

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 6)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 6)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 6)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_6.nc")

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 7)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 7)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 7)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_7.nc")

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 8)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 8)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 8)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_8.nc")

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 9)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 9)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 9)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_9.nc")

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 10)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 10)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 10)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_10.nc")

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 11)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 11)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 11)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_11.nc")

# loadNcdf {hyfo} is used to load netCDF files
Temporal_Temperature <- loadNcdf(filePath="Temporal_Temperature.nc", varname="Precipitation", month = 12)
PRECIS_Temperature <- loadNcdf(filePath="PRECIS_temp_25x25_km_30N.nc", varname="Precipitation", month = 12)
IMN_Temperature <- loadNcdf(filePath="IMN_temp_25x25_km_30N.nc", varname="Precipitation", month = 12)

# biasCorrect {hyfo} is used for DELTA-BC
Temporal_Temperature_BC <- biasCorrect(frc = Temporal_Temperature,    # RCP to be calibrated
                                       hindcast = PRECIS_Temperature, # PRECIS raw data
                                       obs = IMN_Temperature,         # IMN observations
                                       method = "scaling",
                                       extrapolate = "no",
                                       preci = FALSE,
                                       scaleType = 'add')

# writeNcdf {hyfo} is used to export the bias-corrected NetCDFs for QGIS
writeNcdf(gridData = Temporal_Temperature_BC, filePath = "test_12.nc")

#------------------------------------------------------------------------------------------------------
# END of SUB-BLOCK
#------------------------------------------------------------------------------------------------------

# A temporal modelled rasterstack object is created
tempMOD1 <- stack("test_1.nc")
tempMOD2 <- stack("test_2.nc")
tempMOD3 <- stack("test_3.nc")
tempMOD4 <- stack("test_4.nc")
tempMOD5 <- stack("test_5.nc")
tempMOD6 <- stack("test_6.nc")
tempMOD7 <- stack("test_7.nc")
tempMOD8 <- stack("test_8.nc")
tempMOD9 <- stack("test_9.nc")
tempMOD10 <- stack("test_10.nc")
tempMOD11 <- stack("test_11.nc")
tempMOD12 <- stack("test_12.nc")

# CRTM05 projection is assigned to rasterbrick object
crs(tempMOD1) <- CRTM05
crs(tempMOD2) <- CRTM05
crs(tempMOD3) <- CRTM05
crs(tempMOD4) <- CRTM05
crs(tempMOD5) <- CRTM05
crs(tempMOD6) <- CRTM05
crs(tempMOD7) <- CRTM05
crs(tempMOD8) <- CRTM05
crs(tempMOD9) <- CRTM05
crs(tempMOD10) <- CRTM05
crs(tempMOD11) <- CRTM05
crs(tempMOD12) <- CRTM05

# Spatial bias is calculated (modelled - observed!!)
l.observed <- C_brk25km # or C_brk25km_30 or C_brk25km
#l.modelled <- tempMOD

# Monthly variables are loaded for 2071-2100 ONLY !!!!!!!
vs.JAN <- c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,289,301,313,325,337)
vs.FEB <- c(2,14,26,38,50,62,74,86,98,110,122,134,146,158,170,182,194,206,218,230,242,254,266,278,290,302,314,326,338)
vs.MAR <- c(3,15,27,39,51,63,75,87,99,111,123,135,147,159,171,183,195,207,219,231,243,255,267,279,291,303,315,327,339)
vs.APR <- c(4,16,28,40,52,64,76,88,100,112,124,136,148,160,172,184,196,208,220,232,244,256,268,280,292,304,316,328,340)
vs.MAY <- c(5,17,29,41,53,65,77,89,101,113,125,137,149,161,173,185,197,209,221,233,245,257,269,281,293,305,317,329,341)
vs.JUN <- c(6,18,30,42,54,66,78,90,102,114,126,138,150,162,174,186,198,210,222,234,246,258,270,282,294,306,318,330,342)
vs.JUL <- c(7,19,31,43,55,67,79,91,103,115,127,139,151,163,175,187,199,211,223,235,247,259,271,283,295,307,319,331,343)
vs.AUG <- c(8,20,32,44,56,68,80,92,104,116,128,140,152,164,176,188,200,212,224,236,248,260,272,284,296,308,320,332,344)
vs.SEP <- c(9,21,33,45,57,69,81,93,105,117,129,141,153,165,177,189,201,213,225,237,249,261,273,285,297,309,321,333,345)
vs.OCT <- c(10,22,34,46,58,70,82,94,106,118,130,142,154,166,178,190,202,214,226,238,250,262,274,286,298,310,322,334,346)
vs.NOV <- c(11,23,35,47,59,71,83,95,107,119,131,143,155,167,179,191,203,215,227,239,251,263,275,287,299,311,323,335,347)
vs.DEC <- c(12,24,36,48,60,72,84,96,108,120,132,144,156,168,180,192,204,216,228,240,252,264,276,288,300,312,324,336)

# Monthly variables are loaded
vs.JAN <- c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,289,301,313,325,337,349)
vs.FEB <- c(2,14,26,38,50,62,74,86,98,110,122,134,146,158,170,182,194,206,218,230,242,254,266,278,290,302,314,326,338,350)
vs.MAR <- c(3,15,27,39,51,63,75,87,99,111,123,135,147,159,171,183,195,207,219,231,243,255,267,279,291,303,315,327,339,351)
vs.APR <- c(4,16,28,40,52,64,76,88,100,112,124,136,148,160,172,184,196,208,220,232,244,256,268,280,292,304,316,328,340,352)
vs.MAY <- c(5,17,29,41,53,65,77,89,101,113,125,137,149,161,173,185,197,209,221,233,245,257,269,281,293,305,317,329,341,353)
vs.JUN <- c(6,18,30,42,54,66,78,90,102,114,126,138,150,162,174,186,198,210,222,234,246,258,270,282,294,306,318,330,342,354)
vs.JUL <- c(7,19,31,43,55,67,79,91,103,115,127,139,151,163,175,187,199,211,223,235,247,259,271,283,295,307,319,331,343,355)
vs.AUG <- c(8,20,32,44,56,68,80,92,104,116,128,140,152,164,176,188,200,212,224,236,248,260,272,284,296,308,320,332,344,356)
vs.SEP <- c(9,21,33,45,57,69,81,93,105,117,129,141,153,165,177,189,201,213,225,237,249,261,273,285,297,309,321,333,345,357)
vs.OCT <- c(10,22,34,46,58,70,82,94,106,118,130,142,154,166,178,190,202,214,226,238,250,262,274,286,298,310,322,334,346,358)
vs.NOV <- c(11,23,35,47,59,71,83,95,107,119,131,143,155,167,179,191,203,215,227,239,251,263,275,287,299,311,323,335,347,359)
vs.DEC <- c(12,24,36,48,60,72,84,96,108,120,132,144,156,168,180,192,204,216,228,240,252,264,276,288,300,312,324,336,348,360)

# Mean spatial bias is calculated (modelled - observed!!)
of.BIAS.JAN.C <- calc(tempMOD1, mean)
of.BIAS.JAN <- (of.BIAS.JAN.C - l.observed[[1]])

of.BIAS.FEB.C <- calc(tempMOD2, mean)
of.BIAS.FEB <- (of.BIAS.FEB.C - l.observed[[2]])

of.BIAS.MAR.C <- calc(tempMOD3, mean)
of.BIAS.MAR <- (of.BIAS.MAR.C - l.observed[[3]])

of.BIAS.APR.C <- calc(tempMOD4, mean)
of.BIAS.APR <- (of.BIAS.APR.C - l.observed[[4]])

of.BIAS.MAY.C <- calc(tempMOD5, mean)
of.BIAS.MAY <- (of.BIAS.MAY.C- l.observed[[5]])

of.BIAS.JUN.C <- calc(tempMOD6, mean)
of.BIAS.JUN <- (of.BIAS.JUN.C - l.observed[[6]])

of.BIAS.JUL.C <- calc(tempMOD7, mean)
of.BIAS.JUL <- (of.BIAS.JUL.C - l.observed[[7]])

of.BIAS.AUG.C <- calc(tempMOD8, mean)
of.BIAS.AUG <- (of.BIAS.AUG.C - l.observed[[8]])

of.BIAS.SEP.C <- calc(tempMOD9, mean)
of.BIAS.SEP <- (of.BIAS.SEP.C - l.observed[[9]])

of.BIAS.OCT.C <- calc(tempMOD10, mean)
of.BIAS.OCT <- (of.BIAS.OCT.C - l.observed[[10]])

of.BIAS.NOV.C <- calc(tempMOD11, mean)
of.BIAS.NOV <- (of.BIAS.NOV.C - l.observed[[11]])

of.BIAS.DEC.C <- calc(tempMOD12, mean)
of.BIAS.DEC <- (of.BIAS.DEC.C - l.observed[[12]])

# A Bias rasterstack object is created
bias_C_fut <- stack(of.BIAS.JAN,
                    of.BIAS.FEB,
                    of.BIAS.MAR,
                    of.BIAS.APR,
                    of.BIAS.MAY,
                    of.BIAS.JUN,
                    of.BIAS.JUL,
                    of.BIAS.AUG,
                    of.BIAS.SEP,
                    of.BIAS.OCT,
                    of.BIAS.NOV,
                    of.BIAS.DEC)

# A BC rasterstack object is created
EXPORT_BC <- stack(of.BIAS.JAN.C,
                   of.BIAS.FEB.C,
                   of.BIAS.MAR.C,
                   of.BIAS.APR.C,
                   of.BIAS.MAY.C,
                   of.BIAS.JUN.C,
                   of.BIAS.JUL.C,
                   of.BIAS.AUG.C,
                   of.BIAS.SEP.C,
                   of.BIAS.OCT.C,
                   of.BIAS.NOV.C,
                   of.BIAS.DEC.C)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Seasonal Climatic Regions Masking
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#----------------------------------------------------------------------------------------------------------------
# SUBBLOCK: ASC CR Climatic Regions Masking
#----------------------------------------------------------------------------------------------------------------
# North_Pacific = 3
# Central_Pacific = 2
# South_Pacific = 5
# Central = 4
# North = 1
# Caribbean = 6

# A sequence vector is created
seq.Vector <- c(3,2,5,4,1,6)

# The main loop is initiated
for (i in 1:6) {
  
  num.V <-seq.Vector[i]
  
  # Climatic-regions CR blank raster is loaded
  blank <- paste("mod_region.asc", sep = "/")
  
  # The blank map is transformed into SpatialGridDataFrame
  blank.grid <- read.asciigrid(blank, as.image = FALSE)
  
  # The map is transformed into RasterLayer
  blank.raster <- raster(blank.grid)
  
  # CRTM05 projection is assigned to rasterbrick object
  crs(blank.raster) <- CRTM05
  
  # Climatic region is selected
  values(blank.raster)[values(blank.raster) != num.V] <- NA
  
  # A simple plot is requested
  plot(blank.raster)
  
  #----------------------------------------------------------------------------------------------------------------
  # Season DJF
  #----------------------------------------------------------------------------------------------------------------
  rt.PC_DJF <- mask(bias_C_fut[[c(12,1,2)]], blank.raster)
  rt.PC_DJF <- cellStats(rt.PC_DJF, mean, na.rm = TRUE)
  rt.PC_DJF <- mean(rt.PC_DJF)
  
  #----------------------------------------------------------------------------------------------------------------
  # Season MAM
  #----------------------------------------------------------------------------------------------------------------
  rt.PC_MAM <- mask(bias_C_fut[[c(3,4,5)]], blank.raster)
  rt.PC_MAM <- cellStats(rt.PC_MAM, mean, na.rm = TRUE)
  rt.PC_MAM <- mean(rt.PC_MAM)
  
  #----------------------------------------------------------------------------------------------------------------
  # Season JJA
  #----------------------------------------------------------------------------------------------------------------
  rt.PC_JJA <- mask(bias_C_fut[[c(6,7,8)]], blank.raster)
  rt.PC_JJA <- cellStats(rt.PC_JJA, mean, na.rm = TRUE)
  rt.PC_JJA <- mean(rt.PC_JJA)
  
  #----------------------------------------------------------------------------------------------------------------
  # Season SON
  #----------------------------------------------------------------------------------------------------------------
  rt.PC_SON <- mask(bias_C_fut[[c(9,10,11)]], blank.raster)
  rt.PC_SON <- cellStats(rt.PC_SON, mean, na.rm = TRUE)
  rt.PC_SON <- mean(rt.PC_SON)
  
  # A compiled seasonal data.frame is created
  df.compiled <- c(rt.PC_DJF, rt.PC_MAM, rt.PC_JJA, rt.PC_SON)
  
  print(df.compiled)
  
} # Main loop is closed

# If (+) increase
# If (-) decrease

# BC Rasterbrick object is exported to GeoTIFF format
writeRaster(EXPORT_BC, "EXPORT_BC_RCP26_2011_2040_Mean.tif", format="GTiff", overwrite = TRUE)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: ggplot2
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#----------------------------------------------------------------------------------------------------------------
# Temperature
#----------------------------------------------------------------------------------------------------------------

# Input data is loaded and a data.frame is created
df.base <- read.table("temperature_future_rcps_SCALINGmean.txt", header = TRUE)

# A "levels" monthly vector is created for visualization puposes
v.select.name <- c("North_Pacific", "Central_Pacific", "South_Pacific","Central", "North", "Caribbean")

df.base$RCP <- as.factor(df.base$RCP)

df.base.melted <- melt(df.base)

# Months are ordered to avoid ggplot plotting conflicts
df.base.melted$Region <- ordered(df.base.melted$Region, v.select.name)

# Figure Temperature VIE
ggplot() +
  geom_point(aes(x = Region,y = value,shape = RCP,colour = Period),
             data=df.base.melted,size = 3.5, position = position_jitter(width = 0.15)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5.0,min.n = 5.0)) +
  theme_bw() +
  geom_hline(data=df.base,yintercept = 0.0,size = 0.75,linetype = 2,colour = 'gray') +
  geom_hline(data=df.base,yintercept = 1.5,size = 0.75,linetype = 2,colour = 'cornflowerblue') +
  geom_hline(data=df.base,yintercept = 2.0,size = 0.75,linetype = 2,colour = 'chocolate') +
  facet_grid(facets = variable ~ .) +
  theme(axis.text.x = element_text(angle = 0.0),
        text=element_text(size=16,  family="serif")) +
  xlab("Climatic Region") +
  ylab("Mean Daily Temperature (C)") 

# A color palette is created based on http://colorbrewer2.org/
scale03 <- c("#4daf4a", "#377eb8", "#e41a1c")

# Figure Temperature VIE
ggplot() +
  geom_hline(data=df.base.melted,yintercept = 1.5,size = 0.75,linetype = 2,colour = 'black') +
  geom_hline(data=df.base.melted,yintercept = 2.0,size = 0.75,linetype = 2,colour = 'black') +
  geom_point(aes(x = Period,y = value,shape = RCP, colour = RCP),
             data=df.base.melted, size = 4.0, position = position_jitter(width = 0.0001)) +
  facet_grid(facets = variable ~ Region + RCP, scales = 'free_y') +
  geom_line(aes(x = Period,y = value,colour = RCP,group = RCP),
            data=df.base.melted,size = 1.10) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_color_manual(values = scale03) +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 90.0), 
        text=element_text(size=13.5,  family="serif"), legend.position="bottom") +
  xlab("30-year period") +
  ylab("1961-1990_Baseline Change in Mean Daily Temperature [C]")

#----------------------------------------------------------------------------------------------------------------
# ET0
#----------------------------------------------------------------------------------------------------------------

# Input data is loaded and a data.frame is created
df.base2 <- read.table("ET0_future_rcps.txt", header = TRUE)

# A "levels" monthly vector is created for visualization puposes
v.select.name <- c("North_Pacific", "Central_Pacific", "South_Pacific","Central", "North", "Caribbean")

df.base2$RCP <- as.factor(df.base2$RCP)

df.base2.melted <- melt(df.base2)

# Months are ordered to avoid ggplot plotting conflicts
df.base2.melted$Region <- ordered(df.base2.melted$Region, v.select.name)

# Figure ET0 VIE
ggplot() +
  geom_point(aes(x = Region,y = value,shape = RCP,colour = Period),
             data=df.base2.melted,size = 3.5, position = position_jitter(width = 0.15)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5.0,min.n = 5.0)) +
  theme_bw() +
  geom_hline(data=df.base2,yintercept = 0.0,size = 0.75,linetype = 2,colour = 'gray') +
  geom_hline(data=df.base2,yintercept = 0.5,size = 0.75,linetype = 2,colour = 'cornflowerblue') +
  geom_hline(data=df.base2,yintercept = 1.0,size = 0.75,linetype = 2,colour = 'chocolate') +
  facet_grid(facets = variable ~ .) +
  theme(axis.text.x = element_text(angle = 0.0),
        text=element_text(size=16,  family="serif")) +
  xlab("Climatic Region") +
  ylab("Mean Monthly-Daily Et0 (mm/day)")

# Figure ET0 VIE
ggplot() +
  geom_hline(data=df.base2.melted,yintercept = 0.5,size = 0.75,linetype = 2,colour = 'black') +
  geom_hline(data=df.base2.melted,yintercept = 1.0,size = 0.75,linetype = 2,colour = 'black') +
  geom_point(aes(x = Period,y = value,shape = RCP, colour = RCP),
             data=df.base2.melted, size = 4.0, position = position_jitter(width = 0.0001)) +
  facet_grid(facets = variable ~ Region + RCP, scales = 'free_y') +
  geom_line(aes(x = Period,y = value,colour = RCP,group = RCP),
            data=df.base2.melted,size = 1.10) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_color_manual(values = scale03) +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 90.0), 
        text=element_text(size=13.5,  family="serif"), legend.position="bottom") +
  xlab("30-year period") +
  ylab("1961-1990_Baseline Change in Mean Monthly-Daily Et0 [mm/day]")

#--------------------------------------------------------------------------------------------------------------------------------------
# END of code
#--------------------------------------------------------------------------------------------------------------------------------------








