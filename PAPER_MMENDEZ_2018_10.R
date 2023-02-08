# ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# Construction Engineering School
# https://www.tec.ac.cr
# Eng. MSc. Maikel Mendez Morales
# Email: maikel.mendez@gmail.com; mamendez@itcr.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://github.com/maikelonu
# Skype: maikel.mendez
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

#----------------------------------------------------------------------------------------------------
# INFO: This script is intended for temperature cross validation of the following interpolation methods:
# IDW:       Inverse Distance Weighting
# CUBIST:    Cubist regression surface
# ORK:       Ordinary Kriging
# KED:       Kriging with External Drift
# RF:        Random Forest
#----------------------------------------------------------------------------------------------------
# INPUT FILES:
# To be defined
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# OUTPUT FILES:
# To be defined
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# MANUSCRIPT FIGURES:
# To be defined
#----------------------------------------------------------------------------------------------------

# Workspace is cleared
#rm(list = ls())

# CRAN libraries are loaded
require(raster)
require(automap)
require(caret)
require(Cubist)
require(DescTools)
require(devtools)
require(dplyr)
require(FME)
require(geoR)
require(ggplot2)
require(ggthemes)
require(glmulti)
require(gridExtra)
require(GSIF)
require(gstat)
require(lattice)
require(latticeExtra)
require(mlbench)
require(openair)
require(plotrix)
require(plyr)
require(psych)
require(randomForest)
require(ranger)
require(rasterVis)
require(reshape)
require(reshape2)
require(rgdal)
require(scales)
require(sp)
require(VSURF)

# Working directory is defined
setwd("/mnt/BACKUP/R_ITC/XLSX")  # UBUNTU-LINUX
setwd("C:/DATOS/R_ITC/XLSX")     # MS-WINDOWS

# Seed is set to make partitions reproducible
def.seed <- 123
set.seed(def.seed)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Preliminary Variable Selection Using Random Forests
# /////////////////////////////////////////////////////////////////////////////////////

# Input data is loaded and a data.frame is created
# Covariants matrix created in SAGA-GIS and based on ALOS-DEM 1kmx1km
# Relevant predictive field is temperature "TEMP"
df.base <- read.table("vsurf_complete.txt", header = TRUE)

# Pairwise relationships are visualised with a pairs plot in order to analyse
# potencial correlation among variables
pairs.panels(df.base[, -c(1,2,11,15)], hist.col="#ca0020", lm=TRUE, jiggle=TRUE)

# Each layer coming form a different raster layer

#[1]  "X": CRTM05 X position [m]
#[2]  "Y": CRTM05 Y position [m]       
#[3]  "SLOPE": terrain slope [radians] (Tarboton 1997)
#[4]  "ALOS_FILL": ALOS AW3D-30m DEM [m] (Tadono et al. 2016)
#[5]  "ASPECT": [radians] (Tarboton 1997)
#[6]  "TOPOGRAPHIC": Topographic Wetness Index [-] (Beven and Kirkby 1979)
#[7]  "LS_FACTOR":  Slope length (LS) factor as used by the Universal Soil Loss Equation (USLE) [-] (Boehner and Selige 2016)
#[8]  "DIRECT_INSO": Potential incoming solar radiation (insolation) [kW/m2] (Hofierka and Suri 2002)
#[9]  "CONVEXITY": Terrain Surface Convexity for subsequent terrain classification [-] (Iwahashi and Pike 2007)
#[10] "WIND_EXPOSI": Wind Exposition Index calculates the average 'Wind Effect Index' for all directions using an angular step (Gerlitz et al. 2015)
#[11] "BUFFER": Proximity-Grid with euclidean distance to feature cells (observations) [m]
#[12] "TEMP": Mean monthly temperature from IMN data base for the period 1961-1990 [C] (Rojas 1985)
#[13] "COAST": Coastline Proximity-Grid with euclidean distance to feature cells (observations) [m]
#[14] "PREC": Mean monthly precipitation [mm/month] (Mendez et al. 2019)

# If NEEDED!!!, a monthly subset can be performed, to analyse monthly behavior
# instead of a lumped yearly behavior
df.base <- subset(df.base, MONTH == "JAN")

# Irrelevant variables are deleted
#df.base <- df.base[, c(1,2,3,4,12,14)]
df.base <- df.base[, -c(11,15)]

# -------------------------------------------------------------------------------------
# VSURF {VSURF} is used for variable selection running in parallel
# in order to select meaningful covariants....
# -------------------------------------------------------------------------------------
vbase <- VSURF(x=df.base[,c(-11)], # TEMP
               y=df.base[,c(11)],  # All other remaining variables
               ntree = 2000,
               na.action = na.omit,
               nfor.pred = 50,
               nsd = 2,
               parallel = TRUE,
               ncores = 8) # number of threads
# -------------------------------------------------------------------------------------

# A summary is requested
names(vbase)
summary(vbase)

#  Variables selected after "prediction step" .
#  This is the most important!!!
vbase$varselect.pred
# [1] 4; output suggest that [4], meaning Z, meaning ALOS_FILL is the ONLY
# relevant variables, other variables, even when CORR is high, are redundant

# Variables selected after "interpretation step"
vbase$varselect.interp
# [1] 4  8  9 11 10  3  2  7  1  6  5 12

# Variables selected after "thresholding step"
vbase$varselect.thres
# [1]  4  8  9 11 10  3  2  7  1  6  5 12

# VSURF object is plotted
plot(vbase)

# VSURF object is saved in binary R format
# save(vbase, file = "vbase.RData")

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Variable Selection Using {glmulti}
# For exploration purposes ONLY.... not to be used...yet...
# /////////////////////////////////////////////////////////////////////////////////////

# a lm {stats} model is created
mod.geospatial01 <- glm(TEMP ~ ., data = df.base, family = gaussian)

# A summary is requested over mod.geospatial01
summary(mod.geospatial01)

# glmulti {glmulti} function is applied to mod.geospatial01
best.geo <- glmulti(mod.geospatial01,
                    level =1,
                    crit = "aicc",
                    method = "g", # genetic algorithms
                    plotty = FALSE)

# An object summary is requested
summary(best.geo)

# weightable {glmulti} function is called
# The best model convination is requested
weightable(best.geo)
# TEMP ~ 1 + X + Y + SLOPE + ALOS_FILL + ASPECT + TOPOGRAPHIC + LS_FACTOR + DIRECT_INSO
# aicc      weights
# 1   2681.644 1.086850e-01
# As seen, if GLM models where to also be included in Machine Learning interpolation,
# the best model includes much more variables than RF, which only takes "Z"

# plot {graphics} function is called
plot(best.geo, type = "s")

# Best glm model is created
final.geo <- lm(TEMP ~ 1 + X + Y + SLOPE + ALOS_FILL + ASPECT + TOPOGRAPHIC + LS_FACTOR + DIRECT_INSO, 
                data = df.base)

# An object summary is requested
summary(final.geo)
# Adjusted R-squared:  0.9323; is very very high, but this is assuming linear relations
# among the covariants

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: External IMN Binary Temperature 
# /////////////////////////////////////////////////////////////////////////////////////

# Binary R objects are externally called
load("df.temperatureCR.RData")

# Imported data.frames are renamed. 
# It should be changed in the EXCEL_2018.R to stop using these names !!!
df.import <- df.temperatureCR

# 10-fold cross-validation is selected based on JAN ONLY !!!!
# The remainig months are not needed, this is just a fake-artefect
df.K <- subset(df.import, MONTH == "JAN")

# Folds are defined
folds <- cut(seq(1,nrow(df.K)),breaks=10,labels=FALSE)

# Seed is set
set.seed(def.seed)

# Folds vector is randomized
folds <- sample(folds)

# Folds vector is printed
print(folds)

# fold vector is created
k.folds <- rep(folds, 12)

# A k-fold variables is created
df.import$kfold <- k.folds

# data.frame object is saved as *.csv
write.csv(df.temperatureCR, file = "df.temperatureCR.csv")

# A "levels" monthly vector is created for visualization puposes
v.select.name <- c("JAN", "FEB", "MAR","APR", "MAY", "JUN", "JUL","AUG", "SET", "OCT", "NOV","DIC")

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: RASTER ANALYSIS. ILWIS-GIS. This block is very independent from the ramaining blocks
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Georeference attributes coming from ILWIS-GIS
xmin <- 286080.50
xmax <- 661080.50
ymin <- 888734.41
ymax <- 1241734.41
goef.col <- 375.00
geof.row <- 353.00
spa.res <- 1000.00

# An empty data.frame is created with a selected resolution
# The extent of this SpatialGrid is valid only for domain under analysis
georef <- expand.grid(X = seq(xmin, xmax, by = spa.res),
                      Y = seq(ymin, ymax, by = spa.res),
                      KEEP.OUT.ATTRS = F)

# The empty data.frame is converted into SpatialPoints
coordinates(georef) <- ~X + Y

# An empty SpatialGrid is created with a selected resolution
# The extent of this SpatialGrid is valid only for domain under analysis
basin.grid <- SpatialGrid(GridTopology(c(X = xmin, Y = ymin),
                                       c(spa.res, spa.res),
                                       c(goef.col, geof.row)))

# The empty SpatialGrid is converted into SpatialPoints
basin.grid <- SpatialPoints(basin.grid)

# The empty SpatialGrid is converted into SpatialPixels
gridded(basin.grid) <- T

#-------------------------------------------------------------------------------------
# The *.ASC map is imported from GIS to be used as a blank or mask-layer
# This must be done manually as every *.ASC is different:
#-------------------------------------------------------------------------------------
blank <- paste("mask_nacional_filter.asc", sep = "/")
#-------------------------------------------------------------------------------------

# The blank map is transformed into SpatialGridDataFrame
blank.grid <- read.asciigrid(blank, as.image = FALSE)

# The map is transformed into RasterLayer
blank.raster <- raster(blank.grid)

# The map is printed as verification
image(blank.raster, main = "*.ASC Blanking Map", col = topo.colors(64))

# KSE complementary data; Z, ASPECT, SLOPE. Prepared externally by IWILS-GIS
inter.spt <- read.table("SAGA_multiple_xyz_ALOS_1km.xyz", header = T) # "spatial.dat"

# A data.frame clone is created
inter.spt.clone <- inter.spt

# Check manually to ensure all variables are correct
summary(inter.spt)

# Empty KED complementary data SpatialGrid is converted into SpatialPixelsDataFrame
# This example has a 1000 x 1000 m spatial resolution. It varies with the waterbasin under analysis
gridded(inter.spt) = ~ X + Y

# The map is printed for verification purposes
plot(inter.spt)

# Containers lists objects are defined
final.list <- NULL
final.list <- list()

final.list.IDW <- NULL
final.list.IDW <- list()

final.list.cubist <- NULL
final.list.cubist <- list()

final.list.OK <- NULL
final.list.OK <- list()

final.list.RF <- NULL
final.list.RF <- list()

final.list.KED <- NULL
final.list.KED <- list()

depo.SELECTED.out <- NULL
depo.SELECTED.out <- list()

depo.KF1 <- NULL
depo.KF1 <- list()

depo.KF2 <- NULL
depo.KF2 <- list()

# Vector containers are defined and reset
DELTA.MAE.idw3 <- NULL
DELTA.MAE.idw3 <- vector()

DELTA.RMSE.idw3 <- NULL
DELTA.RMSE.idw3 <- vector()

DELTA.R2.idw3 <- NULL
DELTA.R2.idw3 <- vector()

DELTA.MAE.cubist <- NULL
DELTA.MAE.cubist <- vector()

DELTA.RMSE.cubist <- NULL
DELTA.RMSE.cubist <- vector()

DELTA.R2.cubist <- NULL
DELTA.R2.cubist <- vector()

DELTA.MAE.ok <- NULL
DELTA.MAE.ok <- vector()

DELTA.RMSE.ok <- NULL
DELTA.RMSE.ok <- vector()

DELTA.R2.ok <- NULL
DELTA.R2.ok <- vector()

DELTA.MAE.ked <- NULL
DELTA.MAE.ked <- vector()

DELTA.RMSE.ked <- NULL
DELTA.RMSE.ked <- vector()

DELTA.R2.ked <- NULL
DELTA.R2.ked <- vector()

DELTA.MAE.RF <- NULL
DELTA.MAE.RF <- vector()

DELTA.RMSE.RF <- NULL
DELTA.RMSE.RF <- vector()

DELTA.R2.RF <- NULL
DELTA.R2.RF <- vector()

R2.RF.train <- NULL
R2.RF.train <- vector()

error.RF.train <- NULL
error.RF.train <- vector()

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: INTERPOLATION; for all interpolation methods
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

# -------------------------------------------------------------------------------------
# Main list-indexing loop is initiated
# This loop is intended for k-fold validation
# -------------------------------------------------------------------------------------
for (u in 1 :10) {
  
  # K-fold sequence is defined
  k.selected <- u
  
  # -------------------------------------------------------------------------------------
  # Secondary list-indexing loop is initiated
  # This loop is intended for monthly analysis
  # -------------------------------------------------------------------------------------
  for (b in 1 :12) {
    
    # Vector containers are defined and reset
    amount.SHORT <- NULL
    amount.SHORT <- vector()
    amount.LARGE <- NULL
    amount.LARGE <- vector()
    
    # Vector containers are defined and reset
    df.SUB.LARGE <- NULL
    df.SUB.SHORT <- NULL
    lambda.SHORT <- NULL
    lambda.LARGE <- NULL
    df.SUBSET <- NULL
    df.SUBSET <- data.frame()
    df.EXPORT <- NULL
    df.EXPORT <- data.frame()
    
    # data.frame containers are defined and reset
    df.temp01 <- NULL
    df.taylor.IDW <- NULL
    df.taylor.cubist <- NULL
    df.taylor.OK <- NULL
    df.taylor.KED <- NULL
    df.taylor.RF <- NULL
    
    # -------------------------------------------------------------------------------------
    # Subsetting is performed based on attribute "MONTH"
    # This must be done manually as every *.ASC is different:
    # df.SUBSET is cloned for KED final interpolation
    # -------------------------------------------------------------------------------------
    df.SUBSET <- subset(df.import, ID == b)
    df.EXPORT <- subset(df.import, ID == b)
    # -------------------------------------------------------------------------------------
    
    # Stations above 1900 masl are preserved
    df.SUBSET.Z <- subset(df.SUBSET, Z  > 1900)
    df.SUBSET <- subset(df.SUBSET, Z  < 1900)
    df.k2 <- df.SUBSET
    
    # -------------------------------------------------------------------------------------
    # Fixing Sample Sizes
    # -------------------------------------------------------------------------------------
    
    # 75% of the sample size is isolated
    #smp_size <- floor(0.75 * nrow(df.SUBSET))
    
    # If HTR sub-network is below 4, HSR must be equal to HSR as the minimum is 4 stations
    #if (smp_size < 4) {
    #  smp_size <- 3
    #}
    
    # A random training sample is defined
    #train_ind <- sample(seq_len(nrow(df.SUBSET)), size = smp_size)
    
    # Training and Validation datasets are defined
    #df.SUB.LARGE <- df.SUBSET[train_ind, ]  # Training
    #df.SUB.SHORT <- df.SUBSET[-train_ind, ] # Validation
    
    # Training and Validation datasets are defined
    df.SUB.LARGE <- subset(df.SUBSET, kfold != k.selected)  # Training
    df.SUB.SHORT <- subset(df.SUBSET, kfold == k.selected) # Validation
    
    # Stations above 1900 masl are added
    df.SUB.LARGE <- rbind(df.SUB.LARGE, df.SUBSET.Z)
    
    # -------------------------------------------------------------------------------------
    # INTERPOLATION data set *********************************************
    # -------------------------------------------------------------------------------------
    
    # ONLY relevant variables are selected for variogram fitting
    df.LARGE.K <- df.SUB.LARGE # NOT transformed
    
    # Coordinates are added
    coordinates(df.LARGE.K) = ~ X + Y
    
    # -------------------------------------------------------------------------------------
    # SUB-BLOCK: IDW3
    # -------------------------------------------------------------------------------------
    
    # Total Interpolation Methods
    idw3.LARGE <- gstat::idw(TEMP ~ 1,
                             loc = df.LARGE.K,
                             newdata = inter.spt,
                             idp = 3,
                             na.action = na.pass)  # Inverse Distance Weighting Exp 3
    
    # SpatialPixelsDataFrames are converted into RasterLayers
    idw3.all.raster <- raster(idw3.LARGE)
    
    # The map is printed as verification
    #plot(idw3.all.raster + blank.raster*0)
    
    # -------------------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------------------
    # SUB-BLOCK: Auto-OK
    # -------------------------------------------------------------------------------------
    
    # Total Interpolation Methods
    ok.LARGE <- autoKrige(formula = (TEMP^1) ~ 1,
                          input_data = df.LARGE.K,
                          #block = c(10000, 10000),
                          model = c("Sph", "Exp", "Gau", "Mat", "Ste"),
                          #model = c("Exp","Sph", "Gau", "Exc", "Mat", "Ste", "Cir", "Lin", "Bes","Pen", "Per", "Wav","Hol", "Log", "Pow","Spl","Leg"),
                          verbose = FALSE,
                          fix.values = c(NA,NA,NA),
                          start_vals = c(NA,NA,NA),
                          miscFitOptions = list(min.np.bin = 8, merge.small.bins = TRUE),
                          new_data = inter.spt) 
    
    # SpatialPixelsDataFrames are converted into RasterLayers
    ok.all.raster <- raster(ok.LARGE[[1]][1])
    
    # The units are back-transformed
    ok.all.raster <- (ok.all.raster)^(1/1)
    
    # The map is printed as verification
    #plot(ok.all.raster + blank.raster*0)
    
    # -------------------------------------------------------------------------------------
    # SUB-BLOCK: Auto-KED
    # -------------------------------------------------------------------------------------
    
    # Total Interpolation Methods
    ked.LARGE <- autoKrige(formula = (TEMP^1) ~ Z,
                           input_data = df.LARGE.K,
                           #block = c(10000, 10000),
                           model = c("Sph", "Exp", "Gau", "Mat", "Ste"),
                           #model = c("Exp","Sph", "Gau", "Exc", "Mat", "Ste", "Cir", "Lin", "Bes","Pen", "Per", "Wav","Hol", "Log", "Pow","Spl","Leg"),
                           verbose = FALSE,
                           fix.values = c(NA,NA,NA),
                           start_vals = c(NA,NA,NA),
                           miscFitOptions = list(min.np.bin = 8, merge.small.bins = TRUE),
                           new_data = inter.spt) 
    
    # SpatialPixelsDataFrames are converted into RasterLayers
    ked.all.raster <- raster(ked.LARGE[[1]][1])
    
    # The units are back-transformed
    ked.all.raster <- (ked.all.raster)^(1/1)
    
    # The map is printed as verification
    #plot(ked.all.raster + blank.raster*0)
    
    # -------------------------------------------------------------------------------------
    # SUB-BLOCK: Random Forest
    # -------------------------------------------------------------------------------------
    
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # RF: Spatial prediction 2D continuous variable using buffer distances
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    # A buffer distance is derived and a SpatialPixelsDataFrame is created
    df.LARGE.temp <- df.LARGE.K[1,]
    
    grid.dist0.CR <- GSIF::buffer.dist(df.LARGE.temp["TEMP"], # "SpatialPointsDataFrame"
                                       inter.spt,#inter.spt[1], # "SpatialPixelsDataFrame" (coordinates)
                                       as.factor(1:nrow(df.LARGE.temp))) # SpatialPointsDataFrame (number of variables)
    
    #A character vector is created
    #dn0.CR <- paste(names(grid.dist0.CR), collapse="+")
    
    # A formula is created for TEMP ONLY !!!
    #fm0.CR <- as.formula(paste("TEMP ~ ", dn0.CR))
    
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # RF: Spatial prediction 2D variable with covariates
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    # A buffer distance is derived and a SpatialPixelsDataFrame is created
    # grids.spc <- GSIF::spc(inter.spt, as.formula("~ 
    #                                  Slope +
    #                                  Z +
    #                                  Aspect +
    #                                  TWI +
    #                                  LS_Factor +
    #                                  Direct_Insolation +
    #                                  Convexity +
    #                                  Wind_Exposition +
    #                                  XX +
    #                                  YY")) 
    
    grids.spc <- GSIF::spc(inter.spt, 
                           as.formula("~ Z + XX + YY"))
    
    # A formula is created for TEMP
    #fm1 <- as.formula(paste("TEMP ~ ", dn0.CR, " + ", 
    #                        paste(names(grids.spc@predicted), collapse = "+")))
    
    # An alternative formula is created using PCA components ONLY
    # NO buffer distances included
    fm1 <- as.formula(paste("TEMP ~ ", 
                            paste(names(grids.spc@predicted), collapse = "+")))
    
    # {sp} SpatialPointsDataFrame and SpatialPixelsDataFrame are overlayed
    # and a new data.frame is created
    # ov.BASE_P <- over(df.LARGE.K["TEMP"], grid.dist0.CR)
    ov.BASE_P <- df.LARGE.K["TEMP"]
    
    # Two data.frames are cbinded and a new data.frame is created
    rm.BASE_P <- cbind(df.LARGE.K@data["TEMP"], ov.BASE_P)
    
    # {sp} SpatialPointsDataFrame and SpatialPixelsDataFrame are overlayed
    # and a new data.frame is created
    ov.BASE_P1 <- over(df.LARGE.K["TEMP"], grids.spc@predicted)
    
    # Two data.frames are cbinded and a new data.frame is created
    rm.BASE_P1 <- do.call(cbind,
                          list(df.LARGE.K@data["TEMP"],
                               #ov.BASE_P,
                               ov.BASE_P1))
    
    # raster dimensionality is requested
    print(nrow(rm.BASE_P1))
    
    # ranger {ranger} is used for random forests
    m.BASE_P1 <- ranger(fm1,           # formula based on number of variables
                        rm.BASE_P1,    # data.frame
                        quantreg = TRUE,
                        num.trees = 500,
                        importance = "impurity")
    
    # Variable importance is requested
    m.BASE_P1$variable.importance
    
    # A ranger object is requested
    print(m.BASE_P1)
    
    # Ranger object statistics are saved
    R2.RF.train <- m.BASE_P1$r.squared
    error.RF.train <- m.BASE_P1$prediction.error
    
    # predict is used to execute a prediction, ranger vs. data.frame
    # and a new list is created
    grid.dist0.CR1 <- grid.dist0.CR
    
    # PCA are added to SpatialPixelsDataFrame
    grid.dist0.CR1$PC1 <- grids.spc@predicted@data$PC1
    grid.dist0.CR1$PC2 <- grids.spc@predicted@data$PC2
    grid.dist0.CR1$PC3 <- grids.spc@predicted@data$PC3
    #grid.dist0.CR1$PC4 <- grids.spc@predicted@data$PC4
    #grid.dist0.CR1$PC5 <- grids.spc@predicted@data$PC5
    #grid.dist0.CR1$PC6 <- grids.spc@predicted@data$PC6
    #grid.dist0.CR1$PC7 <- grids.spc@predicted@data$PC7
    #grid.dist0.CR1$PC8 <- grids.spc@predicted@data$PC8
    #grid.dist0.CR1$PC9 <- grids.spc@predicted@data$PC9
    #grid.dist0.CR1$PC10 <- grids.spc@predicted@data$PC10
    
    # predict is used to execute a prediction, ranger vs. data.frame
    # and a new list is created
    BASE_P.rfd1 <- predict(m.BASE_P1, grid.dist0.CR1@data)
    
    # A new prediction varible is created within the SpatialPixelsDataFrame
    inter.spt$temp_rfd1 <- BASE_P.rfd1[[1]]
    #meuse.grid$zinc_rfd_range <- (zinc.rfd[,3]-zinc.rfd[,1])/2
    
    # A new prediction raster is created "[represents the position of temp_rfd1 !!!!!!!!]"
    rf.all.raster <- raster(inter.spt[,11])
    
    # The map is printed as verification
    #plot(rf.all.raster + blank.raster*0)
    
    # -------------------------------------------------------------------------------------
    # SUB-BLOCK: Cubist
    # -------------------------------------------------------------------------------------
    
    # cubist.default {Cubist} is used to fit a Cubist model
    mod.cubist01 <- cubist(x = rm.BASE_P1[, -1], y = rm.BASE_P1$TEMP)
    
    # Cubist model summary is requested
    print(summary(mod.cubist01))
    
    # A spatial cubist prediction is executes
    spa.cubist01 <- predict(mod.cubist01, grid.dist0.CR1@data, neighbors = 5)
    
    # A new prediction varible is created within the SpatialPixelsDataFrame
    inter.spt$CUBIST <- spa.cubist01
    
    # SpatialPixelsDataFrames are converted into RasterLayers
    cubist.all.raster <- raster(inter.spt[,12])
    
    # The map is printed as verification
    # plot(cubist.all.raster + blank.raster*0)
    
    # -------------------------------------------------------------------------------------
    # SUB-BLOCK: Generation of final GeoTIFF Export
    # -------------------------------------------------------------------------------------
    
    # Coordinates are added
    coordinates(df.EXPORT) = ~ X + Y
    
    # ****************************************************
    # Final interpolation using KED ---- option 01
    # ****************************************************
    ked.LARGE.EXPORT <- autoKrige(formula = (TEMP^1) ~ Z,
                                  input_data = df.EXPORT,
                                  #block = c(10000, 10000),
                                  model = c("Sph", "Exp", "Gau", "Mat", "Ste"),
                                  #model = c("Exp","Sph", "Gau", "Exc", "Mat", "Ste", "Cir", "Lin", "Bes","Pen", "Per", "Wav","Hol", "Log", "Pow","Spl","Leg"),
                                  verbose = FALSE,
                                  fix.values = c(NA,NA,NA),
                                  start_vals = c(NA,NA,NA),
                                  miscFitOptions = list(min.np.bin = 8, merge.small.bins = TRUE),
                                  new_data = inter.spt) 
    
    # SpatialPixelsDataFrames are converted into RasterLayers
    ked.all.raster.export <- raster(ked.LARGE.EXPORT[[1]][1])
    
    # The units are back-transformed
    ked.all.raster.export <- (ked.all.raster.export)^(1/1)
    
    # Raster objects are blanked
    ked.map.export <- ked.all.raster.export + (blank.raster*0)

    # -------------------------------------------------------------------------------------
    # Raster objects are individually saved within a list container
    # Operator MUST choose between KED, RF o CUBIST to fillup depo.SELECTED.out
    depo.SELECTED.out[[b]] <- ked.map.export
    # -------------------------------------------------------------------------------------
    
    # The map is printed as verification
    # plot(ked.map.export)
    
    # -------------------------------------------------------------------------------------
    # VALIDATION data set *********************************************
    # -------------------------------------------------------------------------------------
    
    # ONLY relevant variables are selected for variogram fitting
    df.SHORT.K <- df.SUB.SHORT # NOT transformed
    
    # Coordinates are added
    coordinates(df.SHORT.K) = ~ X + Y
    
    # The grid is converted into RasterLayer
    raster.georef <- raster(basin.grid)  
    
    # BASE_P field is extracted
    ext.raster <- df.SHORT.K[2]  
    
    # Point values are rasterized
    point.raster <- rasterize(ext.raster, idw3.all.raster) # raster.georef)  
    
    # The layers are disaggregated from the point.raster
    upper.layer <- unstack(point.raster)
    lower.layer <- upper.layer[[2]]
    
    # ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # The R2 is calculated
    # ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    # ---------------------------------------------
    # lower.layer[] comparison layer is extracted
    v.lower <- as.vector(lower.layer[])
    v.lower <- v.lower[!is.na(v.lower)]
    # ---------------------------------------------
    
    df.R2 <- NULL
    v.idw3 <- (idw3.all.raster) - (lower.layer*0)
    v.idw3 <- as.vector(v.idw3)
    v.idw3 <- v.idw3[!is.na(v.idw3)]
    shortest <- min(length(v.idw3), length(v.lower))
    y.SEQ <- tail(v.lower, shortest)
    x.SEQ <- tail(v.idw3, shortest)
    df.R2 <- data.frame(x.SEQ, y.SEQ)
    names(df.R2) <- c("X","Y")
    R2.idw3 <- round((summary(lm(Y ~ X, data=df.R2))$r.squared),4)
    data.OBS.IDW <- y.SEQ
    data.MOD.IDW <- x.SEQ
    df.taylor.IDW <- data.frame(data.OBS.IDW, data.MOD.IDW, rep(b,length(x.SEQ)), c("IDW"), rep(u,length(x.SEQ)))
    names(df.taylor.IDW) <- c("OBS", "MOD", "MONTH_NUM", "METHOD","KF")
    
    df.R2 <- NULL
    v.cubist <- (cubist.all.raster) - (lower.layer*0)
    v.cubist <- as.vector(v.cubist)
    v.cubist <- v.cubist[!is.na(v.cubist)]
    shortest <- min(length(v.cubist), length(v.lower))
    y.SEQ <- tail(v.lower, shortest)
    x.SEQ <- tail(v.cubist, shortest)
    df.R2 <- data.frame(x.SEQ, y.SEQ)
    names(df.R2) <- c("X","Y")
    R2.cubist <- round((summary(lm(Y ~ X, data=df.R2))$r.squared),4)
    data.OBS.cubist <- y.SEQ
    data.MOD.cubist <- x.SEQ
    df.taylor.cubist <- data.frame(data.OBS.cubist, data.MOD.cubist, rep(b,length(x.SEQ)), c("cubist"),rep(u,length(x.SEQ)))
    names(df.taylor.cubist) <- c("OBS", "MOD", "MONTH_NUM", "METHOD","KF")
    
    df.R2 <- NULL
    v.ok.all.raster <- (ok.all.raster) - (lower.layer*0)
    v.ok.all.raster <- as.vector(v.ok.all.raster)
    v.ok.all.raster <- v.ok.all.raster[!is.na(v.ok.all.raster)]
    shortest <- min(length(v.ok.all.raster), length(v.lower))
    y.SEQ <- tail(v.lower, shortest)
    x.SEQ <- tail(v.ok.all.raster, shortest)
    df.R2 <- data.frame(x.SEQ, y.SEQ)
    names(df.R2) <- c("X","Y")
    R2.ok <- round((summary(lm(Y ~ X, data=df.R2))$r.squared),4)
    data.OBS.OK <- y.SEQ
    data.MOD.OK <- x.SEQ
    df.taylor.OK <- data.frame(data.OBS.OK, data.MOD.OK, rep(b,length(x.SEQ)), c("OK"),rep(u,length(x.SEQ)))
    names(df.taylor.OK) <- c("OBS", "MOD", "MONTH_NUM", "METHOD","KF")
    
    df.R2 <- NULL
    v.ked.all.raster <- (ked.all.raster) - (lower.layer*0)
    v.ked.all.raster <- as.vector(v.ked.all.raster)
    v.ked.all.raster <- v.ked.all.raster[!is.na(v.ked.all.raster)]
    shortest <- min(length(v.ked.all.raster), length(v.lower))
    y.SEQ <- tail(v.lower, shortest)
    x.SEQ <- tail(v.ked.all.raster, shortest)
    df.R2 <- data.frame(x.SEQ, y.SEQ)
    names(df.R2) <- c("X","Y")
    R2.ked <- round((summary(lm(Y ~ X, data=df.R2))$r.squared),4)
    data.OBS.KED <- y.SEQ
    data.MOD.KED <- x.SEQ
    df.taylor.KED <- data.frame(data.OBS.KED, data.MOD.KED, rep(b,length(x.SEQ)), c("KED"),rep(u,length(x.SEQ)))
    names(df.taylor.KED) <- c("OBS", "MOD", "MONTH_NUM", "METHOD","KF")
    
    df.R2 <- NULL
    v.rf.all.raster <- (rf.all.raster) - (lower.layer*0)
    v.rf.all.raster <- as.vector(v.rf.all.raster)
    v.rf.all.raster <- v.rf.all.raster[!is.na(v.rf.all.raster)]
    shortest <- min(length(v.rf.all.raster), length(v.lower))
    y.SEQ <- tail(v.lower, shortest)
    x.SEQ <- tail(v.rf.all.raster, shortest)
    df.R2 <- data.frame(x.SEQ, y.SEQ)
    names(df.R2) <- c("X","Y")
    R2.rf <- round((summary(lm(Y ~ X, data=df.R2))$r.squared),4)
    data.OBS.RF <- y.SEQ
    data.MOD.RF <- x.SEQ
    df.taylor.RF <- data.frame(data.OBS.RF, data.MOD.RF, rep(b,length(x.SEQ)), c("RF"),rep(u,length(x.SEQ)))
    names(df.taylor.RF) <- c("OBS", "MOD", "MONTH_NUM", "METHOD","KF")
    
    # ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # The Mean Absolute Error is calculated
    # ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    abs.idw3 <- abs(lower.layer - (idw3.all.raster))
    sum.idw3 <- cellStats(abs.idw3, (sum))
    n.idw3 <- sum(!is.na(abs.idw3[]))
    MAE.idw3 <- sum.idw3/n.idw3
    
    abs.cubist <- abs(lower.layer - (cubist.all.raster))
    sum.cubist <- cellStats(abs.cubist, (sum))
    n.cubist <- sum(!is.na(abs.cubist[]))
    MAE.cubist <- sum.cubist/n.cubist
    
    abs.ok <- abs(lower.layer - (ok.all.raster))
    sum.ok <- cellStats(abs.ok, (sum))
    n.ok <- sum(!is.na(abs.ok[]))
    MAE.ok <- sum.ok/n.ok
    
    abs.ked <- abs(lower.layer - (ked.all.raster))
    sum.ked <- cellStats(abs.ked, (sum))
    n.ked <- sum(!is.na(abs.ked[]))
    MAE.ked <- sum.ked/n.ked
    
    abs.rf <- abs(lower.layer - (rf.all.raster))
    sum.rf <- cellStats(abs.rf, (sum))
    n.rf <- sum(!is.na(abs.rf[]))
    MAE.rf <- sum.rf/n.rf
    
    # ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # The Root Mean Square Error is calculated
    # ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    square.idw3 <- ((lower.layer - (idw3.all.raster)) ^ 2)
    sum.idw3 <- cellStats(square.idw3, (sum))
    n.idw3 <- sum(!is.na(abs.idw3[]))
    RMSE.idw3 <- (sum.idw3/n.idw3)^(0.5)
    
    square.cubist <- ((lower.layer - (cubist.all.raster)) ^ 2)
    sum.cubist <- cellStats(square.cubist, (sum))
    n.cubist <- sum(!is.na(abs.cubist[]))
    RMSE.cubist <- (sum.cubist/n.cubist)^(0.5)
    
    square.ok <- ((lower.layer - (ok.all.raster)) ^ 2)
    sum.ok <- cellStats(square.ok, (sum))
    n.ok <- sum(!is.na(abs.ok[]))
    RMSE.ok <- (sum.ok/n.ok)^(0.5)
    
    square.ked <- ((lower.layer - (ked.all.raster)) ^ 2)
    sum.ked <- cellStats(square.ked, (sum))
    n.ked <- sum(!is.na(abs.ked[]))
    RMSE.ked <- (sum.ked/n.ked)^(0.5)
    
    square.rf <- ((lower.layer - (rf.all.raster)) ^ 2)
    sum.rf <- cellStats(square.rf, (sum))
    n.rf <- sum(!is.na(abs.rf[]))
    RMSE.rf <- (sum.rf/n.rf)^(0.5)
    
    # Vector containers are defined and reset
    DELTA.MAE.idw3 <- MAE.idw3
    DELTA.RMSE.idw3 <- RMSE.idw3
    DELTA.R2.idw3 <- R2.idw3
    
    DELTA.MAE.cubist <- MAE.cubist
    DELTA.RMSE.cubist <- RMSE.cubist
    DELTA.R2.cubist <- R2.cubist
    
    DELTA.MAE.ok <- MAE.ok
    DELTA.RMSE.ok <- RMSE.ok
    DELTA.R2.ok <- R2.ok
    
    DELTA.MAE.ked <- MAE.ked
    DELTA.RMSE.ked <- RMSE.ked
    DELTA.R2.ked <- R2.ked
    
    DELTA.MAE.rf <- MAE.rf
    DELTA.RMSE.rf <- RMSE.rf
    DELTA.R2.rf <- R2.rf
    
    # KED OFs are printed for verification
    print(R2.ked)
    print(R2.rf)
    print(R2.cubist)
    print(u)
    
    # Loop year is printed
    print(unique(df.SUB.LARGE$MONTH))
    
    # A month numerical values is added
    num.month <- b
    num.kf <- u
    
    # A temporal data.frame is created
    df.temp01 <- data.frame(num.month,
                            DELTA.MAE.idw3, DELTA.RMSE.idw3, DELTA.R2.idw3,
                            DELTA.MAE.cubist, DELTA.RMSE.cubist, DELTA.R2.cubist, 
                            DELTA.MAE.ok, DELTA.RMSE.ok, DELTA.R2.ok, 
                            DELTA.MAE.ked, DELTA.RMSE.ked, DELTA.R2.ked,
                            DELTA.MAE.rf, DELTA.RMSE.rf, DELTA.R2.rf,
                            R2.RF.train, error.RF.train,num.kf)
    
    # List containers are filled
    final.list[[b]] <- df.temp01
    final.list.IDW[[b]] <- df.taylor.IDW
    final.list.cubist[[b]] <- df.taylor.cubist
    final.list.OK[[b]] <- df.taylor.OK
    final.list.KED[[b]] <- df.taylor.KED
    final.list.RF[[b]] <- df.taylor.RF
    
  } # SECONDARY LOOP IS CLOSED
  
  # List containing data frames are converted into one unique data.frame
  df.final <- ldply(final.list, data.frame)
  df.final.IDW <- ldply(final.list.IDW, data.frame)
  df.final.cubist <- ldply(final.list.cubist, data.frame)
  df.final.OK <- ldply(final.list.OK, data.frame)
  df.final.KED <- ldply(final.list.KED, data.frame)
  df.final.RF <- ldply(final.list.RF, data.frame)
  
  # A compiled Obs vs. Mod data.frame is created 
  df.compile.comp <- rbind(df.final.IDW,
                           df.final.cubist,
                           df.final.OK,
                           df.final.KED,
                           df.final.RF)
  
  # Monthly variable is transformed to character
  df.compile.comp$MONTH_NUM <- as.character(df.compile.comp$MONTH_NUM)
  
  # Character labels are included
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "1"] <- "JAN"
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "2"] <- "FEB"
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "3"] <- "MAR"
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "4"] <- "APR"
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "5"] <- "MAY"
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "6"] <- "JUN"
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "7"] <- "JUL"
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "8"] <- "AUG"
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "9"] <- "SEP"
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "10"] <- "OCT"
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "11"] <- "NOV"
  df.compile.comp$MONTH_NUM[df.compile.comp$MONTH_NUM == "12"] <- "DIC"
  
  # Monthly variable is transformed to factor
  df.compile.comp$MONTH_NUM <- as.factor(df.compile.comp$MONTH_NUM)
  
  # A "levels" monthly vector is created for visualization puposes
  v.select.name <- c("JAN", "FEB", "MAR","APR", "MAY", "JUN", "JUL","AUG", "SEP", "OCT", "NOV","DIC")
  
  # Months are ordered to avoid ggplot plotting conflicts
  df.compile.comp$Month.ordered <- ordered(df.compile.comp$MONTH_NUM, v.select.name)
  
  # data.frames are saved based on temperature ONLY
  df.final.TEMP <- df.final
  df.compile.comp.TEMP <- df.compile.comp
  
  # List containers are filled
  depo.KF1[[u]] <- df.final.TEMP
  depo.KF2[[u]] <- df.compile.comp.TEMP
  
} # MAIN LOOP IS CLOSED

# data.frames are saved in binary R format
#save(depo.KF1, file = "depo.KF1.RData")
#save(depo.KF2, file = "depo.KF2.RData")

# List containing data frames are converted into one unique data.frame
df.final.TEMP <- ldply(depo.KF1, data.frame)
df.compile.comp.TEMP <- ldply(depo.KF2, data.frame)

# data.frames are saved in binary R format
#save(df.final.TEMP, file = "df.final.TEMP.RData")
#save(df.compile.comp.TEMP, file = "df.compile.comp.TEMP.RData")

# data.frame object is saved as *.csv
#write.csv(df.final.TEMP, file = "df_final_TEMP_K0.csv")
#write.csv(df.compile.comp.TEMP, file = "df_compile_comp_TEMP_K0.csv")

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Post-Processing
# /////////////////////////////////////////////////////////////////////////////////////

load("df.final.TEMP.RData")
load("df.compile.comp.TEMP.RData")

# cubist and IDW3 are deleted from the selection
df.compile.comp.TEMP <- subset(df.compile.comp.TEMP, METHOD != "IDW")

# A "levels" monthly vector is created for visualization puposes
v.select.name <- c("JAN", "FEB", "MAR","APR", "MAY", "JUN", "JUL","AUG", "SEP", "OCT", "NOV","DIC")

# Months are ordered to avoid ggplot plotting conflicts
df.compile.comp.TEMP$Month.ordered <- ordered(df.compile.comp.TEMP$Month.ordered, v.select.name)

# A new scale color is defined
scale01 <- c("#d95f02", "#1f78b4")
scale02 <- c("#d95f02", "#1f78b4", "#1b9e77")
scale03 <- c("#404040", "#ca0020", "#f4a582", "#bababa")

# A temperature Taylor diagram is created
TaylorDiagram(df.compile.comp.TEMP,
              obs = "OBS",
              mod = "MOD",
              group = "METHOD",
              type = "Month.ordered",
              normalise = FALSE,
              cols = scale03,
              fontsize= 12)

# A normalized Taylor diagram is created
TaylorDiagram(df.compile.comp.TEMP,
              obs = "OBS",
              mod = "MOD",
              group = "METHOD",
              type = "Month.ordered",
              normalise = TRUE,
              cols = scale03,
              fontsize= 14)

# Only relevant variables are selected
df.compile.TEMP.of <- df.final.TEMP[, c(1,17)]

# -------------------------------------------------------------------------------------
# ggplot graphics
# -------------------------------------------------------------------------------------

# A labelling variable is defined
df.final.TEMP$MONTH_CHR <- c("TEXT")

# Labelling loop is initialized
for (i in 1 : length(df.final.TEMP$num.month)) {
  
  if (df.final.TEMP$num.month[i] == 1) { df.final.TEMP$MONTH_CHR[i] <- c("JAN") } else if (df.final.TEMP$num.month[i] == 2) {
    df.final.TEMP$MONTH_CHR[i] <- c("FEB") } else if (df.final.TEMP$num.month[i] == 3) {
      df.final.TEMP$MONTH_CHR[i] <- c("MAR") } else if (df.final.TEMP$num.month[i] == 4) {
        df.final.TEMP$MONTH_CHR[i] <- c("APR") } else if (df.final.TEMP$num.month[i] == 5) {
          df.final.TEMP$MONTH_CHR[i] <- c("MAY") } else if (df.final.TEMP$num.month[i] == 6) {
            df.final.TEMP$MONTH_CHR[i] <- c("JUN") } else if (df.final.TEMP$num.month[i] == 7) {
              df.final.TEMP$MONTH_CHR[i] <- c("JUL") } else if (df.final.TEMP$num.month[i] == 8) {
                df.final.TEMP$MONTH_CHR[i] <- c("AUG") } else if (df.final.TEMP$num.month[i] == 9) {
                  df.final.TEMP$MONTH_CHR[i] <- c("SEP") } else if (df.final.TEMP$num.month[i] == 10) {
                    df.final.TEMP$MONTH_CHR[i] <- c("OCT") } else if (df.final.TEMP$num.month[i] == 11) {
                      df.final.TEMP$MONTH_CHR[i] <- c("NOV") } else if (df.final.TEMP$num.month[i] == 12) {
                        df.final.TEMP$MONTH_CHR[i] <- c("DIC") }
  
} # Labelling loop is closed

# Months are ordered to avoid ggplot plotting conflicts
df.final.TEMP$MONTH_CHR <- ordered(df.final.TEMP$MONTH_CHR, v.select.name)

# A RMSE plot is prepared
df.final.TEMP.RMSE <- df.final.TEMP [, c(6,9,12,15,20)]
df.final.TEMP.RMSE <- melt(df.final.TEMP.RMSE)

ggplot() +
  geom_boxplot(aes(y = value,x = variable,colour = variable),data=df.final.TEMP.RMSE,size = 0.75) +
  facet_wrap(facets = ~MONTH_CHR, ncol = 3) +
  theme_bw() +
  theme(axis.text = element_text(angle = 90.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(min.n = 5.0),limits = c(0,3.5))

# A MAE plot is prepared
df.final.TEMP.MAE <- df.final.TEMP [, c(5,8,11,14,20)]
df.final.TEMP.MAE <- melt(df.final.TEMP.MAE)

ggplot() +
  geom_boxplot(aes(y = value,x = variable,colour = variable),data=df.final.TEMP.MAE,size = 0.75) +
  facet_wrap(facets = ~MONTH_CHR, ncol = 3) +
  theme_bw() +
  theme(axis.text = element_text(angle = 90.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(min.n = 5.0),limits = c(0,4))

# A R2 plot is prepared
df.final.TEMP.R2 <- df.final.TEMP [, c(7,10,13,16,20)]
df.final.TEMP.R2 <- melt(df.final.TEMP.R2)

ggplot() +
  geom_boxplot(aes(y = value,x = variable,colour = variable),data=df.final.TEMP.R2,size = 0.75) +
  facet_wrap(facets = ~MONTH_CHR, ncol = 3) +
  theme_bw() +
  theme(axis.text = element_text(angle = 90.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(min.n = 5.0),limits = c(0.4,1))

# A RF plot is prepared
df.final.TEMP.RF <- df.final.TEMP [, c(17,20)]
df.final.TEMP.RF <- melt(df.final.TEMP.RF)

ggplot() +
  geom_boxplot(aes(y = value,x = variable),data=df.final.TEMP.RF,size = 0.75) +
  facet_wrap(facets = ~MONTH_CHR, ncol = 3) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  theme_bw()

# A R2 (OOB) simple boxplot is created
ggplot() +
  geom_boxplot(aes(y = R2.RF.train,x = as.factor(num.month), color=as.factor(num.month)),
               data=df.compile.TEMP.of, size = 1.25,outlier.colour = "black",
               outlier.shape = 8,outlier.size = 2.5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0))+
  #facet_wrap(facets = ~REGION) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90.0),
        text=element_text(size=20,  family="serif")) +
  xlab("Month of the year") +
  ylab("Training RF [R2 (OOB)]") 
#scale_color_manual(values=scale01)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Post-Processing GeoTIFF export
# /////////////////////////////////////////////////////////////////////////////////////

# A raster stack object is created
temp1x1km <- stack(depo.SELECTED.out)

# A simple plot is created on raster level [1]
plot(temp1x1km[[1]])

# Dummy NAs values are assigned
NAvalue(temp1x1km) <- (-9999)

# CRTM05 projection for CR is created
CRTM05 <- CRS("+proj=tmerc +lat_0=0 +lon_0=-84 +k=0.9999 +x_0=500000 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# CRTM05 projection is assigned to raster layer
crs(temp1x1km) <- CRTM05

# Rasterbrick object is exported to GeoTIFF format
writeRaster(temp1x1km, "TEMP_1960_1990_MONTHLY.tif", format="GTiff", overwrite = FALSE)

# CR climati Region ShapeFile is loaded
shp <- shapefile("zonas_climaticas_def_crtm05")

# CRTM05 projection is assigned to rasterbrick object
crs(shp) <- CRTM05

# A "levels" monthly vector is created for visualization puposes
v.select.name <- c("JAN", "FEB", "MAR","APR", "MAY", "JUN", "JUL","AUG", "SEP", "OCT", "NOV","DEC")

# Rasterbrick layer names are replaced
names(temp1x1km) <- v.select.name

# A color palette is created based on http://colorbrewer2.org/
colr.percent <- colorRampPalette(c("#006837",
                                   "#1a9850",
                                   "#66bd63",
                                   "#a6d96a",
                                   "#d9ef8b",
                                   "#fee08b",
                                   "#fdae61",
                                   "#f46d43",
                                   "#d73027",
                                   "#a50026"))

# A color palette is created based on http://colorbrewer2.org/
colr.rain <- colorRampPalette(c("#fff7fb",
                                "#ece7f2",
                                "#d0d1e6",
                                "#a6bddb",
                                "#74a9cf",
                                "#3690c0",
                                "#0570b0",
                                "#045a8d",
                                "#023858"))

# A color palette is created based on http://colorbrewer2.org/
colr.rain2 <- colorRampPalette(c("#ffffd9",
                                 "#edf8b1",
                                 "#c7e9b4",
                                 "#7fcdbb",
                                 "#41b6c4",
                                 "#1d91c0",
                                 "#225ea8",
                                 "#253494",
                                 "#081d58"))

# A color palette is created based on http://colorbrewer2.org/
colr.rain3 <- colorRampPalette(c("#053061",
                                 "#2166ac",
                                 "#4393c3",
                                 "#92c5de",
                                 "#d1e5f0",
                                 "#fddbc7",
                                 "#f4a582",
                                 "#d6604d",
                                 "#b2182b",
                                 "#67001f"))

# Color scale breaks are defined
brk <- do.breaks(c(0, 40), 10)
brk

# Color scale breaks are defined
brk2 <- do.breaks(c(0, 40), 10)
brk2

# Mean precipitation plot is created
levelplot(temp1x1km, layout=c(3, 4),
          col.regions = colr.rain3,
          #at = brk,
          colorkey = list(col = colr.rain3, at = brk)) +
  latticeExtra::layer(sp.polygons(shp, lwd = 1))

# ---------------------------------------------------------
# CV (SD) precipitation plot is created
#levelplot(ratio.RASTER, layout=c(4, 3),
#          col.regions = colr.percent,
#          at = brk2,
#          colorkey = list(col = colr.percent, at = brk2)) +
#  latticeExtra::layer(sp.polygons(shp, lwd = 3))
# ---------------------------------------------------------


#-------------------------------------------------------------------------------------
# References
#-------------------------------------------------------------------------------------

#Maximum Triangle Slope - Tarboton, D.G. (1997):'A new method for the determination of flow directions and upslope areas in grid digital elevation models',
#Water Resources Research, Vol.33, No.2, p.309-319

#T. Tadono, H. Nagai, H. Ishida, F. Oda, S. Naito, K. Minakawa, and H. Iwamoto. GENERATION OF THE 30 M-MESH GLOBAL DIGITAL SURFACE MODEL BY ALOS PRISM. ISPRS - International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, XLI-B4:157-162, 6 2016. ISSN 2194-9034. doi: 10.5194/isprs-archives-XLI-B4-157-2016. URL http://www.int-arch-photogramm-remote-sens-spatial-inf-sci.net/XLI-B4/157/2016/

#Beven, K.J., Kirkby, M.J. (1979): A physically-based variable contributing area model of basin hydrology' Hydrology Science Bulletin 24(1), p.43-69

#Boehner, J., Selige, T. (2006): Spatial Prediction of Soil Attributes Using Terrain Analysis and Climate Regionalisation'
#In: Boehner, J., McCloy, K.R., Strobl, J.: 'SAGA - Analysis and Modelling Applications', Goettinger Geographische Abhandlungen, Vol.115, p.13-27

#Hofierka, J., Suri, M. (2002): The solar radiation model for Open source GIS: implementation and applications. International GRASS users conference in Trento, Italy, September 2002

#Iwahashi, J. & Pike, R.J. (2007): Automated classifications of topography from DEMs by an unsupervised nested-means algorithm and a three-part geometric signature. Geomorphology, Vol. 86, pp. 409-440

#Gerlitz, L., Conrad, O., B?hner, J. (2015): Large scale atmospheric forcing and topographic modification of precipitation rates over High Asia - a neural network based approach. Earth System Dynamics, 6, 1-21. doi:10.5194/esd-6-1-2015.

#Rojas, O.E; IICA, San Jos? (Costa Rica). Proyecto de Agroclimatolog?a. Estudio agroclim?tico de Costa Rica.
#Material type: materialTypeLabelBook Call number: IICA-PM 617 Series: Publicaci?n Miscel?nea (IICA) no. 617. Publisher: San Jos?, C.R.: s.n., 1985Description: 183 p.ISSN: 0534-5391.

#-------------------------------------------------------------------------------------
# END OF CODE
#-------------------------------------------------------------------------------------
