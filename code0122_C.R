
## load the library 
library(raster)
library(sp)
library(randomForest)
library(magrittr)
library(rgdal)
library(gstat)
library(ggplot2)
library(mlr)
library(SemiPar)
library(Hmisc)
library(foreign)
library(maptools)
library(prettymapr)
library(mlrMBO)
library(parallelMap)
library(caret)
library(automap)
library(reshape2)

## start the parallel 
parallelStartSocket(5)

WGS84 <- CRS("+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

study_area <- shapefile("~/WP2/data/study_area.shp")
water <- shapefile("~/WP2/data/water.shp")

## load the veg, soil, land use, groundwater subarea,
## surface water subare, catchment
Soil <- raster("~/WP2/data/soil1.ovr")
Veg <- raster("~/WP2/data/vegetation.ovr")
Land_use <- raster("~/WP2/data/landuse1.ovr")
Cat <- raster("~/WP2/data/catch_name2.tif.ovr")
DEM<-raster("~/WP2/data/topo_ProjectRaster3.tif.ovr")

## define the function 
## preprocess 
study_area <- spTransform(study_area, WGS84)
extent <- c(study_area@bbox[1, 1:2], study_area@bbox[2, 1:2])

water <- spTransform(water, WGS84)

pre <- function(x) {
  projection(x) <- WGS84
  extent(x) <- extent
  x <- raster::mask(x, study_area)
  return(x)
}

read_points <- function(read_data) {
  SP <- SpatialPoints(read_data[, 2:3], proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
  SP <- spTransform(SP, WGS84)
  SP@bbox <- study_area@bbox
  if (length(zerodist(SP)) >= 1) {
    SP <- SP[-(zerodist(SP)[, 1]),]
  }
  #   plot(study_area_withW)
  #  points(SP@coords)
  return(SP)
}

read_pointDataframes <- function(read_data) {
  SP <- SpatialPoints(read_data[, 2:3], proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
  SPD <- SpatialPointsDataFrame(SP, read_data)
  SPD <- spTransform(SPD, WGS84)
  SPD@bbox <- study_area@bbox
  if (length(zerodist(SPD)) >= 1) {
    SPD <- SPD[-(zerodist(SPD)[, 1]),]
  }
  # plot(study_area_withW)
  #points(SPD@coords)
  return(SPD)
}


reclass <- function(df, i, j) {
  df[, "DON"][df[, "DON"] <= i] <- "Low"
  df[, "DON"][df[, "DON"] < j] <- "Medium"
  df[, "DON"][(df[, "DON"] != "Low") & (df[, "DON"] != "Medium")] <- "High"
  df[, "DON"] <- factor(df[, "DON"], levels = c("Low", "Medium", "High"))
  return(df)
}


reclass4<-function(df,i,j){
  for (t in c(1,2)){
    df[, t][df[, t] <=i] <- "Low"
    df[, t][df[, t] < j] <- "Medium"
    df[, t][(df[, t] != "Low") & (df[, t] != "Medium")] <- "High"
    df[, t] <- factor(df[, t], levels = c("Low", "Medium", "High"))
  }
  return(df)
}

# Add X and Y to training 
add_S1S2 <- function(dataset) {
  dataset$s1 <- coordinates(dataset)[, 1]
  dataset$s2 <- coordinates(dataset)[, 2]
  return(dataset)
}

get_landscape<-function(df){
  landscape_all<-data.frame()
  for (ii in seq(1,length(df))){
    aa<-as.data.frame(df[[ii]])
    aa<-subset(aa,aa$Soil!="NA")
    Soil=tail(names(sort(table(aa[,1]))),1)
    Veg=tail(names(sort(table(aa[,2]))),1)
    Landuse=tail(names(sort(table(aa[,3]))),1)
    Catchment=tail(names(sort(table(aa[,4]))),1)
    GW_depth=mean(aa[,5])
    Distance=mean(aa[,6])
    Distance_GWC=mean(aa[,7])
    sing_land<-data.frame(Soil,Veg,Landuse,Catchment,GW_depth,Distance,Distance_GWC)
    landscape_all<-rbind(landscape_all,sing_land)
  }
  return(landscape_all)
}

## preprocess the landscape raster
Soil <- pre(Soil)
Veg <- pre(Veg)
Land_use <- pre(Land_use)
Cat <- pre(Cat)

v_Veg<-values(Veg)
v_Veg[v_Veg %in% c(2,3,4)]=1
v_Veg[v_Veg %in% c(8,9)]=8
v_Veg[v_Veg %in% c(12,13)]=12
v_Veg[v_Veg %in% c(18,19,20)]=18
values(Veg)<-v_Veg

v_land<-values(Land_use)
v_land[v_land %in% c(1,2,5,6,7,11,12,13)]=1
v_land[v_land %in% c(3,4)]=3
v_land[v_land %in% c(8,10)]=8
values(Land_use)<-v_land

#v_soil<-values(Soil)
#v_soil[v_soil %in% c(1,2)]=1
#v_soil[v_soil %in% c(4,5)]=4
#v_soil[v_soil %in% c(6,7)]=6
#v_soil[v_soil %in% c(11,12)]=11
#v_soil[v_soil %in% c(13,14)]=13
#values(Soil)<-v_soil

# Create an empty grid where n is the total number of cells
r <- raster(study_area)
res(r) <- res(Soil) # 10 km if your CRS's units are in km
base_grid <- as(r, 'SpatialGrid')
#plot(base_grid)

## M2, using RF to predict the DON
depth <- read.csv("~/WP2/data/sampling_depth.csv",header=T) %>% read_pointDataframes(.)

# Define the 1st order polynomial equation
f_depth <- as.formula(sampling_d ~ 1)
# Add X and Y to training 
depth<-add_S1S2(depth)
# variogram on the de-trended data.
var.depth <- variogram(f_depth, depth)
#plot(var.depth)
dat.fit_depth <- fit.variogram(var.depth,vgm(c("Sph","Exp")))
# created in the earlier step)
depth_k <- krige(f_depth, depth, base_grid, dat.fit_depth) %>% raster(.) %>% raster::mask(., study_area)
#plot(depth_k)
depth_k@data@names<-"GW_depth"

#Now make the map
### distance
water <- raster::rasterize(water, depth_k)
water_distance <- raster::mask(distance(water),study_area)
water_distance@data@names<-"Distance_to_water"

GW_center<-data.frame(Latitude=c(6495000,6475000,6460000,6448000,6403000),Longitude=rep(402000,5),values=1)
GW_center <- SpatialPoints(GW_center[, c(2:1)], proj4string = WGS84)
GW_center@bbox <- study_area@bbox
base_GWC<-water 
values(base_GWC)<-1
Distance_GWC<-distanceFromPoints(base_GWC,GW_center)
Distance_GWC@data@names<-"Distance_GWC"

## load the data 
landscapes<-stack(Soil,Veg,Land_use,Cat,depth_k,water_distance,Distance_GWC)
names(landscapes) <- c("Soil", "Veg", "Landuse","Catchment", "GW_depth", "Distance","Distance_GWC")


## load the data 
set.seed(666)

all_results<-data.frame()
all_data<-read.csv("~/WP2/data/all_data1127.csv",header = T)

results<-data.frame()

a1=0.5
a2=2.5

set.seed(91)
trainIndex <- createDataPartition(all_data$DON, p = .75, list = FALSE, times = 1)

training<-all_data[ trainIndex,]
testing<-all_data[-trainIndex,]

#Make a distance matrix
## cross validation
## 10-fold cross-validation
rdesc = makeResampleDesc("CV", iters = 10,stratify = TRUE)
## Classification tree, set it up for predicting probabilities
reg_rf = makeLearner("regr.randomForest")

## set the model parameters for random forest
para_rf = makeParamSet(
  makeIntegerParam("ntree", lower = 400, upper = 1000),
  makeIntegerParam("mtry", lower = 2, upper = 7),
  makeIntegerParam("nodesize", lower = 1, upper = 6)
)

## define the search stratgy
ctrl = makeTuneControlIrace(maxExperiments = 200L)

model_build <- function(dataset, n_target) {
  #set.seed(719)
  ## define the regression task for DON 
  WP3_target = makeRegrTask(id = "WP3_target", data = dataset, target = n_target)
  rin = makeResampleInstance(rdesc, task = WP3_target)
  res_rf = mlr::tuneParams(reg_rf, WP3_target, resampling = rdesc, par.set = para_rf, control = ctrl,
                           show.info = FALSE)
  lrn_rf = setHyperPars(reg_rf, par.vals = res_rf$x)
  ## train the final model 
  #set.seed(719)
  rf <- mlr::train(lrn_rf, WP3_target)
  return(rf)
}
  ## load the point data 
  training_df <- read_pointDataframes(training)
  testing_df <-  read_pointDataframes(testing) 
  
  training_points<- read_points(training)
  testing_points <- read_points(testing)
  
  ## map1, using kringing for DON interpolation
  # Add X and Y to training 
  training_df<-add_S1S2(training_df)
  testing_df<-add_S1S2(testing_df)
  
  ## M2, using RF to predict the DON
  a=100
  b=200
  capture_zone_land<-function(df){
    num<-nrow(df)
    landscape_data<-data.frame()
    for (r in seq(1,num)){
      p1_long<-df@coords[r,1]
      p1_lat<-df@coords[r,2]
      pg<-spPolygons(rbind(c(p1_long,p1_lat),c(p1_long+a,p1_lat+b),c(p1_long+2*a,p1_lat+b),
                           c(p1_long+2*a,p1_lat-b),c(p1_long+a,p1_lat-b),c(p1_long,p1_lat)))  
      projection(pg)<- WGS84
      p1_landscape<-raster::extract(landscapes,pg)
      p1_landscape<-get_landscape(p1_landscape)
      landscape_data<-rbind(landscape_data,p1_landscape)
    }
    return(landscape_data)
  }
  
  landscape_train <- capture_zone_land(training_df)
  landscape_test <- capture_zone_land(testing_df)
  
  M2_train <- cbind(as.data.frame(landscape_train), training_df@data[c("DON","date_","Collect_Month")])
  M2_test <- cbind(as.data.frame(landscape_test), testing_df@data[c("DON","date_","Collect_Month")])
  
  names(M2_train) <- colnames(M2_test)
  
  common_landscape<-function(land){
    land_dataset<-data.frame(table(M2_train[,land]))
    land_common<-subset(land_dataset,land_dataset[,2]==max(land_dataset[,2]))[1]
    return(as.matrix(land_common))
  }
  
  soil_max = common_landscape("Soil")[1]
  veg_max=common_landscape("Veg")[1]
  landuse_max = common_landscape("Landuse")[1]
  cat_max = common_landscape("Catchment")[1]
  
  max_list<-list(soil_max,veg_max,landuse_max,cat_max)
  
  for (ii in seq(1,4)){
    M2_train[,ii]<-factor(M2_train[,ii],levels = unique(values(landscapes[[ii]]))[-1])
    M2_test[,ii]<-factor(M2_test[,ii],levels=unique(values(landscapes[[ii]]))[-1])
    M2_test [(which(!(M2_test[,ii] %in% M2_train[,ii]))),ii]<-as.numeric(max_list[[ii]])
     
    M2_train[,ii]<-droplevels(M2_train[,ii])
    M2_test[,ii]<-factor(M2_test[,ii],levels = levels(M2_train[,ii]))
    }
  
   #M2_train<-reclass(M2_train,a1,a2)
   #M2_test<-reclass(M2_test,a1,a2)
   
#   M2_train$Collect_Month<-factor(M2_train$Collect_Month,levels=c("1","2","3","4","5","6","7","8","9","10","11"))
 #  M2_test$Collect_Month<-factor(M2_test$Collect_Month,levels=c("1","2","3","4","5","6","7","8","9","10","11"))
   
   
#  M2_train$DON<-log10(M2_train$DON)
#  M2_test$DON<-log10(M2_test$DON)
  
  #for(i in c("GW_depth","Distance","date_","DOC")) {
    
  #  min_train<-min(M2_train[,i])
  #  max_train<-max(M2_train[,i])
    
 #   M2_train[,i]<-(M2_train[,i]-min_train)/(max_train-min_train)
 #   M2_test[,i]<-(M2_test[,i]-min_train)/(max_train-min_train)

#  }
  
 # for(i in c("Distance_GWC","Collect_Month")){
    
 #   min_train<-min(M2_train[,i])
 #   max_train<-max(M2_train[,i])
    
 #   M2_train[,i]<-(max_train-M2_train[,i])/(max_train-min_train)
  #  M2_test[,i]<-(max_train-M2_test[,i])/(max_train-min_train)

#  }

  WP2Train<-M2_train[,-c(4)]
  WP2Test<-M2_test[,-c(4)]
  
  rf_DON_m2 <- model_build(WP2Train,"DON")

## test in testing set
test_rf = predict(rf_DON_m2, newdata = WP2Test)
## ConfusionMatrix
#print(calculateConfusionMatrix(test_rf))
## get the prediction performance
train_rf = predict(rf_DON_m2, newdata = WP2Train)

print(postResample(test_rf$data$response,test_rf$data$truth))
print(postResample(train_rf$data$response,train_rf$data$truth))







  