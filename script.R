#### Development of urban area in Da Nang and Quang Nam, Vietnam 
#### Change Detection between 1988 and 2020 
####  - a reproduction of Bachelor thesis by Katrin Wernicke - 
# download data here (until May 3rd 2021): 
# https://gigamove.rz.rwth-aachen.de/d/id/eWpSSM6GuT83fB	

## load packages 
library(raster)
library(sp)
library(RStoolbox)
library(rgdal)
library(sf)
library(ggplot2)
library(ggspatial)

### 1 Import of preprocessed and clipped Landsat Data ####
## dataset 2020: ImageCollection of Landsat 8 data, reduced in Google Earth Engine
## dataset 2010, 2000, 1988: ImageCollection of Landsat 5 data, reduced in Google Earth Engine

# import and reprojecting
imProj <- function(importData, year){
  dataset<- brick(paste0("data/", importData, ".tif"))
  dataset <- projectRaster(dataset, crs=crs("+init=epsg:3405"),res = 30)
  }

dataset2020 <- imProj("dataset2020", 2020)
dataset2010 <- imProj("dataset2010", 2010)
dataset2000 <- imProj("dataset2000", 2000)
dataset1988 <- imProj("dataset1988", 1988)

### 2 Vegetation Indices and Principal Component Analysis #####

## the following function calculates the following Indices, performs a Principal Component Analyis
## and stacks all raster layers into a layerstack: 
indicesPca<- function(dataset, bGREEN, bRED, bNIR, bSWIR, year){
  # NDVI Normalized Differenced Vegetation Index: (NIR-RED)/(NIR+RED)
  ndvi <- spectralIndices(dataset, red=bRED, nir=bNIR,indices="NDVI")
  # NDBI Normalized Differenced Built-Up Index (Zha et al. 2003): (SWIR-NIR)/(SWIR+NIR) 
  ndbi <- (dataset[[bSWIR]]-dataset[[bNIR]])/(dataset[[bSWIR]]+dataset[[bNIR]])
  # NDWI Normalized Differenced Water Index (McFeeters 1996): (GREEN-NIR)/(GREEN+NIR) 
  ndwi <- spectralIndices(dataset, green=bGREEN, nir=bNIR, indices="NDWI")
  # NDMI Normalized Differenced Moisture Index (Gao 1996)(also NDWI2):(NIR-SWIR)/(NIR+SWIR) 
  ndmi <- spectralIndices(dataset, red=bRED, nir=bNIR, swir2=bSWIR, indices="NDWI2")
  # Principal Component Analysis
  pca <- rasterPCA(dataset)
  # extract first principal component
  pc1 <- pca$map[[1]]
  ## create layerstack
  layerstack <- raster::stack(dataset, pc1, ndbi, ndvi, ndwi, ndmi)
  # save as raster file
  #writeRaster(x=layerstack, filename=paste0("results/layerstack", year, ".tif"), format="GTiff", overwrite=TRUE)
  return(layerstack)
}
# return layerstacks
layerstack2020 <- indicesPca(dataset2020,3,4,5,6,"2020")
layerstack2010 <- indicesPca(dataset2010,2,3,4,5,"2010")
layerstack2000 <- indicesPca(dataset2000,2,3,4,5,"2000")
layerstack1988 <- indicesPca(dataset1988,2,3,4,5,"1988")

### 3 Supervised Classification ####
classification <- function(layerstack, year){
  # define classes, for better contrast 2 vegetation classes
  keys <- c("water",  "soil", "veghell", "urban", "vegdunkel","sand")
  for (i in 1:length(keys)){
    # import training data (generated in SNAP, for each class one file)
    path <- paste0("data/aoi",year,"_", keys[i],"_Polygon.shp", collapse = NULL)
    current_poly <- readOGR(path)
    current_poly <- spTransform(current_poly, crs(layerstack))
    current_poly$class <- keys[i]
    current_poly$id <- i
    if (i == 1){
      my_poly <- current_poly
    } else {
      my_poly <- rbind(my_poly, current_poly)
    }
  }
  # classification, default is set to Random Forest
  sc <- superClass(layerstack, trainData = my_poly, responseCol ="id")
  # rename/reclassify both "veghell" and "vegdunkel" into class "vegetation" 
  sc_reclass <- reclassify(sc$map, cbind(5,3))
  # reclassify urban (2) to urban (4)
  sc_reclass <- reclassify(sc_reclass, cbind(6,2))
  # save classification results in a raster file
  #writeRaster(x=sc_reclass, filename=paste0("results/myClassification", year, ".tif"), format="GTiff", overwrite=TRUE)
  return(sc_reclass)
}

# return classification 
classified2020 <- classification(layerstack2020, "2020")
classified2010 <- classification(layerstack2010, "2010")
classified2000 <- classification(layerstack2000, "2000")
classified1988 <- classification(layerstack1988, "1988")

## BONUS: area calculation per classification and per class
area2020 <- tapply(area(classified2020), classified2020[], sum)
area2010 <- tapply(area(classified2010), classified2010[], sum)
area2000 <- tapply(area(classified2000), classified2000[], sum)
area1988 <- tapply(area(classified1988), classified1988[], sum)
Carea <- as.data.frame(rbind(area1988, area2000,area2010,area2020))

### 4 Accuracy Assessment ####
accuracy <- function(classified, year){
  # import validation points
  valP <- st_read(paste0("data/", year, "stratified_random_sampling_classified.gpkg"))
  # Accuracy Assessment
  val <- validateMap(classified, valP, "Classified", mode="classification")
  return(val)
}
# return accuracy per year
accuracy2020 <- accuracy(classified2020, "2020")
accuracy2010 <- accuracy(classified2010, "2010")
accuracy2000 <- accuracy(classified2000, "2000")
accuracy1988 <- accuracy(classified1988, "1988")

### 5 Change Detection
# extract urban class from classification
get_urban <- function(raster, i) {
  urban <- raster == 4
  return(reclassify(urban, cbind(1,i)))
}
# combine urban classes form all years in one raster to present changes
urbanChange<- mosaic(get_urban(classified2020,1),get_urban(classified2010,2),get_urban(classified2000,3),get_urban(classified1988,4), fun=max)

#save raster
#writeRaster(urbanChange, filename="results/urbanChange.tif",format="GTiff", overwrite=TRUE)

### 6 Plot results ####
ggplot()+
  ggR(urbanChange, geom_raster=T, ggLayer=T, alpha = 0.5)+
  coord_sf(crs=st_crs(urbanChange), datum=st_crs(3405))+
  scale_fill_gradientn(colors = c("white", "darkred", "red", "orange", "yellow"), 
                       name= "year",labels = c( "NA" , "2020", "2010", "2000", "1988"), 
                       guide = guide_colorbar(
                         barheight = unit(20, units = "mm"),
                         barwidth = unit(3, units = "mm"),
                         draw.ulim = F,
                         title.position = 'top',
                         title.hjust = 0.5,
                         label.hjust = 0.5
                       ))+
  theme(
    legend.text.align = 0,
    legend.text = element_text(size = 7, hjust = 0, color = "black"), #color = "#4e4d47"),
    plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "black"),  
    legend.title = element_text(size = 8, colour = "black"),
    plot.caption = element_text(size = 7, 
                                hjust = 0.5, 
                                color = "#939184")) +
  labs(x="Easting", y="Northing", 
      title = "Development of urban area 1988 - 2020", 
      subtitle = "in Da Nang and Quang Nam, Vietnam", 
      caption = "Data: Landsat 8, Landsat 5 ; Author: Katrin Wernicke;\nDate: 2021-04-19; Projection: UTM, WGS84")+
  annotation_scale()
