# https://r-pkgs.org/git.html


#### Development of urban extent in Da Nang and Quang Nam, Vietnam 
#### Change Detection between 1988 and 2020 

### 1 Import of preprocessed and clipped Landsat Data ####
# preprocessing could also be done in here with GEE package
# import raster composit from 2020

imProj <- function(importData, year){
  dataset<- brick(paste0("data/", importData, ".tif"))
  dataset <- projectRaster(dataset, crs=crs("+init=epsg:3405"),res = 30)
  }

dataset2020 <- imProj("dataset2020", 2020)
dataset2010 <- imProj("dataset2010", 2010)
dataset2000 <- imProj("dataset2000", 2000)
dataset1988 <- imProj("dataset1988", 1988)

### 2 Vegetation Indices and Principal Component Analysis #####

# 
## the following function calculates the following Indices, performs a Principal Component Analyis
## and stacks all raster layers into a layerstack: 
# NDVI Normalized Differenced Vegetation Index
# NDBI Normalized Differenced Built-Up Index
# NDWI Normalized Differenced Water Index
# NDMI Normalized Differenced Moisture Index

indicesPca<- function(dataset, bGREEN, bRED, bNIR, bSWIR, year){
  # NDVI (NIR-RED)/(NIR+RED)
  ndvi <- spectralIndices(dataset, red=bRED, nir=bNIR,indices="NDVI")
  # NDBI (SWIR-NIR)/(SWIR+NIR) Zha et al. 2003
  ndbi <- (dataset[[bSWIR]]-dataset[[bNIR]])/(dataset[[bSWIR]]+dataset[[bNIR]])
  # NDWI (GREEN-NIR)/(GREEN+NIR) McFeeters 1996
  ndwi <- spectralIndices(dataset, green=bGREEN, nir=bNIR, indices="NDWI")
  # NDMI (NDWI2) (NIR-SWIR)/(NIR+SWIR) Gao 1996
  ndmi <- spectralIndices(dataset, red=bRED, nir=bNIR, swir2=bSWIR, indices="NDWI2")
  # Principal Component Analysis
  pca <- rasterPCA(dataset)
  # extract first principal component
  pc1 <- pca$map[[1]]
  ## create layerstack
  layerstack <- raster::stack(dataset, pc1, ndbi, ndvi, ndwi, ndmi)
  # save as raster file
  writeRaster(x=layerstack, filename=paste0("results/layerstack", year, ".tif"), format="GTiff", overwrite=TRUE)
  return(layerstack)
}

## return layerstacks
layerstack2020 <- indicesPca(dataset2020,3,4,5,6,"2020")
layerstack2010 <- indicesPca(dataset2010,2,3,4,5,"2010")
layerstack2000 <- indicesPca(dataset2000,2,3,4,5,"2000")
layerstack1988 <- indicesPca(dataset1988,2,3,4,5,"1988")

### 3 Supervised Classification ####

classification <- function(layerstack, year){
  # define classes, for better contrast 2 vegetation classes
  keys <- c("water",  "soil", "veghell", "urban", "vegdunkel","sand")
  for (i in 1:length(keys)){
    # import training data (generated in SNAP, for each class one shp)
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
  #writeOGR(my_poly, dsn= "results/", layer = paste0("training_data", year, ".shp", collapse=NULL), driver ="ESRI Shapefile")
  # classification, default is set to Random Forest
  sc <- superClass(layerstack, trainData = my_poly, responseCol ="id")
  # rename/reclassify both veghell und vegdunkel zu einer Klasse
  sc_reclass <- reclassify(sc$map, cbind(5,3))
  # reclassify urban (2) to urban (4)
  sc_reclass <- reclassify(sc_reclass, cbind(6,2))
  writeRaster(x=sc_reclass, filename=paste0("results/myClassification", year, ".tif"), format="GTiff", overwrite=TRUE)
  return(sc_reclass)
}

# return classification 
classified2020 <- classification(layerstack2020, "2020")
classified2010 <- classification(layerstack2010, "2010")
classified2000 <- classification(layerstack2000, "2000")
classified1988 <- classification(layerstack1988, "1988")
### Accuracy Assessment####

accuracy <- function(classified, year){
  # import validation points
  valP <- st_read(paste0("data/", year, "stratified_random_sampling_classified.gpkg"))
  # Accuracy Assessment
  val <- validateMap(classified, valP, "Classified", mode="classification")
  return(val)
}
accuracy2020 <- accuracy(classified2020, "2020")
accuracy2010 <- accuracy(classified2010, "2010")
accuracy2000 <- accuracy(classified2000, "2000")
accuracy1988 <- accuracy(classified1988, "1988")

## change detection
# extract urban class from classification

get_urban <- function(raster, i) {
  urban <- raster == 4
  return(reclassify(urban, cbind(1,i)))
}

urbanChange<- mosaic(get_urban(classified2020,1),get_urban(classified2010,2),get_urban(classified2000,3),get_urban(classified1988,4), fun=max)
writeRaster(urbanChange, filename="results/urbanChange.tif",format="GTiff", overwrite=TRUE)
## plot results

plot(mm)

ggplot()+
  ggR(urbanChange, geom_raster=T, ggLayer=T)+
  scale_fill_manual(values=c("yellow", "orange", "red", "brown"), name= "Land Cover")+
  coord_sf(crs=st_crs(urbanChange), datum=st_crs(3405))

