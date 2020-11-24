################################################################################
#
#                   St-Lawrence sedimentation analysis
#
#                     Part 1: Data preparation
#
################################################################################

library(sp)
library(rgdal)
library(raster)
library(maps)
library(rgeos)
library(spatstat)

source("0_Functions.R")

#------------------------
# Channel data reading
#------------------------

sectors <- readOGR(dsn = "Data/Limites du chenal/Secteurs de sondage.shp")
sectors <- spTransform(sectors, CRS("+init=epsg:2950"))

chainage <- readOGR(dsn = "Data/Limites du chenal/Points de chainage.shp")
chainage <- spTransform(chainage, CRS("+init=epsg:2950"))

centre <- readOGR(dsn = "Data/Limites du chenal/Ligne de Centre du Chenal.shp")
centre <- spTransform(centre,CRS("+init=epsg:2950"))

land_borders <- readOGR(dsn = "Data/10m_physical/ne_10m_land.shp")

#----- Select only sectors upriver to Trois-Rivières -----
# Trois-Rivières coordinates (EPSG 2950)
TRcoords <- c(376377, 5134896)
# We select only sectors to the west of Trois-Rivières
mask.extent <- extent(c(extent(sectors)[1], TRcoords[1], extent(sectors)[3:4]))
selected.sectors <- intersect(sectors, mask.extent) 


#------------------------
# Sedimentation data reading by sector /!\ Very long
#------------------------
sects <- unique(substr(selected.sectors$name,1,3))
ns <- length(sects)
# Objects to store sedimentation raster and sounding dates 
sedimentRaster <- vector("list",ns); names(sedimentRaster) <- sects
sedimentDates <- vector("list",ns); names(sedimentDates) <- sects
# Loop on sectors
for (i in 1:ns){
  # Read bathymetry data and compute difference to obtain sedimentation raster 
  result <- try(compute_sedimentation(sects[i], path = "E:/Bathymetry", res = 5, 
    fun = c("max","min"), pb_center = "center", years = 3:18,
    polygons = sectors[grep(sects[i], sectors$name),])
  )
  if (!inherits(result, "try-error")){
    sedimentRaster[[sects[i]]] <- result$raster
    origin(sedimentRaster[[sects[i]]]) <- 0 # Useful for rasters merging
    sedimentDates[[sects[i]]] <- result$dates
  }
}

#### We save the sedimentation rasters
# save(sedimentRaster, sedimentDates, file="Data/SedimentationRasters.Rdata")

#------------------------
# Merging all sectors sedimentation rasters together
#------------------------

# load("Data/SedimentationRasters.Rdata")  

# Merging the rasters
mergedRaster <- do.call(merge, unname(sedimentRaster))
names(mergedRaster) <- sprintf("A%s",formatC(2:18,width=2,flag="0"))


#------------------------
# Aggregating the merged raster into time series representing zones
#------------------------
increment <- 500 # The length of zones (in metres)

# Function that aggregates the raster
final.data <- aggregate_sedimentation(mergedRaster, increment, centre) #/!\ Long

# Remove the zone column and first column for which we didn't have autumn soundings
sedimentation.matrix <- matrix(NA, length(final.data$pt), 
    ncol(final.data$series) - 2)
sedimentation.matrix[final.data$series[,1],] <- final.data$series[,-(1:2)]

# Reorder the matrix to follow the spatial ordering (starting at Montreal)
orderFromOrigin <- order(rowSums(final.data$pt@coords))
final.pt <- final.data$pt[orderFromOrigin]
sedimentation.matrix <- sedimentation.matrix[orderFromOrigin,]
rownames(sedimentation.matrix) <- seq_along(final.data$pt)

# Remove zones with 5 or less years of data
moreThanFive <- apply(sedimentation.matrix, 1, function(x) sum(!is.na(x))>5)
sedimentSeries <- sedimentation.matrix[moreThanFive,]
locations <- final.pt[rownames(sedimentSeries)]

# Name years
colnames(sedimentSeries) <- sprintf("20%s",
    substr(colnames(final.data$series[,-(1:2)]), 2, 3))

# Sample size
nSeries <- nrow(sedimentSeries)
nYear <- ncol(sedimentSeries)

#------------------------
# Compute date difference for each value
#------------------------
dateDiff <- matrix(NA, nSeries, nYear)
dimnames(dateDiff) <- dimnames(sedimentSeries)
for (i in seq_len(nSeries)){
    secti <- intersect(selected.sectors, locations[i])
    sectNam <- substr(secti$name,1,3)
    dati <- sedimentDates[[sectNam]]
    if (!is.null(dati)) dateDiff[i,] <- dati$Diff[-1]
}
dateDiff[1,] <- dateDiff[2,]

#----- SAVE FINAL DATA -----
save(locations, sedimentSeries, dateDiff,  
  file = "Data/ProcessedSediment.RData")