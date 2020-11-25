################################################################################
#
#                   St-Lawrence sedimentation analysis
#
#                     Part 2: Descriptive statistics
#
################################################################################

library(rgdal)
library(raster)
library(maptools)
library(fields)

source("0_Functions.R")

#-----------------------
#  Data loading
#-----------------------

#----- Siltation layer thickness data
load("Data/ProcessedSediment.RData")

n <- nrow(sedimentSeries)
nt <- ncol(sedimentSeries)

#----- St-Lawrence borders
land_borders <- readOGR(dsn = "Data/10m_physical/ne_10m_land.shp")
SLlon <- c(-74,-72)
SLlat <- c(45,46.5)
stLawrence <- intersect(land_borders, extent(c(SLlon, SLlat)))
stLawrence <- spTransform(stLawrence, CRS("+init=epsg:2950"))
yshift <- -1500
xshift <- 0
stLawrence <- elide(stLawrence, shift = c(xshift, yshift))

#-----------------------
#  Number of available data per location (supplementary materials)
#-----------------------

nObsLoc <- apply(sedimentSeries, 1, function(x) sum(!is.na(x)))

x11(width = 15)
par(mfrow = c(1,2), oma = c(0, 0, 0, 5))

barplot(100 * table(nObsLoc) / n, xlab = "Number of data years", col = 4, 
    ylab = "Proportion (%)")
text(par("usr")[1], par("usr")[4], "A", cex = 3, xpd = T, pos = 3)

plot.pts(locations, borders = stLawrence, values = nObsLoc, 
  increment = 500, pal = tim.colors, nc = 10, 
  args.plot = list(axes = F, pch = 16, ylab = "", xlab = "", cex = .7), 
  args.colorbar = list(zlim = c(6, 15), type.labels = "classes", 
    zlab = "Number of data years", axis.pars = list(cex = 1), 
    label.pars = list(cex = 1.3, pos = 3)), 
  args.secteurs = list(pos = 4, cex = 0.5), args.limits = list(xpd = T, 
    labels = c('Montréal', 'Trois-Rivières'), cex= 1.5, adj = c(1,-1)), 
  args.borders = list(border="darkgrey")
)
text(par("usr")[1], par("usr")[4], "B", cex = 3, xpd = T, pos = 3)

dev.print(png, filename = "Figures/SupFigure1.png", res = 200, units = "in")

#-----------------------
#  Summary statistics (section 3.1)
#-----------------------

# Basic summary
summary(c(sedimentSeries))
sd(sedimentSeries, na.rm = T)

# Zone-specific maxima
zonemax <- apply(sedimentSeries, 1, max, na.rm = T)
range(zonemax)

# Where is the grand maximum
grandmax <- which(sedimentSeries == max(sedimentSeries, na.rm = T), arr.ind = T)
plot(stLawrence)
points(locations[grandmax[1]], pch = 16)

# Where is the minimum
grandmin <- which(sedimentSeries == min(sedimentSeries, na.rm = T), arr.ind = T)
plot(stLawrence)
points(locations[grandmin[1]], pch = 16)

# Number of negative values
sum(sedimentSeries < 0, na.rm = T); mean(sedimentSeries < 0, na.rm = T)

# Number of maxima above .30m
sum(zonemax > 0.3, na.rm = T); mean(zonemax > 0.3, na.rm = T)
