################################################################################
#
#                   St-Lawrence sedimentation analysis
#
#                     Part 3: Mixed-effects model
#
################################################################################

library(nlme)
library(splines)
library(rgdal)
library(raster)
library(maptools)
library(fields)
library(colorspace)

source("0_Functions.R")

#-----------------------
#  Data loading
#-----------------------

#----- Sedimenation data
load("Data/ProcessedSediment.RData")

#----- Create data
n <- nrow(sedimentSeries)
nt <- ncol(sedimentSeries)

# Create data.frame
regdat <- data.frame(sed = as.vector(t(sedimentSeries)), 
  ddif = scale(as.vector(t(dateDiff)), scale = F), 
  year = scale(rep(1:nt, n), scale = F),
  loc = rep(as.numeric(rownames(sedimentSeries)), each = nt)
)

# As a first approximation for those to remove
# regdat$ddif[regdat$ddif > 345] <- NA

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
#  Model
#-----------------------
  
#----- Apply mixed effect model
mixedres <- lme(sed ~ year + ddif, random = ~ year | loc, 
  data = regdat, na.action = na.omit, control = list(opt = "optim"))

# Summary of fitted model
summary(mixedres) 

# Confidence intervals of estimates
intervals(mixedres)

#----- Plot results on a map
x11(width = 20, height = 10)
par(mar = c(5,4,4,7) + .1, mfrow = c(1,2))

#Intercept
plot.pts(locations, borders = stLawrence, values = coef(mixedres)[,1], 
  increment = 500, pal = c("blue", "grey", "red"), 
  args.plot = list(axes = F, pch = 16, ylab = "", xlab = "", cex = .7), 
  args.colorbar = list(zlim = fixef(mixedres)[1] + max(abs(ranef(mixedres)[,1])) * c(-1, 1), 
    type.labels = "continuous", 
    zlab = expression(paste("Mean ", alpha + a[i],  " (m)")), 
    formatC.pars = list(digits = 2, format = "f"), axis.pars = list(cex = 1.2), 
    label.pars = list(cex = 1.3, pos = 3)), 
  args.secteurs = list(pos = 4, cex = 0.5), args.limits = list(xpd = T, 
    labels = c('Montréal', 'Trois-Rivières'), cex= 1.5, adj = c(1,-1)), 
  args.borders = list(border="darkgrey")
)
text(par("usr")[1], par("usr")[4], "A", cex = 3, xpd = T)

# Trend
plot.pts(locations, borders = stLawrence, values = coef(mixedres)[,"year"], 
  increment = 500, pal = diverge_hcl(12, h = c(246, 40), c = 96), 
  args.plot = list(axes = F, pch = 16, ylab = "", xlab = "", cex = .7), 
  args.colorbar = list(zlim = max(abs(coef(mixedres)[,"year"])) * c(-1, 1), 
    type.labels = "continuous", 
    zlab = expression(paste("Trend ", beta + b[i],  " (m / year)")), 
    formatC.pars = list(digits = 3, format = "f"), axis.pars = list(cex = 1), 
    label.pars = list(cex = 1.3, pos = 3)), 
  args.secteurs = list(pos = 4, cex = 0.5), args.limits = list(xpd = T, 
    labels = c('Montréal', 'Trois-Rivières'), cex= 1.5, adj = c(1,-1)), 
  args.borders = list(border="darkgrey")
)
text(par("usr")[1], par("usr")[4], "B", cex = 3, xpd = T)

dev.print(png, filename = "Figures/Figure2.png", res = 200, units = "in")
setEPS()
dev.print(postscript, file = "Figures/Figure2.eps")

#----- Model assessment
norandomtrend <- update(mixedres, random = ~ 1 | loc)
norandomeff <- lm(sed ~ year + ddif, data = regdat, na.action = na.omit)

comparison <- anova(mixedres, norandomtrend, norandomeff)

write.table(comparison, file = "Figures/Table1.csv", sep = ",")

#----- Diagnostic
x11(width = 10, height = 10)
layout(matrix(c(1,1,2,2,0,3,3,0), nrow = 2, byrow = T))

# Residual plot vs year
plot(regdat$year[-mixedres$na.action], residuals(mixedres),
  ylab = "Residuals", xlab = "year", xaxt = "n")
lines(lowess(residuals(mixedres) ~ regdat$year[-mixedres$na.action]),
  col = "red", lwd = 2)
axis(1, at = pretty(1:nt) - mean(1:nt), labels = pretty(1:nt) + 2002)
text(par("usr")[1], par("usr")[4], "A", cex = 3, xpd = T, adj = c(1.5, -.5))

# Residual plot vs date difference
plot(regdat$ddif[-mixedres$na.action], residuals(mixedres),
  ylab = "Residuals", xlab = "Date difference", xaxt = "n")
lines(lowess(residuals(mixedres) ~ regdat$ddif[-mixedres$na.action]),
  col = "red", lwd = 2)
axis(1, at = pretty(dateDiff) - mean(dateDiff, na.rm = T), 
    labels = pretty(dateDiff))
text(par("usr")[1], par("usr")[4], "B", cex = 3, xpd = T, adj = c(1.5, -.5))
  
# Distribution of residuals
hist(residuals(mixedres), breaks = 20, xlab = "Residuals", main = "", col = 4)
text(par("usr")[1], par("usr")[4], "C", cex = 3, xpd = T, adj = c(1.5, -.5))

dev.print(png, filename = "Figures/SupFigure2.png", res = 200, units = "in")
