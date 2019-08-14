# set working directory according to your computer
setwd("C:/Users/u0133977/Dropbox/Hypomirror/Chao_Plant/Data and Code")
setwd("C:/Users/gjbowen/Dropbox/Hypomirror/Chao_Plant/Data and Code")

# load packages
library(raster)
library(assignR)

# read precip isoscapes
precip_d2h_mean = raster("precip_NA_d2h/predkrig.tiff")
precip_d2h_sd = raster("precip_NA_d2h/stdkrig.tiff")
precip_d2h = stack(precip_d2h_mean, precip_d2h_sd)

precip_d18o_mean = raster("precip_NA_d18o/predkrig.tiff")
precip_d18o_sd = raster("precip_NA_d18o/stdkrig.tiff")
precip_d18o = stack(precip_d18o_mean, precip_d18o_sd)

# read plant isoscapes
plant_d18o_mean = raster("plant_NA_d18o/p_mean.tiff")
plant_d18o_hi = raster("plant_NA_d18o/p_hi.tiff")
plant_d18o_lo = raster("plant_NA_d18o/p_lo.tiff")

plant_d2h_mean = raster("plant_NA_d2h/p_mean.tiff")
plant_d2h_hi = raster("plant_NA_d2h/p_hi.tiff")
plant_d2h_lo = raster("plant_NA_d2h/p_lo.tiff")

# calc plant_d18o - use hi-lo range as 3 sigma range
plant_d18o_sd = (plant_d18o_hi - plant_d18o_lo) / 6
plant_d18o_sd = sqrt(plant_d18o_sd ^ 2 + precip_d18o_sd ^ 2)
plant_d18o = stack(plant_d18o_mean, plant_d18o_sd)
names(plant_d18o) = c("mean", "sd")

# calc plant_d2h - use hi-lo range as 3 sigma range
plant_d2h_sd = (plant_d2h_hi - plant_d2h_lo) / 6
plant_d2h_sd = sqrt(plant_d2h_sd ^ 2 + precip_d2h_sd ^ 2)
plant_d2h <- stack(plant_d2h_mean, plant_d2h_sd)
names(plant_d2h) = c("mean", "sd")

#Put everything on AEA projection, first isoscapes
aea = CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 
          +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
precip_d2h = projectRaster(precip_d2h, crs=aea)
precip_d18o = projectRaster(precip_d18o, crs=aea)
plant_d2h = projectRaster(plant_d2h, crs=aea)
plant_d18o = projectRaster(plant_d18o, crs=aea)

#Have a look
plot(plant_d2h)
plot(plant_d18o)
plot(precip_d2h)
plot(precip_d18o)

plot(plant_d2h$mean - precip_d2h$predkrig)
plot(plant_d2h$sd - precip_d2h$stdkrig)
plot(plant_d18o$mean - precip_d18o$predkrig)
plot(plant_d18o$sd - precip_d18o$stdkrig)

# read bird feather isotope data and make SPDF objects
bird = read.csv("Hobson2015.csv")
bird_d2h = SpatialPointsDataFrame(matrix(c(bird$Longitude, bird$Latitude),ncol=2), 
                                  data.frame("d2H" = bird$d2H), 
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
bird_d18o = SpatialPointsDataFrame(matrix(c(bird$Longitude, bird$Latitude),ncol=2), 
                                  data.frame("d18O" = bird$d18O), 
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#project vector files
bird_d2h = spTransform(bird_d2h, aea)
bird_d18o = spTransform(bird_d18o, aea)
naMap = spTransform(naMap, aea)

#### plot  isoscapes ####
png("Fig2.png", units = "cm", res = 600, width = 20, height = 25)
par(mfrow=c(4,2), mai=c(0,0,0,1))
plot(precip_d2h[[1]], axes=FALSE, box=FALSE, xlim=c(-2.5e6, 2.5e6), ylim=c(0,4e6), legend=FALSE)
plot(precip_d2h[[1]], legend.only=TRUE, smallplot=c(0.75,0.77,0.2,0.8), 
     legend.args=list(text=expression("Precipitation "*delta^{2}*"H"), side=4, line=3.5))
text(-2.35e6, 3.95e6, "(a)", cex=1.5)
points(bird_d18o)

plot(precip_d2h[[2]], axes=FALSE, box=FALSE, xlim=c(-2.5e6, 2.5e6), ylim=c(0,4e6), legend=FALSE)
plot(precip_d2h[[2]], legend.only=TRUE, smallplot=c(0.75,0.77,0.2,0.8), 
     legend.args=list(text=expression("Precipitation "*delta^{2}*"H 1"*sigma), side=4, line=3.5))
text(-2.35e6, 3.95e6, "(b)", cex=1.5)

plot(precip_d18o[[1]], axes=FALSE, box=FALSE, xlim=c(-2.5e6, 2.5e6), ylim=c(0,4e6), legend=FALSE)
plot(precip_d18o[[1]], legend.only=TRUE, smallplot=c(0.75,0.77,0.2,0.8), 
     legend.args=list(text=expression("Precipitation "*delta^{18}*"O"), side=4, line=3.5))
text(-2.35e6, 3.95e6, "(c)", cex=1.5)

plot(precip_d18o[[2]], axes=FALSE, box=FALSE, xlim=c(-2.5e6, 2.5e6), ylim=c(0,4e6), legend=FALSE)
plot(precip_d18o[[2]], legend.only=TRUE, smallplot=c(0.75,0.77,0.2,0.8), 
     legend.args=list(text=expression("Precipitation "*delta^{18}*"O 1"*sigma), side=4, line=3.5))
text(-2.35e6, 3.95e6, "(d)", cex=1.5)

plot(plant_d2h[[1]], axes=FALSE, box=FALSE, xlim=c(-2.5e6, 2.5e6), ylim=c(0,4e6), legend=FALSE)
plot(plant_d2h[[1]], legend.only=TRUE, smallplot=c(0.75,0.77,0.2,0.8), 
     legend.args=list(text=expression("Plant "*delta^{2}*"H"), side=4, line=3.5))
text(-2.35e6, 3.95e6, "(e)", cex=1.5)

plot(plant_d2h[[2]], axes=FALSE, box=FALSE, xlim=c(-2.5e6, 2.5e6), ylim=c(0,4e6), legend=FALSE)
plot(plant_d2h[[2]], legend.only=TRUE, smallplot=c(0.75,0.77,0.2,0.8), 
     legend.args=list(text=expression("Plant "*delta^{2}*"H 1"*sigma), side=4, line=3.5))
text(-2.35e6, 3.95e6, "(f)", cex=1.5)

plot(plant_d18o[[1]], axes=FALSE, box=FALSE, xlim=c(-2.5e6, 2.5e6), ylim=c(0,4e6), legend=FALSE)
plot(plant_d18o[[1]], legend.only=TRUE, smallplot=c(0.75,0.77,0.2,0.8), 
     legend.args=list(text=expression("Plant "*delta^{18}*"O"), side=4, line=3.5))
text(-2.35e6, 3.95e6, "(g)", cex=1.5)

plot(plant_d18o[[2]], axes=FALSE, box=FALSE, xlim=c(-2.5e6, 2.5e6), ylim=c(0,4e6), legend=FALSE)
plot(plant_d18o[[2]], legend.only=TRUE, smallplot=c(0.75,0.77,0.2,0.8), 
     legend.args=list(text=expression("Plant "*delta^{18}*"O 1"*sigma), side=4, line=3.5))
text(-2.35e6, 3.95e6, "(h)", cex=1.5)

dev.off()

#### do the calRasters on the full datasets ####
cr_precip_d2h = calRaster(bird_d2h, precip_d2h, sdMethod = 1, genplot = FALSE, savePDF = FALSE)
cr_plant_d2h = calRaster(bird_d2h, plant_d2h, sdMethod = 1, genplot = FALSE, savePDF = FALSE)
cr_precip_d18o = calRaster(bird_d18o, precip_d18o, sdMethod = 1, genplot = FALSE, savePDF = FALSE)
cr_plant_d18o = calRaster(bird_d18o, plant_d18o, sdMethod = 1, genplot = FALSE, savePDF = FALSE)

#### plot regression relationships for full datasets ####
png("../Fig3.png", units="in", res=600, width = 7, height = 6.5)
par(mar = c(5,5,0.5,0.5))
layout(matrix(c(1,2,3,4), nrow = 2))

plot(cr_precip_d2h$lm.data, pch=21, bg="white", 
     xlab=expression("Precipitation "*delta^{2}*"H"), 
     ylab=expression("Feather "*delta^{2}*"H"))
mod = cr_precip_d2h$lm.model
abline(mod)
xl = par("usr")[1] + (par("usr")[2] - par("usr")[1]) / 15
yl = par("usr")[4] - (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, "(a)")
eqn1 = paste0("y = ", round(mod$coefficients[2], 2), " * x ",
            ifelse(sign(mod$coefficients[1]) == 1, "+ ", "- "),
            abs(round(mod$coefficients[1],1)))
yl = par("usr")[4] - 3 * (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, eqn1, pos=4)
r2 = round(summary(mod)$adj.r.squared,2)
eqn2 = bquote(Adj ~ R^2 == .(r2))
yl = par("usr")[4] - 4 * (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, eqn2, pos=4)

plot(cr_precip_d18o$lm.data, pch=21, bg="white", 
     xlab=expression("Precipitation "*delta^{18}*"O"), 
     ylab=expression("Feather "*delta^{18}*"O"))
mod = cr_precip_d18o$lm.model
abline(mod)
xl = par("usr")[1] + (par("usr")[2] - par("usr")[1]) / 15
yl = par("usr")[4] - (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, "(c)")
eqn1 = paste0("y = ", round(mod$coefficients[2], 2), " * x ",
              ifelse(sign(mod$coefficients[1]) == 1, "+ ", "- "),
              abs(round(mod$coefficients[1],1)))
yl = par("usr")[4] - 3 * (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, eqn1, pos=4)
r2 = round(summary(mod)$adj.r.squared,2)
eqn2 = bquote(Adj ~ R^2 == .(r2))
yl = par("usr")[4] - 4 * (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, eqn2, pos=4)


plot(cr_plant_d2h$lm.data, pch=21, bg="white", 
     xlab=expression("Plant "*delta^{2}*"H"), 
     ylab=expression("Feather "*delta^{2}*"H"))
mod = cr_plant_d2h$lm.model
abline(mod)
xl = par("usr")[1] + (par("usr")[2] - par("usr")[1]) / 15
yl = par("usr")[4] - (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, "(b)")
eqn1 = paste0("y = ", round(mod$coefficients[2], 2), " * x ",
              ifelse(sign(mod$coefficients[1]) == 1, "+ ", "- "),
              abs(round(mod$coefficients[1],1)))
xl = par("usr")[1] + 2 * (par("usr")[2] - par("usr")[1]) / 15
yl = par("usr")[4] - 2 * (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, eqn1, pos=4)
r2 = round(summary(mod)$adj.r.squared,2)
eqn2 = bquote(Adj ~ R^2 == .(r2))
yl = par("usr")[4] - 3 * (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, eqn2, pos=4)

plot(cr_plant_d18o$lm.data, pch=21, bg="white", 
     xlab=expression("Plant "*delta^{18}*"O"), 
     ylab=expression("Feather "*delta^{18}*"O"))
mod = cr_plant_d18o$lm.model
abline(mod)
xl = par("usr")[1] + (par("usr")[2] - par("usr")[1]) / 15
yl = par("usr")[4] - (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, "(d)")
eqn1 = paste0("y = ", round(mod$coefficients[2], 2), " * x ",
              ifelse(sign(mod$coefficients[1]) == 1, "+ ", "- "),
              abs(round(mod$coefficients[1],1)))
xl = par("usr")[1] + 0.5 * (par("usr")[2] - par("usr")[1]) / 15
yl = par("usr")[4] - 2.5 * (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, eqn1, pos=4)
r2 = round(summary(mod)$adj.r.squared,2)
eqn2 = bquote(Adj ~ R^2 == .(r2))
yl = par("usr")[4] - 3.5 * (par("usr")[4] - par("usr")[3]) / 13
text(xl, yl, eqn2, pos=4)

dev.off()

#### calc QA ####
niter = 100
vali = 40
precip_d2h_QA <- QA(isoscape = precip_d2h, known = bird_d2h, valiStation = vali, 
                    valiTime = niter, setSeed = T)
precip_d18o_QA <- QA(precip_d18o, known = bird_d18o, valiStation = vali, 
                     valiTime = niter, setSeed = T)
plant_d2h_QA <- QA(plant_d2h, known = bird_d2h, valiStation = vali, 
                   valiTime = niter, setSeed = T)
plant_d18o_QA <- QA(plant_d18o, known = bird_d18o, valiStation = vali, 
                     valiTime = niter, setSeed = T)

# save output
save(precip_d2h_QA, file="precip_d2h_QA.RData")
save(precip_d18o_QA, file="precip_d18o_QA.RData")
save(plant_d2h_QA, file="plant_d2h_QA.RData")
save(plant_d18o_QA, file="plant_d18o_QA.RData")

# to continue from here...
load("precip_d2h_QA.RData")
load("precip_d18o_QA.RData")
load("plant_d2h_QA.RData")
load("plant_d18o_QA.RData")

#### plot all 4 models ####
qas = list("Precip H" = precip_d2h_QA, "Plant H" = plant_d2h_QA,
           "Precip O" = precip_d18o_QA, "Plant O" = plant_d18o_QA)
plot.QA(qas, savePDF = T)
