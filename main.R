## Louise Searle
## January 13 2015
## Lesson 7 Exercise.

## Libraries.
library(raster)
library(maptools)
## Load data.
load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/vcfGewata.rda")
load("data/trainingPoly.rda")

## Create raster brick of bands.
gewata <- brick(GewataB1,GewataB2,GewataB3,GewataB4,GewataB5,GewataB7)
## View band pairs scatterplots.
pairs(gewata)

## Replace unwanted values with NA.
vcfGewata[vcfGewata > 100] <- NA
## Add VCF layer to raster brick.
covariates <- addLayer(gewata, vcfGewata)
## Correct names.
names(covariates) <- c("band1","band2","band3","band4","band5","band7","VCF")

## Extract covariates values to matrix.
covValues <- getValues(covariates)
## Remove NA values from matrix.
covValues <- na.omit(covValues)
## Convert matrix to a data frame.
covariatesDF <- as.data.frame(covValues)
## Check data frame.
head(covariatesDF, n = 10)

## Create linear model.
modelVCF <- lm(formula= VCF ~ band1 + band2 + band3 + band4 + band5 + band7, data=covariatesDF)
summary(modelVCF)

## Predict raster values.
predVCF <- predict(covariates, model=modelVCF, na.rm=TRUE)

## Calculate actual difference.
diffVCF <- vcfGewata - predVCF
## Remove out of range values for predicted values for ease of plotting.
copyVCF <- predVCF
copyVCF[copyVCF > 100] <- NA
copyVCF[copyVCF < 0] <- NA
copyDiff <- vcfGewata - copyVCF
## Plot actual, predicted and difference VCF.
opar <- par(mfrow=c(1,1))
plot(copyVCF, main='VCF Predicted')
plot(vcfGewata, main ='VCF Actual')
plot(copyDiff, main='VCF Difference', col=bpy.colors(10))

## Square difference values per cell.
residVCF <- diffVCF^2
## For all cells, find the square root of the mean value.
rmseVCF <- sqrt(cellStats(residVCF, mean, na.rm=T))
## Check RMSE.
rmseVCF

## Rasterise landuse polygons, with class values.
classes <- rasterize(trainingPoly, diffVCF, field='Class')
## Calculate mean difference per zone.
diffClass <- zonal(diffVCF, classes, fun=mean, na.rm=T)
## Convert to data frame.
diffClass <- data.frame(diffClass)
## Add landuse names to data frame, according to zone code.
diffClass$class <- levels(trainingPoly$Class)
## Assign correct names.
names(diffClass) <- c('Zone', 'Difference', 'Class')
## Check difference.
diffClass

## Calculate mean of square of differences per cell, per class.
meanClass <- zonal(residVCF, classes, fun=mean, na.rm=T)
## Convert to data frame.
rmseClass <- data.frame(meanClass)
## Calculate square root of mean.
rmseClass$value <- sqrt(rmseClass$value)
## Add landuse names to data frame, according to zone code.
rmseClass$class <- levels(trainingPoly$Class)
## Assign correct names.
names(rmseClass) <- c('Zone', 'RMSE', 'Class')
## Check table.
rmseClass
