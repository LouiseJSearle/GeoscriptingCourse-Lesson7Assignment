## Louise Searle
## January 13 2015
## Lesson 7 Exercise.

library(raster)
library(maptools)

## load data
load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/vcfGewata.rda")
load("data/trainingPoly.rda")

gewata <- brick(GewataB1,GewataB2,GewataB3,GewataB4,GewataB5,GewataB7)
pairs(gewata)

vcfGewata[vcfGewata > 100] <- NA
covariates <- addLayer(gewata, vcfGewata)
names(covariates) <- c("band1","band2","band3","band4","band5","band7","VCF")

covValues <- getValues(covariates)
covValues <- na.omit(covValues)
covariatesDF <- as.data.frame(covValues)

# gewataDF <- data.frame(band1=GewataB1@data@values, 
#                        band2=GewataB2@data@values, 
#                        band3=GewataB3@data@values, 
#                        band4=GewataB4@data@values, 
#                        band5=GewataB5@data@values, 
#                        band7=GewataB7@data@values,
#                        VCF=vcfGewata@data@values)

modelVCF <- lm(formula= VCF ~ band1 + band2 + band3 + band4 + band5 + band7, data=covariatesDF)

summary(modelVCF)
# Call:
#      lm(formula = VCF ~ band1 + band2 + band3 + band4 + band5 + band7, 
#         data = gewataDF)
# 
# Residuals:
#      Min      1Q  Median      3Q     Max 
# -58.205  -4.636   0.715   5.211 258.436 
# 
# Coefficients:
#                   Estimate   Std. Error t value   Pr(>|t|)    
# (Intercept)       8.549e+01  5.723e-02 1493.805   <2e-16 ***
#      band1        9.291e-02  1.889e-04  491.817   <2e-16 ***
#      band2       -1.505e-01  2.353e-04 -639.622   <2e-16 ***
#      band3       -2.528e-03  1.698e-04  -14.889   <2e-16 ***
#      band4        1.661e-02  2.790e-05  595.378   <2e-16 ***
#      band5       -2.006e-02  7.036e-05 -285.119   <2e-16 ***
#      band7        2.549e-05  8.969e-05    0.284    0.776    
# ---
#      Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8.606 on 1808277 degrees of freedom
# (13712 observations deleted due to missingness)
# Multiple R-squared:  0.8582,     Adjusted R-squared:  0.8582 
# F-statistic: 1.824e+06 on 6 and 1808277 DF,  p-value: < 2.2e-16

predVCF <- predict(covariates, model=modelVCF, na.rm=TRUE)
predVCF[predVCF < 0] <- 0

diffVCF <- vcfGewata - predVCF
residVCF <- diffVCF^2
rmseVCF <- sqrt(cellStats(residVCF, mean, na.rm=T))

opar <- par(mfrow=c(1,3))
plot(predVCF)
plot(vcfGewata)
plot(diffVCF)

classes <- rasterize(trainingPoly, diffVCF, field='Class')

diffClass <- zonal(diffVCF, classes, fun=mean, na.rm=T)
plot(diffClass)

meanClass <- zonal(residVCF, classes, fun=mean, na.rm=T)
rmseClass <- data.frame(meanClass)
rmseClass$value <- sqrt(rmseClass$value)
rmseClass$class <- levels(trainingPoly$Class)
names(rmseClass) <- c('Zone', 'RMSE', 'Class')

rmseClassSP <- unionSpatialPolygons(trainingPoly, trainingPoly$Class)
rmseClassSPDF <- SpatialPolygonsDataFrame(classPol, rmseClass, match.ID=F)
spplot(rmseClassSPDF, zcol='RMSE')
