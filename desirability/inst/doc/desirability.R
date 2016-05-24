### R code from vignette source 'desirability.Rnw'

###################################################
### code chunk number 1: makeDesirePlots
###################################################
library(desirability)
pdf("dMax.pdf", width = 7, height = 4)
   par(mfrow=c(1,3))
   
   plot(dMax(1, 3), nonInform = FALSE)
   
   plot(dMax(1, 3, 5),  TRUE, col = 2, nonInform = FALSE)
   text(2.73, .3, "scale = 5", col = 2) 

   plot(dMax(1, 3, 1/5),  TRUE, col = 4, nonInform = FALSE) 
   text(1.3, .8, "scale = .2", col = 4)   
   
   text(2, .62, "scale = 1", col = 1)   
   title("(a)")   

   plot(dMin(1, 3), nonInform = FALSE)
   
   plot(dMin(1, 3, 5),  TRUE, col = 2, nonInform = FALSE)
   text(1.5, 0.1, "scale = 5", col = 2)   
  
   plot(dMin(1, 3, 1/5),  TRUE, col = 4, nonInform = FALSE)  
   text(1.5, 1, "scale = .2", col = 4)  
    
   text(2, .62, "scale = 1", col = 1)   
   title("(b)")   

   plot(dTarget(1, 2, 3), nonInform = FALSE)
   
   plot(dTarget(1, 2, 3, 5),  TRUE, col = 2, nonInform = FALSE)
   text(1.9, .1, "lowScale = 5", col = 2) 

   plot(dTarget(1, 2, 3, 1/5),  TRUE, col = 4, nonInform = FALSE)      
   text(1.3, 0.9, "lowScale = .2", col = 4)   
   
   text(1.35, 0.6, "lowScale = 1", col = 1)   
   title("(c)")   
dev.off()
par(mfrow=c(1,1))


###################################################
### code chunk number 2: responseEq
###################################################
conversionPred <- function(x) 81.09 + 1.0284 * x[1] + 4.043 * x[2] + 6.2037 * x[3] - 
   1.8366 * x[1]^2 + 2.9382 * x[2]^2 - 5.1915 * x[3]^2 +
   2.2150 * x[1] * x[2] + 11.375 * x[1] * x[3] - 3.875 * x[2] * x[3]

activityPred <- function(x) 59.85 + 3.583 * x[1] + 0.2546 * x[2] + 2.2298 * x[3] + 
   0.83479 * x[1]^2 + 0.07484 * x[2]^2 + 0.05716 * x[3]^2 -
   0.3875 * x[1] * x[2] - 0.375 * x[1] * x[3] + 0.3125 * x[2] * x[3]


###################################################
### code chunk number 3: setupData
###################################################
plotGrid <- expand.grid(time = seq(-1.7, 1.7, length = 50), 
                        temperature = seq(-1.7, 1.7, length = 4), 
                        catalyst = seq(-1.7, 1.7, length = 50))
plotGrid$conversionPred <- apply(plotGrid[, 1:3], 1, conversionPred)
plotGrid$activityPred <- apply(plotGrid[, 1:3], 1, activityPred)


###################################################
### code chunk number 4: conversionPlot
###################################################
library(lattice)
textInfo <- trellis.par.get("add.text")
textInfo$cex <- .7
trellis.par.set("add.text", textInfo)
print(contourplot(conversionPred ~ time + catalyst|temperature, plotGrid, aspect = 1, as.table = TRUE))


###################################################
### code chunk number 5: activityPlot
###################################################
print(contourplot(activityPred ~ time + catalyst|temperature, plotGrid, aspect = 1, as.table = TRUE))


###################################################
### code chunk number 6: makeDesire
###################################################
conversionD <- dMax(80, 97)
activityD <- dTarget(55, 57.5, 60)
predOutcomes <- c(conversionPred(c(0,0,0)), activityPred(c(0,0,0)))
print(predOutcomes)
predict(conversionD, predOutcomes[1])
predict(activityD, predOutcomes[2])


###################################################
### code chunk number 7: makeDesire
###################################################
overallD <- dOverall(conversionD, activityD)
print(overallD)
predict(overallD, predOutcomes)


###################################################
### code chunk number 8: setupDPlots
###################################################
dValues <- predict(overallD, plotGrid[, 4:5], all = TRUE)
plotGrid <- cbind(plotGrid, dValues)


###################################################
### code chunk number 9: conversionPlot2
###################################################
print(contourplot(D1 ~ time + catalyst|temperature, plotGrid, aspect = 1, as.table = TRUE))


###################################################
### code chunk number 10: activityPlot2
###################################################
print(contourplot(D2 ~ time + catalyst|temperature, plotGrid, aspect = 1, as.table = TRUE))


###################################################
### code chunk number 11: overallPlot
###################################################
print(contourplot(Overall ~ time + catalyst|temperature, plotGrid, aspect = 1, as.table = TRUE))


###################################################
### code chunk number 12: optFunction
###################################################
rsmOpt <- function(x, dObject, space = "square")
{
  conv <- conversionPred(x)
  acty <- activityPred(x)

  out <- predict(dObject, data.frame(conv = conv, acty = acty))

  if(space == "circular")
    {
      if(sqrt(sum(x^2)) > 1.682) out <- 0
    } else if(space == "square") if(any(abs(x) > 1.682)) out <- 0
  out
}


###################################################
### code chunk number 13: squareRegion
###################################################
searchGrid <- expand.grid(time = seq(-1.5, 1.5, length = 5),
                          temperature = seq(-1.5, 1.5, length = 5),
                          catalyst = seq(-1.5, 1.5, length = 5))

for(i in 1:dim(searchGrid)[1])
{
  tmp <- optim(as.vector(searchGrid[i,]), 
               rsmOpt, 
               dObject = overallD,
               space = "square",
               control = list(fnscale = -1))
  if(i == 1)
    {
      best <- tmp
    } else {
      if(tmp$value > best$value) best <- tmp   
    }
}
print(best) 


###################################################
### code chunk number 14: circularRegion
###################################################
for(i in 1:dim(searchGrid)[1])
{
  tmp <- optim(as.vector(searchGrid[i,]), 
               rsmOpt, 
               space = "circular",      
               dObject = overallD,      
               control = list(fnscale = -1))
  if(i == 1)
    {
      best <- tmp
    } else {
      if(tmp$value > best$value) best <- tmp   
    }
}
print(best) 


###################################################
### code chunk number 15: logistic
###################################################
foo <- function(u) 1/(1+exp(-u))
xInput <- seq(-5,5, length = 20)
logisticD <- dArb(xInput, foo(xInput))


###################################################
### code chunk number 16: logisticPlot
###################################################
pdf("logistic.pdf", width = 5, height = 4)
   plot(logisticD)
dev.off()


###################################################
### code chunk number 17: cats
###################################################
values <- c("value1" = .1, "value2" = .9, "value3" = .2)
groupedDesirabilities <- dCategorical(values)
groupedDesirabilities


###################################################
### code chunk number 18: box
###################################################
pdf("box.pdf", width = 5, height = 4)
plot(dBox(-1.682, 1.682), nonInform = FALSE)
dev.off()


###################################################
### code chunk number 19: groups
###################################################
pdf("groups.pdf", width = 5, height = 4)
plot(groupedDesirabilities, nonInform = FALSE)
dev.off()


