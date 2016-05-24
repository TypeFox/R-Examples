## ----simulate CA---------------------------------------------------------
library(gjam)
sim <- gjamSimData(n=500,S=10,q=3,typeNames='CA')
summary(sim)

## ----show formula--------------------------------------------------------
sim$formula

## ----plot simulated y, fig.show = "hold", fig.width=6--------------------
par(bty='n', mfrow=c(1,2))
h <- hist(c(-1,sim$y),nclass=50,plot=F)
plot(h$counts,h$mids,type='s')
plot(sim$w,sim$y,cex=.2)

## ----fit CA data---------------------------------------------------------
# a few iterations
modelList <- list(ng=100, burnin=10, typeNames = sim$typeNames)
out       <- gjamGibbs(sim$formula, sim$xdata, sim$y, modelList)
summary(out)

## ----summary of chains---------------------------------------------------
summary(out$chains)

## ----summary of fitted model, eval=FALSE---------------------------------
#  summary(out$modelSummary)

## ----show classes--------------------------------------------------------
out$modelSummary$classBySpec

## ----plot CA data, eval=FALSE--------------------------------------------
#  sim       <- gjamSimData(n=500,S=10,q=3,typeNames='CA')
#  modelList <- list(ng=2000, burnin=500, typeNames = sim$typeNames)
#  out       <- gjamGibbs(sim$formula, sim$xdata, sim$y, modelList)
#  plotPars  <- list(trueValues = sim$trueValues,width=3,height=2,
#                    CLUSTERPLOTS=T, SMALLPLOTS=F)
#  fit       <- gjamPlot(output = out, plotPars)

## ----example output, fig.show = "hold", fig.width = 6--------------------
par(bty='n', mfrow=c(1,3))
plot(sim$trueValues$beta, out$modelSummary$betaMu)
plot(sim$trueValues$corSpec, out$modelSummary$corMu)
plot(sim$y,out$modelSummary$yMu, cex=.2)

## ----effort simulation---------------------------------------------------
S       <- 5                             
n       <- 50
effort  <- list( columns = 1:S, values = round(runif(n,.5,5),1) )
sim     <- gjamSimData(n,S,q=5,typeNames='DA',effort=effort)
effort

## ----w vs y, bty='n', fig.width=6----------------------------------------
plot(sim$w,sim$y, cex=.2)

## ----including effort, bty='n', fig.width=6------------------------------
plot(sim$w*effort$values, sim$y, cex=.2)

## ----fitting, eval=F-----------------------------------------------------
#  S         <- 5
#  n         <- 500
#  effort    <- list( columns = 1:S, values = round(runif(n,.5,5),1) )
#  sim       <- gjamSimData(n,S,q=5,typeNames='DA',effort=effort)
#  modelList <- list(ng=2000, burnin=500, typeNames = sim$typeNames,
#                    effort = effort)
#  out       <- gjamGibbs(sim$formula, sim$xdata, sim$y, modelList)
#  plotPars  <- list(trueValues = out$trueValues,width=3,height=2)
#  gjamPlot(output = outDA, plotPars)

## ----composition data, eval=FALSE----------------------------------------
#  sim       <- gjamSimData(n=1000,S=8,q=5,typeNames='CC')
#  modelList <- list(ng=2000, burnin=500, typeNames = sim$typeNames)
#  out       <- gjamGibbs(sim$formula,sim$xdata, sim$y, modelList)
#  plotPars  <- list(trueValues = sim$trueValues, width=3, height=2,
#                    CLUSTERPLOTS=T, SMALLPLOTS=F)
#  gjamPlot(output = out, plotPars)

## ----trait summary-------------------------------------------------------
data(forestTraits)
summary(forestTraits)

## ----trait censor 1------------------------------------------------------
y   <- forestTraits$y
tmp <- gjamCensorY(values=-999,intervals=cbind( c(-Inf,-999) ),
                         y = y, whichcol=1:6)
censor <- list('CA' = tmp$censor)

## ----trait censor 2------------------------------------------------------
tmp    <- gjamCensorY(values=c(0,1), intervals=cbind( c(-Inf,0),c(1,Inf) ),
                       y = y, whichcol=13:14)
censor <- append(censor,list('CA' = tmp$censor))

## ----trait analysis, eval=FALSE------------------------------------------
#  modelList <- list(ng=3000, burnin=500, typeNames = forestTraits$typeNames,
#                    xfactors='soilFactor',censor=censor)
#  out       <- gjamGibbs(~ temp + moisture*(deficit + soilFactor),
#                         xdata = forestTraits$xdata, y = y, modelList = modelList)

## ----trait plots, eval=FALSE---------------------------------------------
#  plotPars  <- list(width=3,height=2, CLUSTERPLOTS=T)
#  gjamPlot(output = out, plotPars)

## ----ordinal, eval=FALSE-------------------------------------------------
#  sim       <- gjamSimData(n=1000,S=10,q=3,typeNames='OC')
#  modelList <- list(ng = 2000, burnin = 500, typeNames = sim$typeNames)
#  out       <- gjamGibbs(sim$formula, sim$xdata, sim$y, modelList)

## ----ordinal partition, eval=FALSE---------------------------------------
#  keep <- strsplit(colnames(out$modelSummary$cutMu),'C-') #only saved columns
#  keep <- matrix(as.numeric(unlist(keep)),ncol=2,byrow=T)[,2]
#  plot(sim$trueValues$cuts[,keep],out$modelSummary$cutMu)

## ----ordinal plots, eval=FALSE-------------------------------------------
#  plotPars  <- list(trueValues = sim$trueValues,width=3,height=2)
#  fit       <- gjamPlot(output = out, plotPars)

## ----many types----------------------------------------------------------
typeNames <- c('OC','OC','OC','CC','CC','CC',
               'CC','CC','CA','CA','PA','PA')         
sim       <- gjamSimData(n=1000,S=length(typeNames),q=3,typeNames=typeNames)
# a few iterations
modelList <- list(ng = 100, burnin = 20, typeNames = sim$typeNames)
out       <- gjamGibbs(sim$formula, sim$xdata, sim$y, modelList)
tmp <- data.frame(sim$typeNames,out$modelSummary$classBySpec[,1:10])
print(tmp)

## ----mixed analysis, eval=FALSE------------------------------------------
#  # repeat with ng = 2000, burnin = 500, then plot
#  modelList <- list(ng = 2000, burnin = 500, typeNames = sim$typeNames)
#  out       <- gjamGibbs(sim$formula, sim$xdata, sim$y, modelList)
#  plotPars  <- list(trueValues = sim$trueValues,width=3,height=2)
#  gjamPlot(output = out, plotPars)

## ----simulate missing data-----------------------------------------------
sim <- gjamSimData(n=500,S=5,q=3,typeNames='OC', nmiss = 20)
which(is.na(sim$xdata),arr.ind=T)

## ----holdouts, eval=FALSE------------------------------------------------
#  sim       <- gjamSimData(n=1000,S=5,q=3,typeNames='CA', nmiss = 20)
#  modelList <- list(ng = 2000, burnin = 500, typeNames = sim$typeNames, holdoutN=50)
#  out       <- gjamGibbs(sim$formula, sim$xdata, sim$y, modelList)
#  
#  plot(out$x[out$missingIndex],out$modelSummary$xpredMu[out$missingIndex])
#  title('missing in x'); abline(0,1)
#  plot(out$x[out$holdoutIndex,-1],out$modelSummary$xpredMu[out$holdoutIndex,-1])
#  title('holdouts in x'); abline(0,1)
#  plot(out$y[out$holdoutIndex,],out$modelSummary$yMu[out$holdoutIndex,])
#  title('holdouts in y'); abline(0,1)

## ----cluster plots, eval=FALSE-------------------------------------------
#  plotPars  <- list(trueValues = sim$trueValues, width=3, height=2, CLUSTERPLOTS=T)
#  gjamPlot(output = out, plotPars)

