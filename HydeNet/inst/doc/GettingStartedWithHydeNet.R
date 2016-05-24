## ---- eval=FALSE---------------------------------------------------------
#  install.packages("HydeNet")

## ---- eval=FALSE---------------------------------------------------------
#  setRepositories(ind=1:2)
#  devtools::install_github("nutterb/HydeNet")

## ------------------------------------------------------------------------
library(HydeNet)

mtcars2 <- transform(mtcars,
                     cyl = factor(cyl),
                     gear=factor(gear),
                     am = factor(am))

## ------------------------------------------------------------------------
carNet <- HydeNetwork(~ cyl
                      + disp | cyl
                      + hp | disp
                      + wt
                      + gear
                      + mpg | disp*hp*wt*gear,
                      data=mtcars2)

## ---- fig.width = 5, eval=FALSE------------------------------------------
#  plot(carNet)

## ------------------------------------------------------------------------
HydeNet:::writeJagsModel(carNet, node = "cyl")
HydeNet:::writeJagsModel(carNet, node = "mpg")

writeNetworkModel(carNet, pretty = TRUE)

## ------------------------------------------------------------------------
carNet1 <- compileJagsModel(carNet)

## ------------------------------------------------------------------------
carNet2 <- compileJagsModel(carNet, data = list(cyl = "8") )

## ------------------------------------------------------------------------
carNet3 <- compileJagsModel(carNet, data=list(cyl="8"))

## ---- fig.width=6, fig.height=6------------------------------------------
post1 <- HydePosterior(carNet1,
                       variable.names = c("cyl","hp","mpg"),
                       n.iter = 10000,
                       bind=FALSE)

post2 <- HydePosterior(carNet2,
                       variable.names = c("cyl","hp","mpg"),
                       n.iter = 10000,
                       bind = FALSE)

str(post1, max.level = 3)
plot(post1$codas[,c("hp","mpg")])

## ------------------------------------------------------------------------
bp1 <- bindPosterior(post1)
bp2 <- bindPosterior(post2)

head(bp1)
head(bp2) #notice cyl = "8" for all samples

plot(density(bp1$hp), ylim=c(0,0.06), main = "hp");
lines(density(bp2$hp), col="red", lty=5)

plot(density(bp1$mpg), main = "mpg");
lines(density(bp2$mpg), col="red", lty=5)

## ---- eval=FALSE---------------------------------------------------------
#  options(Hyde_fitModel=FALSE)

## ---- echo=FALSE---------------------------------------------------------
data(PE, package='HydeNet')
autoNet <- HydeNetwork(~ wells
                       + pe | wells
                       + d.dimer | pregnant*pe
                       + angio | pe
                       + treat | d.dimer*angio
                       + death | pe*treat,
                       data = PE)
autoNet <- setNode(autoNet, treat,
                   nodeFormula = treat ~ poly(d.dimer, 2) + angio,
                   p = fromData())

## ---- echo=FALSE---------------------------------------------------------
print(autoNet, treat)

## ---- echo=FALSE---------------------------------------------------------
library(HydeNet)
options(Hyde_fitModel=TRUE)
data(PE, package='HydeNet')
autoNet <- HydeNetwork(~ wells
                       + pe | wells
                       + d.dimer | pregnant*pe
                       + angio | pe
                       + treat | d.dimer*angio
                       + death | pe*treat,
                       data = PE)
autoNet <- setNode(autoNet, treat,
                   nodeFormula = treat ~ poly(d.dimer, 2) + angio,
                   p = fromData())
print(autoNet, treat)

## ------------------------------------------------------------------------
library(HydeNet)
data(PE, package='HydeNet')
autoNet <- HydeNetwork(~ wells
                       + pe | wells
                       + d.dimer | pregnant*pe
                       + angio | pe
                       + treat | d.dimer*angio
                       + death | pe*treat,
                       data = PE)
writeNetworkModel(autoNet, pretty=TRUE)

## ------------------------------------------------------------------------
library(HydeNet)
options(Hyde_maxDigits=2)
data(PE, package='HydeNet')
autoNet <- HydeNetwork(~ wells
                       + pe | wells
                       + d.dimer | pregnant*pe
                       + angio | pe
                       + treat | d.dimer*angio
                       + death | pe*treat,
                       data = PE)
writeNetworkModel(autoNet, pretty=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  > #* Object based on a list of models
#  > g1 <- lm(wells ~ 1, data=PE)
#  > g2 <- glm(pe ~ wells, data=PE, family="binomial")
#  > g3 <- lm(d.dimer ~ pe + pregnant, data=PE)
#  > g4 <- xtabs(~ pregnant, data=PE)
#  > g5 <- glm(angio ~ pe, data=PE, family="binomial")
#  > g6 <- glm(treat ~ d.dimer + angio, data=PE, family="binomial")
#  > g7 <- glm(death ~ pe + treat, data=PE, family="binomial")
#  
#  > bagOfModels <- list(g1,g2,g3,g4,g5,g6,g7)
#  > bagNet <- HydeNetwork(bagOfModels)
#  
#  > #* Time to print bagNet implicitly
#  > a <- Sys.time()
#  > bagNet
#  > b <- Sys.time()
#  > b-a
#  Time difference of 33.53736 secs
#  
#  > #* Time to print bagNet explicitly
#  > a <- Sys.time()
#  > print(bagNet)
#  > b <- Sys.time()
#  > b-a
#  Time difference of 0 secs

