## ---- eval=c(2,3), echo=2------------------------------------------------
install.packages("HydeNet")
library(HydeNet)
options(Hyde_fitModel = FALSE)

## ---- fig.width=7, eval=1------------------------------------------------
net <- HydeNetwork(~ wells
                   + pe | wells
                   + d.dimer | pregnant*pe
                   + angio | pe
                   + treat | d.dimer*angio
                   + death | pe*treat)

plot(net)

## ---- echo=FALSE---------------------------------------------------------
net

## ------------------------------------------------------------------------
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
glm(treat ~ d.dimer+angio, data=PE, family="binomial")$coef

## ------------------------------------------------------------------------
xtabs(~PE$pregnant) / nrow(PE)

## ------------------------------------------------------------------------
g1 <- lm(wells ~ 1, data=PE)
g2 <- glm(pe ~ wells, data=PE, family="binomial")
g3 <- lm(d.dimer ~ pe + pregnant, data=PE)
g4 <- xtabs(~ pregnant, data=PE)
g5 <- cpt(angio ~ pe, data=PE)
g6 <- glm(treat ~ d.dimer + angio, data=PE, family="binomial")
g7 <- cpt(death ~ pe + treat, data=PE)

bagOfModels <- list(g1,g2,g3,g4,g5,g6,g7)

bagNet <- HydeNetwork(bagOfModels)
writeNetworkModel(bagNet, pretty=TRUE)

## ------------------------------------------------------------------------
net <- setNode(network = net, node = pregnant,
               nodeType = "dbern", p=.4)
net

## ------------------------------------------------------------------------
net <- setNode(net, wells,
               nodeType = "dnorm", 
               mu = 5, tau = 1 / (1.5^2))

net$nodeType$wells
net$nodeParams$wells

## ------------------------------------------------------------------------
net <- setNode(net, wells,
               nodeType = "dcat",
               pi = vectorProbs(p = c(.3, .6, .1), wells) )
               
net$nodeType$wells
net$nodeParams$wells

## ------------------------------------------------------------------------
net <- setNode(net, wells,
               nodeType = "dcat",
               pi = vectorProbs(p = c(37, 162, 48), wells) )
               
net$nodeType$wells
net$nodeParams$wells

## ------------------------------------------------------------------------
net <- setNode(net, wells,
               nodeType = "dcat",
               pi = "pi.wells[1] <- 0.15; pi.wells[2] <- 0.66; pi.wells[3] <- 0.19")

## ---- eval=FALSE---------------------------------------------------------
#  data(jagsDists, package='HydeNet')
#  jagsDists[,c(1,2,5,6)]

## ---- eval=FALSE---------------------------------------------------------
#  net <- setNode(net, XYZ, nodeType = "dweib", nu=2, lambda=5)

## ---- error=TRUE---------------------------------------------------------
net <- setNode(net, d.dimer, nodeType = "dpois", lambda=-10)

## ------------------------------------------------------------------------
net <- setNode(net, d.dimer, nodeType="dnorm",
               mu=fromFormula(), tau=1/30,  #sigma^2 = 30
               nodeFormula = d.dimer ~ 210 + 29*pregnant + 68*pe)

net$nodeType$d.dimer
net$nodeParams$d.dimer
net$nodeFormula$d.dimer

## ---- eval=FALSE---------------------------------------------------------
#  net <- setNode(net, d.dimer, nodeType="dnorm",
#                 mu="210 + 29*pregnant + 68*pe", tau=1/30)

## ---- eval=FALSE---------------------------------------------------------
#  net <- setNode(net, d.dimer, nodeType="dt",
#                 mu="210 + 29*pregnant + 68*pe", tau=1/20, k=2)

## ---- eval=FALSE---------------------------------------------------------
#  data(jagsFunctions, package='HydeNet')
#  jagsFunctions

## ------------------------------------------------------------------------
equation <- "-6.3 + 0.02*d.dimer + 2.9*angio - 0.005*d.dimer*angio"
net <- setNode(net, treat, nodeType="dbern",
               p=paste("ilogit(", equation, ")"), 
               validate=FALSE)

## ------------------------------------------------------------------------
bagNet$nodeType$d.dimer
bagNet$nodeParams$d.dimer
bagNet$nodeFormula$d.dimer

## ------------------------------------------------------------------------
new.DDimer.Model <- lm(d.dimer ~ pe * pregnant, data=PE)
bagNet <- setNodeModels(bagNet, new.DDimer.Model)

writeNetworkModel(bagNet, pretty=TRUE)

## ------------------------------------------------------------------------
h <- cpt(death ~ pe + treat, data=PE)

## ---- eval=FALSE---------------------------------------------------------
#  h <- inputCPT(test ~ disease)

## ---- eval=FALSE---------------------------------------------------------
#  print(h)

## ---- fig.width=3, eval = -10--------------------------------------------
craps <- HydeNetwork(~ d1 + d2 + diceSum | d1*d2
                       + firstRollOutcome | diceSum)

craps <- setNode(craps, d1, nodeType="dcat",
                 pi = vectorProbs(p = rep(1/6,6), d1),
                 validate = FALSE)
craps <- setNode(craps, d2, nodeType="dcat",
                 pi = vectorProbs(p = rep(1/6,6), d2),
                 validate = FALSE)

craps <- setNode(craps, diceSum, nodeType = "determ",
                 define = fromFormula(),
                 nodeFormula = diceSum ~ di1 + di2)

craps <- setNode(craps, firstRollOutcome, nodeType = "determ",
                 define = fromFormula(),
                 nodeFormula = firstRollOutcome ~ 
                                 ifelse(diceSum < 4 | diceSum > 11, -1,
                                    ifelse(diceSum == 7 | diceSum == 11, 1,0)))

plot(craps)

## ------------------------------------------------------------------------
writeNetworkModel(craps, pretty=TRUE)

## ---- fig.width=7, eval=1------------------------------------------------
net2 <- update(net, . ~ . + newTest | pe
                          + treat | newTest
                          - pregnant)


plot(net2)

## ------------------------------------------------------------------------
net2

