## ---- echo=FALSE, eval=TRUE----------------------------------------------
set.seed(10)
suppressMessages(library(spiders))
Predators <- 5*c(11,22,33,44,77)
Traps <- Predators
PreySpecies <- 3
Times <- 5
ST <- Times*PreySpecies
g <- matrix(sqrt(2), nrow=Times, ncol=PreySpecies)     # gamma
l <- matrix(0.5*sqrt(2), nrow=Times, ncol=PreySpecies) # c
fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=F)
Traps <- fdata$caught
Traps$adj <- 1
Spiders <- fdata$eaten

## ------------------------------------------------------------------------
head(Traps)

## ------------------------------------------------------------------------
tail(Traps)

## ------------------------------------------------------------------------
all.equal(unique(Traps[,'time']), unique(Spiders[,'time']))

## ------------------------------------------------------------------------
prefs <- predPref(Spiders, Traps, hypotheses=c(null="C", alt="Cst"))

## ------------------------------------------------------------------------
summary(prefs)

## ------------------------------------------------------------------------
prefs$alt$c
length(prefs$alt$c)

## ------------------------------------------------------------------------
b <- c(1/2, 1/2, 0, 0, 0, -1/2, -1/2, rep(0, 8))

## ------------------------------------------------------------------------
testC(prefs, b, sig.level=0.8)

