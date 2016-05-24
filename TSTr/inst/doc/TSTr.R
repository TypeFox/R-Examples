## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(TSTr.print_min = 4L, TSTr.print_max = 4L)
library(TSTr)

## ------------------------------------------------------------------------
# Create a tree with the names of the US states
states <- sample(state.name)
stateTree <- newTree(states)

str(stateTree)

## ------------------------------------------------------------------------
# Add some Canada regions to the previous stateTree
regions  <- c("Quebec", "Ontario", "Manitoba", "Saskatchewan", "Alberta", "British Columbia")
US.CanadaTree <- addToTree(stateTree, regions)

# Add one more region
US.CanadaTree <- addWord(US.CanadaTree, "Yukon")

## ------------------------------------------------------------------------
# View the final dimensions of the tree
dimTree(US.CanadaTree)

## ------------------------------------------------------------------------
# Search a specific state
searchWord(US.CanadaTree, "Alabama")
searchWord(US.CanadaTree, "Baltimore")

## ------------------------------------------------------------------------
# Complete strings: States and regions that begin with "A" and "Al"
completeWord(US.CanadaTree, "A")
completeWord(US.CanadaTree, "Al")

## ------------------------------------------------------------------------
# Peter Norvig spell corrector.
PNcheck(US.CanadaTree, "Conecticut")
PNcheck(US.CanadaTree, "Sorth Carolina", useUpper = TRUE)

## ------------------------------------------------------------------------
# Symmetric Delete pre-calculation step.
US.CanadaDT <- SDkeeper(states, 2)
US.CanadaTree <- SDkeeper(states, 1, useTST = TRUE)

## ------------------------------------------------------------------------
# Symmetric Delete spell checking.
SDcheck(US.CanadaDT, "rkansas", summarize = TRUE)
SDcheck(US.CanadaTree, "Texas2")

