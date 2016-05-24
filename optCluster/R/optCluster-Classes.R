########### optCluster Class Definitions #############

## Change raggr to S4 class
setOldClass("raggr")
setClass("optCluster",representation(inputData = "matrix",
									clVal = "clValid",
									ranksWeights = "list",
                                  	rankAgg = "raggr"))