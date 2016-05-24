## ----set-options, echo=FALSE, cache=FALSE-----------------------------------------------------------------------------
options(width=120)
#opts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")
#read_chunk("some/script/I/want/to/load.R")

## ---- eval = FALSE----------------------------------------------------------------------------------------------------
#  library(R6)
#  vignette('Introduction', package = 'R6')

## ---- eval = FALSE----------------------------------------------------------------------------------------------------
#  library(data.tree)
#  vignette(package = 'data.tree')

## ---- echo=FALSE------------------------------------------------------------------------------------------------------
ahpFile <- system.file("extdata", "car.ahp", package="ahp")
cat(readChar(ahpFile, file.info(ahpFile)$size))

## ---------------------------------------------------------------------------------------------------------------------
library(ahp)
ahpFile <- system.file("extdata", "car.ahp", package="ahp")
carAhp <- Load(ahpFile)

## ---------------------------------------------------------------------------------------------------------------------
library(data.tree)
print(carAhp, filterFun = isNotLeaf)

## ---------------------------------------------------------------------------------------------------------------------
Calculate(carAhp)
print(carAhp, priority = function(x) x$parent$priority["Total", x$name])

## ---- eval = FALSE----------------------------------------------------------------------------------------------------
#  Visualize(carAhp)

## ---------------------------------------------------------------------------------------------------------------------

Analyze(carAhp)


## ---------------------------------------------------------------------------------------------------------------------
AnalyzeTable(carAhp)

## ---------------------------------------------------------------------------------------------------------------------
AnalyzeTable(carAhp, 
             variable = "priority", 
             sort = "orig",
             pruneFun = function(node, decisionMaker) PruneByCutoff(node, decisionMaker, 0.05),
             weightColor = "skyblue",
             consistencyColor = "red",
             alternativeColor = "green")

## ---- comment = NA----------------------------------------------------------------------------------------------------
ahpFile <- system.file("extdata", "vacation.ahp", package="ahp")
cat(readChar(ahpFile, file.info(ahpFile)$size))

## ---------------------------------------------------------------------------------------------------------------------
ahpFile <- system.file("extdata", "vacation.ahp", package="ahp")
vacationAhp <- Load(ahpFile)
Calculate(vacationAhp)


## ---- eval = FALSE----------------------------------------------------------------------------------------------------
#  Visualize(vacationAhp)

## ---------------------------------------------------------------------------------------------------------------------
AnalyzeTable(vacationAhp)

## ---------------------------------------------------------------------------------------------------------------------
AnalyzeTable(vacationAhp, decisionMaker = "Dad")

## ---------------------------------------------------------------------------------------------------------------------
AnalyzeTable(vacationAhp, decisionMaker = "Mom")

## ---------------------------------------------------------------------------------------------------------------------
AnalyzeTable(vacationAhp, decisionMaker = "Kid")

