## ----echo=FALSE,results='hide'-------------------------------------------
library("proftools")
options(keep.source = TRUE)
knitr::opts_chunk$set(size = "small", fig.align = "center")

## ------------------------------------------------------------------------
srcfile <- system.file("samples", "bootlmEx.R", package = "proftools")

## ----eval = FALSE--------------------------------------------------------
#  profout <- tempfile()
#  Rprof(file = profout, gc.profiling = TRUE, line.profiling = TRUE)
#  source(srcfile)
#  Rprof(NULL)
#  pd <- readProfileData(profout)
#  unlink(profout)

## ----cache = TRUE--------------------------------------------------------
pd <- profileExpr(source(srcfile))

## ------------------------------------------------------------------------
head(funSummary(pd), 10)

## ------------------------------------------------------------------------
head(funSummary(pd, srclines = FALSE), 10)

## ------------------------------------------------------------------------
head(callSummary(pd), 10)

## ------------------------------------------------------------------------
srcSummary(pd)

## ------------------------------------------------------------------------
hotPaths(pd, total.pct = 10.0)

## ------------------------------------------------------------------------
filteredPD <- filterProfileData(pd, select = "withVisible", skip = 4)

## ------------------------------------------------------------------------
hotPaths(filteredPD, total.pct = 10.0) 

## ------------------------------------------------------------------------
glmPD <- filterProfileData(filteredPD, focus = "glm")
hotPaths(glmPD, total.pct = 5.0)

## ----fullCallGraph, fig.cap = "Full call graph of profile data."---------
plotProfileCallGraph(pd)

## ----filteredCallGraph, fig.cap = "Call graph for \\code{glm.fit} call."----
plotProfileCallGraph(filterProfileData(pd, focus = "glm.fit"))

## ----printProfileCallGraph, eval=FALSE-----------------------------------
#  printProfileCallGraph(filterProfileData(pd, focus = "glm.fit"))

## ----echo = FALSE, comment = NA------------------------------------------
printProfileCallGraph(filterProfileData(pd, focus = "glm.fit"))

## ----flameGraph, out.width = "4in", fig.cap = "Flame graph visualizing hot paths for the full profile data."----
flameGraph(pd)

## ----filteredFlameGraph, out.width = "4in", fig.cap = "Flame graph of the filtered profile data."----
flameGraph(filteredPD)

## ----timeGraph, out.width = "4in", fig.cap = "Time graph of the full profile data."----
flameGraph(pd, order = "time")

## ----eval = FALSE--------------------------------------------------------
#  fg <- flameGraph(pd)
#  identify(fg)

## ----calleeTreeMap, out.width = "4in", fig.cap = "Call tree map of the full profile data."----
calleeTreeMap(pd)

