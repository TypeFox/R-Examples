## ---- echo = FALSE, results = "hide", message = F, warning = F-----------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE)
## Initialize geoprocessing environment
library(RSAGA)
data("landslides")

## ---- echo = FALSE, results = "hide", message = F, warning = F, eval = FALSE----
#  env <- rsaga.env(path = "C:/SAGA-GIS/saga_2.2.0_x64")
#  write.sgrd(data = dem, file = "dem", header = dem$header)

## ---- collapse = TRUE----------------------------------------------------
rsaga.env()

## ---- collapse=TRUE, eval = FALSE----------------------------------------
#  env <- rsaga.env(path = "C:/SAGA-GIS/saga_2.2.0_x64")
#  env$version

## ---- results = "hide", eval = FALSE-------------------------------------
#  rsaga.geoprocessor(lib = "ta_morphometry", module = "Slope, Aspect, Curvature",
#                     param = list(ELEVATION = paste(getwd(),"/dem.sgrd", sep = ""),
#                                  SLOPE = paste(getwd(),"/slope.sgrd", sep = "")),
#                     env = env)

## ---- results = "hide", eval = FALSE-------------------------------------
#  rsaga.slope(in.dem = "dem", out.slope = "slope", env = env)

## ---- eval = FALSE-------------------------------------------------------
#  help(rsaga.slope)

## ---- results = "hide", eval = FALSE-------------------------------------
#  # A character string of available libraries:
#  rsaga.get.libraries(path = env$modules)
#  
#  # A list of modules in a library:
#  rsaga.get.modules(libs = "ta_morphometry", env = env)

## ---- results = "hide", eval=FALSE---------------------------------------
#  rsaga.get.usage(lib = "ta_morphometry", module = "Slope, Aspect, Curvature", env = env)
#  
#  # Compare module parameters between versions:
#  rsaga.get.usage(lib = "ta_morphometry", module = "Slope, Aspect, Curvature",
#                  env = rsaga.env(path = "C:/SAGA-GIS/saga_2.1.0_x64"))
#  
#  rsaga.get.usage(lib = "ta_morphometry", module = "Slope, Aspect, Curvature",
#                  env = rsaga.env(path = "C:/SAGA-GIS/saga_2.2.0_x64"))

## ------------------------------------------------------------------------
data(landslides)

## ---- eval = FALSE-------------------------------------------------------
#  write.sgrd(data = dem, file = "dem", header = dem$header,
#             env = env)  # write.sgrd and read.sgrd use SAGA, and should specify 'env'

## ---- eval = FALSE-------------------------------------------------------
#  write.ascii.grid(data = dem, file = "dem", header = dem$header)

## ---- eval = FALSE-------------------------------------------------------
#  rsaga.slope(in.dem = "C:/InData/dem", out.slope = "C:/OutData/slope", env = env)

## ---- collapse = TRUE, results = "hide", warning = FALSE, eval = FALSE----
#  # By individual function calls:
#  rsaga.slope("dem", "slope", method = "poly2zevenbergen", env = env)
#  rsaga.plan.curvature("dem", "cplan", method = "poly2zevenbergen", env = env)
#  rsaga.profile.curvature("dem", "cprof", method = "poly2zevenbergen", env = env)
#  
#  # By one function that calculates each of the terrain parameters:
#  rsaga.slope.asp.curv("dem", out.slope = "slope",
#                       out.cprof = "cprof", out.cplan = "cplan",
#                       method = "poly2zevenbergen",
#                       env = env)

## ---- collapse = TRUE, results = "hide", eval = FALSE--------------------
#  rsaga.topdown.processing("dem", out.carea = "carea", method = "mfd", env = env)

## ---- collapse = TRUE, results = "hide", eval = FALSE--------------------
#  rsaga.grid.calculus(in.grids = "carea", out.grid = "log10_carea",
#                      formula = ~ log(a), env = env)

## ---- eval = FALSE-------------------------------------------------------
#  # Pick grid values; add to landslides data.frame:
#  landslides <- pick.from.saga.grid(landslides, "slope", varname = "slope", env = env,
#                                    X.name = "x", Y.name = "y")
#  landslides <- pick.from.saga.grid(landslides, "cprof", varname = "cprof", env = env,
#                                    X.name = "x", Y.name = "y")
#  landslides <- pick.from.saga.grid(landslides, "cplan", varname = "cplan", env = env,
#                                    X.name = "x", Y.name = "y")
#  landslides <- pick.from.saga.grid(landslides, "carea", varname = "carea", env = env,
#                                    X.name = "x", Y.name = "y")

## ---- results = "hide", eval = FALSE-------------------------------------
#  rsaga.sgrd.to.esri(in.sgrds = c("slope", "cprof", "cplan", "carea"),
#                     out.grids = c("slope", "cprof", "cplan", "carea"),
#                     out.path = getwd(), format = "ascii", env = env)

## ---- eval = FALSE-------------------------------------------------------
#  my.trafo <- function(x) {
#     x$log.carea <- log10(x$carea)
#     return(x)
#  }
#  
#  landslides <- my.trafo(landslides)

## ---- eval = FALSE-------------------------------------------------------
#  library(gam)
#  
#  fit <- gam(lslpts ~ s(slope) + cprof + s(cplan) + s(log.carea),
#             data = landslides, family = binomial)
#  
#  multi.local.function(in.grids = c("slope", "cprof", "cplan", "carea"),
#     in.varnames = c("slope", "cprof", "cplan", "carea"),
#     out.varnames = "lsl_pred",
#     fun = grid.predict,
#     control.predict = list(type = "response"),
#     trafo = my.trafo,
#     fit = fit)

