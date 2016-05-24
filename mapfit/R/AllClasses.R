#### PH Class (S4)

## ph: general PH
## cf1: canonical form 1
## cf1srm: a special class for truncated data

setClass("ph", representation(size="numeric", alpha="numeric", Q="Matrix", xi="numeric", df="numeric"))
setClass("cf1", representation(rate="numeric"), contains="ph")
## setClass("cf1srm", representation(omega="numeric"), contains="cf1")

## hyper Erlang 

setClass("herlang", representation(size="numeric", mixrate="numeric", shape="numeric", rate="numeric"))

ph.moment <- function(k, ph, ...) { cat("Warning: Do not call me. Please check a class.\n") }
setGeneric("ph.moment")

## data class
setClass("phdata", representation(size="numeric", data="data.frame"))
setClass("phdata.wtime", representation(diff="numeric"), contains="phdata")
setClass("phdata.group", contains="phdata")

### map: general MAP
setClass("map", representation(size="numeric", alpha="numeric", D0="Matrix", D1="Matrix", df="numeric"))
setClass("gmmpp", contains="map")

setClass("hmm", representation(size="numeric", alpha="numeric", P="Matrix"))
setClass("erhmm", representation(shape="numeric", rate="numeric"), contains="hmm")

setClass("mapdata", representation(size="numeric", data="data.frame"))
setClass("mapdata.time", contains="mapdata")
setClass("mapdata.group", contains="mapdata")


## generic methods

### for data
mapfit.mean <- function(x, ...) { cat("Warning: Do not call me. Please check a class.\n") }
setGeneric("mapfit.mean")

mean <- function(x, ...) { cat("Warning: Do not call me. Please check a class.\n") }
setGeneric("mapfit.mean")

## initial values
emfit.init <- function(model, data, verbose, ...) { cat("Warning: Do not call me. Please check a class.\n") }
setGeneric("emfit.init")

## estep
emfit.estep <- function(model, data, ...) { cat("Warning: Do not call me. Please check a class.\n") }
setGeneric("emfit.estep")

## mstep
emfit.mstep <- function(model, eres, data, ...) { cat("Warning: Do not call me. Please check a class.\n") }
setGeneric("emfit.mstep")

## degrees of freedom
emfit.df <- function(model, ...) { cat("Warning: Do not call me. Please check a class.\n") }
setGeneric("emfit.df")

## print
emfit.print <- function(model, ...) { cat("Warning: Do not call me. Please check a class.\n") }
setGeneric("emfit.print")

setMethod("print", signature(x = "ph"), function(x, ...) emfit.print(x, ...))
setMethod("print", signature(x = "cf1"), function(x, ...) emfit.print(x, ...))
setMethod("print", signature(x = "herlang"), function(x, ...) emfit.print(x, ...))
setMethod("print", signature(x = "map"), function(x, ...) emfit.print(x, ...))
setMethod("print", signature(x = "erhmm"), function(x, ...) emfit.print(x, ...))


