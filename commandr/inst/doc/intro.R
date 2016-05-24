### R code from vignette source 'intro.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=72)


###################################################
### code chunk number 2: example-mini
###################################################
library(commandr)

## let's name the first stage 'trim' and specify input types 
setStage("trim", intype = "numeric")
## a class called 'StageTrim' is defined automatically
## then we implemnted protocol 1, 2, 3  for this stage
setProtocol("remove", fun = function(x){x[!is.na(x)]}, parent = "trim")
setProtocol("replace", representation = list(val = "numeric"),
            fun = function(x, val = 0){
              x[is.na(x)] <- val
              x
}, parent = "trim")
setProtocol("range", representation = list(low = "numeric", high = "numeric"), 
            fun = function(x, low = 0, high = Inf) x[x >= low & x <= high & !is.na(x)],
            parent = "trim")

## let's name the second stage 'average' and specify input types 
setStage("average", intype = "numeric")
## implement protocol 1 and 2 for stage 1
setProtocol("mean", fun = mean, parent = "average")
setProtocol("quantile", representation = list(probs = "numeric"),
            fun = quantile, parent = "average")

d <- c(1, 2, 3, NA,  5, 6, 7)

## First pipepine: 1. remove missing value 2. compuate mean
## by default, use default protocol for each stage
p <- Pipeline("trim", "average")
perform(p, d)

## Second pipepine: 1. replcae missing value with 100. 2. compuate mean
## make another pipeline easily
p <- Pipeline(Protocol("trim", "replace", val = 100), 
              "average")

perform(p, d)

## Third pipepine: 1. remove missing value and get value above 2. 2. compuate quantile
p <- Pipeline(Protocol("trim", "range", low = 2),
              Protocol("average", "quantile", probs = 0.75),
              displayName = "Filter and Average")
perform(p, d)


###################################################
### code chunk number 3: pipe-accessor
###################################################
## accessor
inType(p)
outType(p)
parameters(p)
protocol(p, "average")
displayName(p)
## find a protocol via protocol name
findProtocols(p, "average")


###################################################
### code chunk number 4: pipe-coerce
###################################################
# make a new example
setStage("DemoCastN2C", intype = "numeric", outtype = "character")
setProtocol("cast", fun = function(x){
               message("Convert from numeric to character")
               as.character(x)
            },
            parent = "DemoCastN2C")

setStage("DemoCastC2F", intype = "character", outtype = "factor")
setProtocol("cast", fun = function(x){
               message("Convert from character to factor")
               as.factor(x)
            },
            parent = "DemoCastC2F")

setStage("DemoCastF2L", intype = "factor", outtype = "list")
setProtocol("cast", fun = function(x){
               message("Convert from factor to list")
               as.list(x)
            },
            parent = "DemoCastF2L")

d <- 1:3
p <- Pipeline(Protocol("DemoCastN2C"),
              Protocol("DemoCastC2F"),
              Protocol("DemoCastF2L"))
p
perform(p, d)
# subsetting
# convert to a factor
p12 <- p[1:2]
p12
perform(p12, d)

#
p23 <- pipeline(p, intype = "character")
p23
perform(p23, as.character(d))

#
p12 <- head(p, 2)
p12
#or
head(p, outtype = "factor")
head(p, role = "DemoCastC2F")

tail(p, 2)
tail(p, intype = "character")
tail(p, intype = "factor")
tail(p, role = "DemoCastC2F")

#combination
p1 <- Pipeline(Protocol("DemoCastN2C"))
p2 <- Pipeline(Protocol("DemoCastC2F"))
p3 <- Pipeline(Protocol("DemoCastF2L"))
c(p1 ,p2)
p[2] <- p2


###################################################
### code chunk number 5: sessionInfo
###################################################
sessionInfo()


