###
# NADA for R by Lopaka Lee.
#
# Version 1.3-0
# Copyright (2004, 2005, 2006) Lopaka Lee
#
# A S-language software module based on 
# methodologies described by Dennis R. Helsel in his book 
# Nondetects and Data Analysis: Statistics for Censored Environmental Data.
#
# NADA is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
# 
# NADA is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.  You should have received a copy of the GNU General
# Public License along with NADA; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
###

### Definitions common to All sections in this package

## Globals

#.onLoad = function(lib, pkg) 
#{
#    require(methods)
#}

# These probabilities are the used as defaults in methods like quantile.
NADAprobs = c(0.05,0.10,0.25,0.50,0.75,0.90,0.95)

## Generics

setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setGeneric("print", function(x, ...) standardGeneric("print"))

setGeneric("summary", function(object, ...) standardGeneric("summary"))

setGeneric("mean", function(x, ...) standardGeneric("mean"))

setGeneric("sd", function(x, na.rm=FALSE) standardGeneric("sd"))

setGeneric("median", function(x, na.rm=FALSE) standardGeneric("median"))

setGeneric("quantile", function(x, ...) standardGeneric("quantile"))

setGeneric("predict", function(object, ...) standardGeneric("predict"))

setGeneric("pexceed", function(object, ...) standardGeneric("pexceed"))

setGeneric("lines", function(x, ...) standardGeneric("lines"))

setGeneric("boxplot", function(x, ...) standardGeneric("boxplot"))

setGeneric("cor", function(x, y = NULL, use = "all.obs",
          method = c("pearson", "kendall", "spearman")) standardGeneric("cor"))

# LCL and UCL return string representations of 
# the lower and upper conf limits of an object (e.g. "0.95LCL")
setGeneric("LCL", function(x) standardGeneric("LCL"))
setGeneric("UCL", function(x) standardGeneric("UCL"))

## Broken for the time being -- use lines
#setGeneric("abline", 
#           function(a, b, h, v, reg, coef, untf, col, lty, lwd, ...) 
#           standardGeneric("abline"))

setGeneric("residuals", function(object, ...) standardGeneric("residuals"))

setGeneric("coef", function(object, ...) standardGeneric("coef"))

## Maybe in the future transform() could be a generic 
#  Remember that transform is different in R minor versions < 3
#if (as.numeric(version$minor) < 3) {
#    if (!isGeneric("transform"))
#      setGeneric("transform", function(x, ...) standardGeneric("transform"))
#} else {
#    if (!isGeneric("transform"))
#      setGeneric("transform", function(`_data`, ...)
#                 standardGeneric("transform"))
#}


## Classes

setClass("NADAList", "list")

## Methods

#setMethod("summary", signature(), function(x, ...))

setMethod("print", signature("NADAList"), function(x, ...) show(x))

setMethod("show", signature("NADAList"), function(object)
{
    tag = names(object)
    for (i in 1:length(object))
      {
        cat(tag[i], "\n")
        print(object[[i]])
        cat("\n")
      }
})

#-->> BEGIN general utility functions

# Returns par("usr") transformed if in log units
cenpar.usr = 
function(log)
{
    usr = par("usr")
    switch(log,
        xy = (usr = 10^usr),
        x  = (usr[1:2] = 10^usr[1:2]),
        y  = (usr[3:4] = 10^usr[3:4])
    )
    return(usr)
}

##
# split_qual extracts qualifed and unqualifed vectors from a vector
# containing concatenated qualifiying characters and numeric values
# like "<0.5".  Only handles one kind of censoring character/symbol.
splitQual =
function(v, qual.symbol = "<")
{
    v = as.character(v)

    obs = as.numeric(sub(qual.symbol, "", v))
    cen = rep(FALSE, length(obs))
    cen[grep(qual.symbol, v)] = TRUE 

    return(list(obs=obs, cen=cen))
}

## pctCen -- Simple function to save some typing
pctCen =
function(obs, censored)
{
    if (!is.logical(censored)) stop("censored must be logical vector!\n")

    return(100*(length(obs[censored])/length(obs)))
}

#-->> END general utility functions

