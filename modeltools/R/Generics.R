#
#  Copyright (C) 2006 Torsten Hothorn, Friedrich Leisch
#  $Id: Generics.R 4119 2008-09-22 14:49:39Z leisch $
#



## generics used in flexmix and flexclust which may also be useful in
## other packages

setGeneric("ICL", function(object, ...) standardGeneric("ICL"))

setGeneric("KLdiv", function(object, ...) standardGeneric("KLdiv"))

setGeneric("Lapply", function(object, FUN, ...) standardGeneric("Lapply"))

setGeneric("clusters", function(object, newdata, ...) standardGeneric("clusters"))

setGeneric("getModel", function(object, ...) standardGeneric("getModel"))

setGeneric("parameters", function(object, ...) standardGeneric("parameters"))

setGeneric("posterior", function(object, newdata, ...) standardGeneric("posterior"))

setGeneric("prior", function(object, ...) standardGeneric("prior"))

setGeneric("refit", function(object, newdata, ...) standardGeneric("refit"))

setGeneric("relabel", function(object, by, ...) standardGeneric("relabel"))


###**********************************************************

setGeneric("info",
function(object, which, ...) standardGeneric("info"))

setMethod("info", signature(object="ANY", which="missing"),
function(object, which, ...)
{
    info(object, which="help")
})
    
infoCheck <- function(object, which, ...)
{
    which %in% info(object, "help")
}



###**********************************************************


