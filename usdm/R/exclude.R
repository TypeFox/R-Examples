# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July. 2015
# Version 1.2
# Licence GPL v3

if (!isGeneric("exclude")) {
  setGeneric("exclude", function(x,vif,...)
    standardGeneric("exclude"))
}  

setMethod ('exclude' ,signature(x='RasterStack', vif='VIF'),
           function (x,vif,...) {
             n <- names(x)
             for (i in 1:length(vif@results[,1])) if (!as.character(vif@results[i,1]) %in% n) stop("One or all variables in VIF are not in the Raster object")
             x[[as.character(vif@results[,1])]]
            }
           )


setMethod ('exclude' ,signature(x='RasterBrick', vif='VIF'),
           function (x,vif,...) {
             n <- names(x)
             for (i in 1:length(vif@results[,1])) if (!as.character(vif@results[i,1]) %in% n) stop("One or all variables in VIF are not in the Raster object")
             brick(x[[as.character(vif@results[,1])]])
           }
)

setMethod ('exclude' ,signature(x='data.frame', vif='VIF'),
           function (x,vif, ...) {
             n <- colnames(x)
             for (i in 1:length(vif@results[,1])) if (!as.character(vif@results[i,1]) %in% n) stop("One or all variables in VIF are not in the data.frame object")
             x[,as.character(vif@results[,1])]
           }
)

setMethod ('exclude' ,signature(x='matrix', vif='VIF'),
           function (x,vif, ...) {
             n <- colnames(x)
             for (i in 1:length(vif@results[,1])) if (!as.character(vif@results[i,1]) %in% n) stop("One or all variables in VIF are not in the matrix object")
             x[,as.character(vif@results[,1])]
           }
)

setMethod ('exclude' ,signature(x='RasterStack', vif='missing'),
           function (x,vif,th) {
             n <- names(x)
             if(missing(th)) th <- 10
             vif <- vifstep(x)
             print(vif)
             if (length(vif@excluded) > 0) x[[as.character(vif@results[,1])]]
             else x
           }
)

setMethod ('exclude' ,signature(x='RasterBrick', vif='missing'),
           function (x,vif, th) {
             n <- names(x)
             if(missing(th)) th <- 10
             vif <- vifstep(x,th=th)
             print(vif)
             if (length(vif@excluded) > 0) brick(x[[as.character(vif@results[,1])]])
             else x
           }
)

setMethod ('exclude' ,signature(x='data.frame', vif='missing'),
           function (x,vif, th) {
             n <- colnames(x)
             if(missing(th)) th <- 10
             vif <- vifstep(x,th=th)
             print(vif)
             if (length(vif@excluded) > 0) x[,as.character(vif@results[,1])]
             else x
           }
)

setMethod ('exclude' ,signature(x='matrix', vif='missing'),
           function (x,vif, th) {
             n <- colnames(x)
             if(missing(th)) th <- 10
             vif <- vifstep(x,th=th)
             print(vif)
             if (length(vif@excluded) > 0) x[,as.character(vif@results[,1])]
             else x
           }
)
