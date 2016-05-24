## THIS FILE NEEDS THOROUGH REVISION DUE TO THE INTRODUCTION OF
## FORWARD/BACKWARD WARPING !!!

## Paril 1 :-), 2011 (RW)
## First working version of predict.ptw.

## Eventually this function should also allow the prediction with a
## different set of warping coefficients, but we'll leave that for now.
## The way to do it is to create a ptw object with the try=TRUE option
## and THEN do the prediction.

predict.ptw <- function(object, newdata, 
                        what = c("response", "time"),
                        RTref = NULL, ...)
{
  what <- match.arg(what)

  switch(what,
         response = {
           if (missing(newdata)) return(object$warped.sample)

           if (!is.matrix(newdata))
             newdata <- matrix(newdata, nrow = 1)
           
           if (object$warp.type == "individual" &
               nrow(newdata) > 1 &
               nrow(newdata) != nrow(object$warp.fun))
             stop("Incorrect number of rows in newdata")

           ## if necessary, replicate the one warping function so that
           ## every row in newdata has one warping function
           if (object$warp.type == "individual") {
             WF <- object$warp.fun
           } else {
             WF <- matrix(object$warp.fun,
                          nrow(object$sample),
                          ncol(object$warp.fun), byrow = TRUE)
             
           }

           ## do the warping
           t(sapply(1:nrow(newdata),
                    function(i) ##interpol(WF[i,], newdata[i,])))
                      approx(newdata[i,], NULL, WF[i,])$y))
         },
         time = {
           correctedTime <- 
             -sweep(object$warp.fun, 2, 2*(1:ncol(object$ref)), FUN = "-")
           if (is.null(RTref)){
             if (is.null(colnames(object$ref))) {
               RTref <- 1:ncol(object$ref)
             } else {
               RTref <- as.numeric(colnames(object$ref))
             }
           }

           if (missing(newdata)) {
             newdata <- RTref
             newdataIndices <- 1:length(RTref)
           } else {
             ## the round statement is necessary to explicitly convert
             ## newdataIndices to an integer value: otherwise the
             ## interpol function later will give back an incorrect
             ## value. Why???             
             newdataIndices <-
                 round((newdata - min(RTref)) * (length(RTref) - 1) /
                       diff(range(RTref)) + 1)
           }

           t(sapply(1:nrow(correctedTime),
                    function(i)
                      approx(RTref, NULL, correctedTime[i, newdataIndices])$y))
##                    interpol(correctedTime[i, newdataIndices],
##                             RTref)))
         })
}
