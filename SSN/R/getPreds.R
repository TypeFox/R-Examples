################################################################################
# getPreds.R
#   The purpose of this function is to get the predictions and standard
#   errors from a glmssn.predict object. For observed sites, the predictions represent
#   the cross-validation predictions and associated standard errors.
################################################################################

getPreds <- function(x,pred.type="cv") {

    if(pred.type == "cv") {
        pred.tmp <- CrossValidationSSN(x)

        pred.mat <- matrix(data = NA, nrow = nrow(pred.tmp), ncol = 3)
        colnames(pred.mat) <- c("pid", "cv.pred", "cv.se")

        pred.mat[,"pid"]<- x$ssn.object@obspoints@SSNPoints[[1]]@point.data[,"pid"]
        pred.mat[,"cv.pred"]<- pred.tmp[,"cv.pred"]
        pred.mat[,"cv.se"]<- pred.tmp[,"cv.se"]

    } else if(pred.type == "pred")  {
        if (length(x$ssn.object@predpoints@ID) == 0) {
            stop("no prediction points exist in glmssn.predict object") }
        else if (length(x$ssn.object@predpoints@ID) == 1){
            predpointsID <- x$arg$predpoints
            slot.i <- 1 }
        else {
            for (i in 1:length(x$ssn.object@predpoints@ID)){
                if (x$args$predpointsID == x$ssn.object@predpoints@ID[[i]]) {
                    predpointsID <- x$ssn.object@predpoints@ID[[i]]
                    slot.i <- i
                }
            }}

        pred.name <- x$args$zcol
        predSE.name <- paste(pred.name, ".predSE", sep = "")
        predSE.name2<- paste(pred.name, "PredSE", sep = "")
        pred.mat <- cbind(x$ssn.object@predpoints@SSNPoints[[slot.i]]@point.data[,"pid"],
                          x$ssn.object@predpoints@SSNPoints[[slot.i]]@point.data[,pred.name],
                          x$ssn.object@predpoints@SSNPoints[[slot.i]]@point.data[,predSE.name])
        colnames(pred.mat) <- c("pid", pred.name, predSE.name2)

    } else {
        stop("invalid prediction type")
    }

    return(pred.mat)
}


