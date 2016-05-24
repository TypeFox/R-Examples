conformity.Check <-
function(ordPmat, CorrMat){
if (ncol(ordPmat)!=ncol(CorrMat) ){stop("dimension of CorrMat does not conform to ordPmat")}
if (ncol(CorrMat)!=nrow(CorrMat) ){stop("The CorrMat must be a square matrix")}
}
