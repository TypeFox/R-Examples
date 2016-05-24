cv.plsRglmmodel.default <- function(dataY,dataX,nt=2,limQ2set=.0975,modele="pls", family=NULL, K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12), weights, method, verbose=TRUE)
  {
  mf <- match.call(expand.dots = FALSE)
m <- match(c("dataY","dataX","nt","limQ2set","modele", "family", "K", "NK", "grouplist", "random", "scaleX", "scaleY", "keepcoeffs", "keepfolds", "keepdataY", "keepMclassed", "tol_Xi", "weights", "method", "verbose"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("PLS_glm_kfoldcv")
cvmodel <- eval(mf, parent.frame())
  class(cvmodel) <- "cv.plsRglmmodel"
  cvmodel$call <- match.call()
  cvmodel
}
