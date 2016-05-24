cv.plsRmodel.default <- function(dataY,dataX,nt=2,limQ2set=.0975,modele="pls", K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12), weights,verbose=TRUE) 
{
if (!(modele %in% c("pls"))) {stop("Use cv.plsRglm to cross-validate PLSRGLRs")}
mf <- match.call(expand.dots = FALSE)
m <- match(c("dataY","dataX","nt","limQ2set","modele","K","NK","grouplist","random","scaleX","scaleY","keepcoeffs","keepfolds","keepdataY","keepMclassed","tol_Xi","weights","verbose"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("PLS_lm_kfoldcv")
cvmodel <- eval(mf, parent.frame())

cvmodel$call <- match.call()
class(cvmodel) <- "cv.plsRmodel"
return(cvmodel)
}
