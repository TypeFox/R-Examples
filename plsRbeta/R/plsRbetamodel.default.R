plsRbetamodel.default <- function(dataY,dataX,nt=2,limQ2set=.0975,dataPredictY=dataX,modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12),weights,method,sparse=FALSE,sparseStop=TRUE,naive=FALSE,link=NULL,link.phi=NULL,type="ML")
{
mf <- match.call(expand.dots = FALSE)
m <- match(c("dataY","dataX","nt","limQ2set","dataPredictY","modele","family","typeVC","EstimXNA","scaleX","scaleY","pvals.expli","alpha.pvals.expli","MClassed","tol_Xi","weights","method","sparse","sparseStop","naive","link","link.phi","type"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("PLS_beta")
estmodel <- eval(mf, parent.frame())

  if (typeVC!="none") {stop("Use plsRbeta_kfoldcv for applying kfold cross validation to beta regression models")}
  class(estmodel) <- "plsRbetamodel"
  estmodel$call <- match.call()
  estmodel
}
