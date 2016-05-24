plsRmodel.default <- function(dataY,dataX,nt=2,limQ2set=.0975,dataPredictY=dataX,modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12),weights,sparse=FALSE,sparseStop=TRUE,naive=FALSE,verbose=TRUE)
{
if(missing(modele)){modele="pls"}
if (!(modele %in% c("pls"))) {stop("Use plsRglm for applying PLSR to glms")}
mf <- match.call(expand.dots = FALSE)
m <- match(c("dataY","dataX","nt","limQ2set","dataPredictY","modele","family","typeVC","EstimXNA","scaleX","scaleY","pvals.expli","alpha.pvals.expli","MClassed","tol_Xi","weights","sparse","sparseStop","naive","verbose"), names(mf), 0L)
mf <- mf[c(1L, m)]
if(is.null(mf$modele)){mf$modele<-"pls"}
mf[[1L]] <- as.name("PLS_lm")
estmodel <- eval(mf, parent.frame())

  estmodel$call <- match.call()
  class(estmodel) <- "plsRmodel"
  estmodel
}
