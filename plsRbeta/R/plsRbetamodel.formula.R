plsRbetamodel.formula <- function(formula,data=NULL,nt=2,limQ2set=.0975,dataPredictY,modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12),weights,subset,start=NULL,etastart,mustart,offset,method="glm.fit",control= list(),contrasts=NULL,sparse=FALSE,sparseStop=TRUE,naive=FALSE,link=NULL,link.phi=NULL,type="ML")
{
if (typeVC!="none") {stop("Use plsRbeta_kfoldcv for applying kfold cross validation to beta regression models")}
if (missing(data)) {data <- environment(formula)}
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula","data","nt","limQ2set","dataPredictY","modele","family","typeVC","EstimXNA","scaleX","scaleY","pvals.expli","alpha.pvals.expli","MClassed","tol_Xi","weights","subset","start","etastart","mustart","offset","method","control","contrasts","sparse","sparseStop","naive","link","link.phi","type"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("PLS_beta_formula")
estmodel <- eval(mf, parent.frame())
  class(estmodel) <- "plsRbetamodel"
  estmodel$call <- match.call()
  estmodel
}
