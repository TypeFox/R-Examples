plsRmodel.formula <- function(formula,data,nt=2,limQ2set=.0975,dataPredictY,modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12), weights,subset,contrasts=NULL,sparse=FALSE,sparseStop=TRUE,naive=FALSE,verbose=TRUE)
{
  if(missing(modele)){modele="pls"}
  if (!(modele %in% c("pls"))) {stop("Use plsRglm for applying PLSR to glms")}
if (missing(data)) {data <- environment(formula)}
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula","data","nt","limQ2set","dataPredictY","modele","family","typeVC","EstimXNA","scaleX","scaleY","pvals.expli","alpha.pvals.expli","MClassed","tol_Xi","weights","subset","contrasts","sparse","sparseStop","naive","verbose"), names(mf), 0L)
mf <- mf[c(1L, m)]
if(is.null(mf$modele)){mf$modele<-"pls"}
mf[[1L]] <- as.name("PLS_lm_formula")
estmodel <- eval(mf, parent.frame())
  class(estmodel) <- "plsRmodel"
  estmodel$call <- match.call()
  estmodel
}
