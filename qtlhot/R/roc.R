CausalRocCurve <- function(outA, outB, outC, outD, outE, 
  method, alpha=0.05)
{
  ntests <- nrow(outA[[1]])
  ctA <- counts(outA, alpha, method)
  ctB <- counts(outB, alpha, method)
  ctC <- counts(outC, alpha, method)
  ctD <- counts(outD, alpha, method)
  ctE <- counts(outE, alpha, method)
  if(method != "cit"){
    tp <- ctA[1,1]+ctB[1,1]
    fn <- ctC[1,1]+ctD[1,1]+ctE[1,1]
    tn <- ctD[1,3]+ctC[1,4]+ctE[1,4]
    fp <- ctA[1,2:4]+ctB[1,2:4]
  }
  if(method == "cit"){
    tp <- ctA[1,1]+ctB[1,1]
    fn <- ctC[1,1]+ctD[1,1]+ctE[1,1]
    tn <- ctD[1,3]+ctC[1,3]+ctE[1,3]
    fp <- ctA[1,2:3]+ctB[1,2:3]
  }
  tpr <- tp/(tp+fn)
  fpr <- fp/(fp+tn)
  data.frame(tpr, fpr, tp, fp, tn, fn)
}
#############################################################################
causal.ROC.curve <- function(...) CausalRocCurve(...)
#############################################################################
GetRocMatrix <- function(outA, outB, outC, outD, outE, alpha,
                           condense.labels = TRUE, verbose = FALSE)
{
  n <- length(alpha)
  
  mymethods <- if(condense.labels)
    c("aic", "j.aic", "p.aic", "np.aic",
      "bic", "j.bic", "p.bic", "np.bic", "cit")
  else
    c("aic","par.cmst.joint.aic","par.cmst.aic","non.par.cmst.aic",
      "bic","par.cmst.joint.bic","par.cmst.bic","non.par.cmst.bic","cit")
  
  Tpr <- matrix(NA,9,n, dimnames=list(mymethods, as.character(alpha)))
  Fpr <- matrix(NA,9,n, dimnames=list(mymethods, as.character(alpha)))
  for(k in 1:n){
    for(j in 1:9){
      aux <- causal.ROC.curve(outA, outB, outC, outD, outE, method=mymethods[j], alpha[k])
      Tpr[j,k] <- aux[1,1]
      Fpr[j,k] <- aux[1,2]
    }
    if(verbose)
      print(k)
  }
  list(Tpr=Tpr, Fpr=Fpr)
}
#############################################################################
get.Roc.Matrix <- function(...) GetRocMatrix(..., condense.labels = TRUE)
