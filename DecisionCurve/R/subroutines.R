

#calculate net benefit for a fitted risk or by fitting a risk model if formula is specified
calculate.nb <- function(y, d, rH, formula, data, family, formula.ind, casecontrol.rho){
  #predictor y,  estimated risk
  #disease indicator d
  #vector of high risk thresholds rH

  if(formula.ind){

    #cohort data, no correction for population level outcome rate
    if(is.null(casecontrol.rho)){
       myglm <- do.call(glm, list("formula" = formula, "data" = data, "family" = family ))
       y <- fitted(myglm)
    }else{
      #offset by the relative observed outcome prevalence and the provided population rho
      obs.rho = mean(d)
      offset = - log((casecontrol.rho)/ (1-casecontrol.rho)) + log((obs.rho)/(1-obs.rho))
      myglm <- do.call(glm, list("formula" = formula, "data" = data, "family" = family, "offset" = rep(offset, nrow(data)) ))
      y = predict(myglm, type = "link") - offset
      y <- exp(y)/(1+exp(y))
    }

  }

  N = length(y)

  #denominator for tpf
  tpf.den <- sum(d ==1)

  #true positive fraction Pr(risk(X) > rH | D = 1)
  tpf <- sum.I(rH, "<", y[d==1])/tpf.den

  #false positive fraction Pr(risk(X) > rH | D = 0 )
  fpf <- sum.I(rH, "<", y[d==0])/ (N - tpf.den)

  #disease prevalence
if(is.null(casecontrol.rho)){
    rho = mean(d ==1)
    prob.high.risk <- sum.I(rH, "<", y)/length(y)

  }else{
    rho = casecontrol.rho
    prob.high.risk <- rho*sum.I(rH, "<", y[d==1])/length(y[d==1]) + (1-rho)*sum.I(rH, "<", y[d==0])/length(y[d==0])
  }


  #net benefit
  nb = tpf*rho - (rH/(1-rH))*(1-rho)*fpf

  #standardized net benefit
  snb = nb/rho

  #detection probability (for impact plots)
  dp <- tpf*rho

  out  = data.frame("threshold" = rH,
                    "FPR" = fpf , "TPR" = tpf,
                    "NB" = nb, "sNB" = snb,
                    "rho" = rho, "prob.high.risk" = prob.high.risk,
                    "DP" = dp)

 # AUC   = sum(sort(tpf, decreasing = FALSE)*(sort(fpf, decreasing = FALSE)-c(sort(fpf, decreasing = FALSE)[-1],0)))


}



add.ci.columns <- function(x){
  n.out = nrow(x)

  x$FPR_lower <- NA; x$FPR_upper <- NA
  x$TPR_lower <-NA; x$TPR_upper <- NA
  x$NB_lower <-NA; x$NB_upper <- NA
  x$sNB_lower <-NA; x$sNB_upper <- NA
  x$rho_lower <- NA; x$rho_upper <- NA
  x$prob.high.risk_lower = NA; x$prob.high.risk_upper =NA
  x$DP_lower = NA; x$DP_upper =NA
  x
}


VTM <- function(vc, dm)
{
  matrix(vc, ncol = length(vc), nrow = dm, byrow = T)
}


#this function is magic!
#for each element of yy it gives the count of how many Yi fall <, >, <=, or >= to it
#output is a vector of length yy.
sum.I<-function(yy,FUN,Yi,Vi=NULL){

  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}

  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')

  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos

  if (!is.null(Vi)) {

    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)

    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)

    return(rbind(0,Vi)[pos+1,])

  } else return(pos)
}


costbenefit_to_threshold <- function(CB){
  CB/(1+CB)
}

threshold_to_costbenefit <- function(rh){
  rh/(1-rh)
}
