

bayes_binom_one_postprob_nstage=function(reviews,eta,zeta,p0,p1,prior.a=0,prior.b=0,h0=p0,h1=p1,round=TRUE,warn=TRUE){

  ln=length(reviews)
  if(length(eta)==1 & length(eta)!=ln){
    warning("Length of eta does not match number of interim analyses, assumed eta is constant at all interim analyses")
    eta=rep(eta,ln)
  }
  if(length(zeta)==1 & length(zeta)!=ln){
    warning("Length of zeta does not match number of interim analyses, assumed zeta is constant at all interim analyses")
    zeta=rep(zeta,ln)
  }
  # vectorised p0 and p1 so need not be constant throughout trial. note Bayesian properties will be based on h0, h1 if provided or last p0 p1 if not provided.
  if(length(p0)!=ln){
    p0=rep(p0,ln)
  }
  if(length(p1)!=ln){
    p1=rep(p1,ln)
  }
  # h0 and h1 are just for frequentist properties so may be different from p0 and p1 for the design
  if(is.null(h0)){
    h0=p0[ln]
  }
  if(is.null(h1)){
    h1=p1[ln]
  }

  # generate success,failure for each interim analysis
  success=rep(NA,ln)
  failure=rep(NA,ln)
  for(i in 1:ln){
    success[i]=min(c(reviews[i]+2,which((1-pbeta(p0[i],prior.a+0:reviews[i],prior.b+reviews[i]-0:reviews[i]))>zeta[i])),na.rm =TRUE)-1
    if(success[i]==(reviews[i]+1)){success[i]=NA}
    failure[i]=max(c(0,which((pbeta(p1[i],prior.a+0:reviews[i],prior.b+reviews[i]-0:reviews[i]))>eta[i])),na.rm =TRUE)-1
    if(failure[i]==-1){failure[i]=NA}
  }

  if(failure[ln]+1!=success[ln]){
    if(warn==TRUE){
      cat("Design criteria do not meet at final interim analysis, the required design will be altered\n")
      cat("Proposed failure criteria:",failure[ln],"\n")
      cat("Proposed success criteria:",success[ln],"\n")
      cat("success should be 1 greater than failure at endpoint. If it is more than 1 then there is a region of uncertainty. If it is less than 1 then fewer patients are required. Use b.binom.single.onestage to find smallest maximum sample size required. Not all values above this value will have a suitable distribution.")
    }
    failure[ln]=success[ln]-1
  }

  return(properties_binom_one(failure,success,reviews,p0=h0,p1=h1,prior.a,prior.b,round))

}
