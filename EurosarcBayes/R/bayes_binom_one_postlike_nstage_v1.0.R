bayes_binom_one_postlike_nstage=function(reviews,prob.success,prob.failure,eta,zeta,p0,p1,prior.a=1e-6,prior.b=1e-6,round=TRUE,warn=TRUE){

  ln=length(reviews)
  if(length(prob.success)==1 & length(prob.success)!=(ln-1)){
    warning("Length of prob.success does not match number of interim analyses, assumed prob.success is constant at all interim analyses")
    prob.success=rep(prob.success,ln-1)
  }
  if(length(prob.failure)==1 & length(prob.failure)!=(ln-1)){
    warning("Length of prob.failure does not match number of interim analyses, assumed prob.failure is constant at all interim analyses")
    prob.failure=rep(prob.failure,ln-1)
  }

  if(prior.a==0 | prior.b==0){
    stop("Prior parameters must be non-zero")
  }



  # generate success,failure for each interim analysis
  success=rep(NA,ln)
  failure=rep(NA,ln)

  #final endpoint
  success[ln]=which((1-pbeta(p0,prior.a+0:reviews[ln],prior.b+reviews[ln]-0:reviews[ln]))>zeta)[1]-1
  failure[ln]=max(which((pbeta(p1,prior.a+0:reviews[ln],prior.b+reviews[ln]-0:reviews[ln]))>eta))-1

  #interim endpoints
  if(ln>1){
    for(i in 1:(ln-1)){

      probs=1-pbetabinom.ab(success[ln]-0:reviews[i]-1,reviews[ln]-reviews[i], prior.a+0:reviews[i],prior.b+reviews[i]-0:reviews[i])
      failure[i]=max(c(0,which((1-probs)>prob.failure[i])),na.rm =TRUE)-1
      success[i]=min(c(reviews[i]+2,which(probs>prob.success[i])),na.rm =TRUE)-1

      if(success[i]==(reviews[i]+1)){success[i]=NA}
      if(failure[i]==-1){failure[i]=NA}

    }
  }

  if(failure[ln]+1!=success[ln]){
    if(warn==TRUE){
      cat("Design criteria do not meet at final interim analysis, the required design will be altered\n")
      cat("Proposed failure criteria:",failure[ln],"\n")
      cat("Proposed success criteria:",success[ln],"\n")
      cat("success should be 1 greater than failure at endpoint. If it is more than 1 then there is a region of uncertainty. If it is less than 1 then fewer patients are required. Use b.binom.one.onestage to find smallest maximum sample size required. Not all values above this value will have a suitable distribution.")
    }
    failure[ln]=success[ln]-1
  }

  return(properties_binom_one(failure,success,reviews,p0,p1,prior.a,prior.b,round))


}
