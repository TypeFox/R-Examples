binom_one_alpha=function(result.success,result.n,p0,failure,success,n){
  # calculate exact alpha or p-value for a given binomial design for a given result.
  # Adjusts for interim analysis which gives a different result if it is not taken into account.
  # design provided (failure,success,n) should be based on the actual trial and not the proposed design.

  if(length(which(result.n==n))==0){
    stop("The final sample size must match design. Use the actual trial to adjust the design appropriately rather than the planned design")
  }

  prob.lower=0
  prob.upper=0


  ###################################################################
  ln=length(n)
  old=c(1,rep(0,max(n)))
  for(i in 1:ln){
    if(i>1){
      patients=n[i-1]
    } else {
      patients=0
    }

    ne=rep(0,max(n)+1)
    # draw patients from distribution to get to stage i.
    prob=dbinom(0:(n[i]-patients),n[i]-patients,p0)
    # update current trial distribution
    for(j in 1:(patients+1)){
      ne[j:(j+n[i]-patients)]=ne[j:(j+n[i]-patients)] + old[j]*prob
    }

    if(result.n==n[i]){
      if(result.success<=failure[i]){
        if(result.success==0){
          return(1-(prob.lower))
        } else {
          return(1-(prob.lower+sum(ne[1:(result.success)])))
        }

      } else {
        return(prob.upper+sum(ne[(result.success+1):(n[i]+1)]))
      }
    }

    if(!is.null(failure)){
      if(!is.na(failure[i])){
        prob.lower=prob.lower+sum(ne[1:(failure[i]+1)])
        ne[1:(failure[i]+1)]=0
      }}

    if(!is.na(success[i])){
      prob.upper=prob.upper+sum(ne[(success[i]+1):(n[i]+1)])
      ne[(success[i]+1):(n[i]+1)]=0
    }
    old=ne
  }
}
