binom_one_power=function(p,failure,success,n){
  # power function for binomial single arm trial for a provided design.

  probsuccess=rep(0,length(p))
  for(k in 1:length(p)){
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
      prob=dbinom(0:(n[i]-patients),n[i]-patients,p[k])
      # update current trial distribution
      for(j in 1:(patients+1)){
        ne[j:(j+n[i]-patients)]=ne[j:(j+n[i]-patients)] + old[j]*prob
      }
      if(!is.null(failure)){
        if(!is.na(failure[i])){
          ne[1:(failure[i]+1)]=0
        }}

      if(!is.na(success[i])){
        probsuccess[k]=probsuccess[k]+sum(ne[(success[i]+1):(n[i]+1)])
        ne[(success[i]+1):(n[i]+1)]=0
      }
      old=ne
    }
  }
  return(probsuccess)

}

# print the power curve for a given trial design
plot_binom_one_power=function(failure,success,n,ndivisions=1000,xlim=c(0,1),xaxs="i",yaxs="i",ylim=c(0,1.1),main="Power curve for a single arm binomial trial design",xlab="Probability of successful treatment",ylab="Probability of successful trial",p=NULL,alpha=NULL,power=NULL,col.error="red",...){
  p0=(xlim[2]-xlim[1])*(0:ndivisions)/ndivisions-xlim[1]
  pow=binom_one_power(p0,failure,success,n)
  plot(p0,pow,type="l",xlim=xlim,xaxs=xaxs,yaxs=yaxs,ylim=ylim,xlab=xlab,ylab=ylab,main=main,...)
  if(!is.null(p)){
    for(i in 1:length(p)){
      segments(p[i],0,p[i],binom_one_power(p[i],failure,success,n))
    }
  }
  if(!is.null(alpha)){
    abline(h=alpha,lty=2,col=col.error)
  }
  if(!is.null(power)){
    abline(h=power,lty=2,col=col.error)
  }
}
