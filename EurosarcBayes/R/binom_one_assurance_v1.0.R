


binom_one_assurance=function(failure,success,n,ass.dist,type="continuous",lower=0,upper=1,...){
  # Calculates bayesian assurance for a specific design and prior expectation of treatment (ass.dist).
  # Effectively a wrapper of binomial.power which integrates the power of a design weighting p according to ass.dist
  # Can take either a discrete distribution or probability density function

  if(type=="discrete"){

    prob.drug=ass.dist[,1]
    if(min(prob.drug)<0 | max(prob.drug)>1){
      stop("Probability of the drug must be between 0 and 1")
    }
    weight=ass.dist[,2]
    if((1-sum(weight))>10^-6){
      stop("Probability does not sum to 1")
    }

    return(sum(weight*binom_one_power(prob.drug,failure,success,n)))

  } else if(type=="continuous"){

    if((1-integrate(ass.dist,lower=lower,upper=upper,...)$value)>10^-6){
      warning("Probability for assurance distribution does not sum to 1")
    }


    intergrand=function(p,n,success,failure,...){
      return(ass.dist(p,...)*binom_one_power(p,failure,success,n))
    }

    prob.success=(integrate(intergrand,failure=failure,success=success,n=n,lower=lower,upper=upper, ...))
    #print(prob.success)
    return(prob.success$value)
  } else {
    stop("Expecting type = 'continuous' or 'discrete'")
  }

}

plot_binomassurance=function(failure,success,n,ass.dist,type="continuous",ndivisions=1000,xlim=c(0,1),xaxs="i",yaxs="i",ylim=NULL,main="Assurance distribution",col="red",col.fill="green",lwd=2,xlab="Probability of successful treatment",ylab="Prior assurance probability",...){

  if(type=="continuous"){
    po=(xlim[2]-xlim[1])*(0:ndivisions)/ndivisions-xlim[1]
    pow=binom_one_power(po,failure,success,n)
    plot.ass.dist=ass.dist(po)

    if(is.null(ylim)){
      ylim=range(pretty(plot.ass.dist))
      ylim[2]=1.1*ylim[2]
    }

    plot(po,plot.ass.dist,type="l",xlim=xlim,xaxs=xaxs,yaxs=yaxs,ylim=ylim,xlab=xlab,ylab=ylab,main=main,...)
    polygon(c(0,po,1),c(0,plot.ass.dist,0),col=col,border=col)
    polygon(c(0,po,1),c(0,pow*plot.ass.dist,0),col=col.fill,border=col.fill)
    lines(po,plot.ass.dist)
    box()

  } else if(type=="discrete"){

    po=ass.dist[,1]
    pow=binom_one_power(po,failure,success,n)
    plot.ass.dist=ass.dist[,2]
    if(is.null(ylim)){
      ylim=c(0,range(pretty(plot.ass.dist))[2])
      ylim[2]=1.1*ylim[2]
    }

    plot(po,plot.ass.dist,col="white",type="l",xlim=xlim,xaxs=xaxs,yaxs=yaxs,ylim=ylim,main=main,xlab=xlab,ylab=ylab,lwd=lwd,...)
    for(i in 1:length(po)){
      segments(po[i],0,po[i],plot.ass.dist[i],col=col,lwd=lwd)
      segments(po[i],0,po[i],pow[i]*plot.ass.dist[i],col=col.fill,lwd=lwd)
    }
    box()
  } else {
    stop("Expecting type = 'distribution' or 'discrete'")
  }

}




