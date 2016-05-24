#' @title Confidence Intervals (for one population proportion)

#' @description An app to investigate how many times a confidence interval for one population proportion
#'  captures the true population parameter. The true population proportion is plotted as a vertical red line and
#'  the user can visually see how changes to the sample, population proportion, sample size, and confidence level
#'  affect the width of the confidence interval.  Summary information is output to the console to tally the number of times
#'  the computed confidence interval covers the true population mean and how many times
#'  it misses.  
#' 
#' @rdname CIProp
#' @usage CIProp()
#' @return Graphical and numerical output
#' @export
#' @author Rebekah Robinson \email{rebekah_robinson@@georgetowncollege.edu}
#' @note Uses manipulate from RStudio
#' @examples
#' \dontrun{
#' if (require(manipulate)) CIProp()
#' }
CIProp=function(){
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  nprev=NULL
  pprev=NULL
  lprev=NULL
  beginning=TRUE
  tally=c(0,0)
  results=matrix(0,2,2,dimnames=list(c("Count","Proportion"),c("Covers","Misses")))
  manipulate(
    n=slider(1,500,initial=50,label="Number of Trials"),
    p=slider(0,1,initial=0.5,step=0.01,label="Value of p"),
    conf.level=slider(60,99,initial=95,label="Confidence Level"),
    sim.reps=picker("One at a time","20 at a time", "100 at a time",label="Number of Repititions"),  
{
  current.sliders=list(n,p,conf.level)
  prev.sliders=list(nprev,pprev,lprev)
  if(!identical(current.sliders,prev.sliders)){
    tally=c(0,0)
  }
  if(sim.reps=="One at a time"){
    mysamp=sample(c(0,1),n,replace=TRUE,prob=c(1-p,p))
    plot(0:1,0:1,col=0,axes=FALSE,xlab=paste("p=",p),ylab="",main=paste("Confidence Interval for p=",p))
    axis(1,xaxp=c(0,1,10))
    abline(v=p,col="red",lty=2)
    if(!beginning){
      success=sum(mysamp)
      phat=success/n
      conf=conf.level/100
      z.input=conf+((1-conf)/2)
      z=qnorm(z.input)
      se=sqrt((phat*(1-phat))/(n))
      margin=z*se
      ci.left=phat-margin
      ci.right=phat+margin
      segments(x0=ci.left,y0=0,x1=ci.right,y1=0,col="blue",lwd=2)
      points(phat,0,pch=19)
      if(ci.left<p & ci.right>p){
        tally[1]=tally[1]+1
      }
      else{
        tally[2]=tally[2]+1
      } 
    }
    lprev<<-conf.level
    nprev<<-n
    pprev<<-p
    if(!beginning){
      if(sum(tally)==0){
        results[1,]=tally
        results[2,]=c(0,0)
      }
      else{
        results[1,]=tally
        results[2,]=round(c((tally[1]/sum(tally)),(tally[2]/sum(tally))),4)
      }
      print(results)
    }
    beginning=FALSE
  } 
  if(sim.reps=="20 at a time"){
    n.reps=20
    s=replicate(n.reps,sample(c(0,1),n,replace=TRUE,prob=c(1-p,p)))
    success=NULL
    for(i in 1:n.reps){
      success[i]=sum(s[,i])
    }
    phat=success/n
    conf=conf.level/100
    z.input=conf+((1-conf)/2)
    z=qnorm(z.input)
    se=sqrt((phat*(1-phat))/(n))
    margin=z*se
    ci.left=phat-margin
    ci.right=phat+margin
    plot(0:1,0:1,axes=FALSE,col=0,xlab=paste("p=",p),ylab="",main=paste("Confidence Interval for p=",p))
    axis(1,xaxp=c(0,1,10))
    abline(v=p,col="red",lty=2)
    for(i in 1:20){
      segments(x0=ci.left[i],y0=(i-1)*(0.05),x1=ci.right[i],y1=(i-1)*0.05,col="blue",lwd=2)
      points(phat[i],(i-1)*0.05,pch=19) 
      if(ci.left[i]<p & ci.right[i]>p){
        tally[1]=tally[1]+1
      }
      else{
        tally[2]=tally[2]+1
      } 
    }
    lprev<<-conf.level
    nprev<<-n
    if(sum(tally)==0){
      results[1,]=tally
      results[2,]=c(0,0)
    }
    else{
      results[1,]=tally
      results[2,]=round(c((tally[1]/sum(tally)),(tally[2]/sum(tally))),4)
    }
    print(results)
  }
  if(sim.reps=="100 at a time"){
    n.reps=100
    s=replicate(n.reps,sample(c(0,1),n,replace=TRUE,prob=c(1-p,p)))
    success=NULL
    for(i in 1:n.reps){
      success[i]=sum(s[,i])
    }
    phat=success/n
    conf=conf.level/100
    z.input=conf+((1-conf)/2)
    z=qnorm(z.input)
    se=sqrt((phat*(1-phat))/(n))
    margin=z*se
    ci.left=phat-margin
    ci.right=phat+margin
    plot(0:1,0:1,axes=FALSE,col=0,xlab=paste("p=",p),ylab="",main=paste("Confidence Interval for p=",p))
    axis(1,xaxp=c(0,1,10))
    abline(v=p,col="red",lty=2)
    for(i in 1:100){
      segments(x0=ci.left[i],y0=(i-1)*(0.01),x1=ci.right[i],y1=(i-1)*0.01,col="blue",lwd=1)
      points(phat[i],(i-1)*0.01,pch=20)
      if(ci.left[i]<p & ci.right[i]>p){
        tally[1]=tally[1]+1
      }
      else{
        tally[2]=tally[2]+1
      } 
    }
    lprev<<-conf.level
    nprev<<-n
    if(sum(tally)==0){
      results[1,]=tally
      results[2,]=c(0,0)
    }
    else{
      results[1,]=tally
      results[2,]=round(c((tally[1]/sum(tally)),(tally[2]/sum(tally))),4)
    }
    print(results)
  }
})
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("conf.level","sim.reps"))
