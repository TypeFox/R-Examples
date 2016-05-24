#' @title Type I and Type II Errors

#' @description An app to explore the concepts of Type I and Type II errors, and the concept of
#' power.  We take samples from a population that is imagined to be normal, and perform the t-procedures
#' for one mean.  The Null Hypotheis is H0:  mu=170.  A slider allows us to vary the true mean mu.
#' 
#' @rdname Type12Errors
#' @usage Type12Errors() 
#' @return Graphical and numerical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @note Uses \code{manipulate}.
#' @examples
#' \dontrun{
#' if (require(manipulate)) Type12Errors()
#' }
Type12Errors <-
function(){
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  #imaginary population is normal(170,20)
  nprev <- NULL
  lprev <- NULL
  tprev <- NULL
  tally <- c(0,0)
  beginning <- TRUE
  ymax <- 0.04
  manipulate( 
    true=slider(165,175,initial=170,step=1,label="True Population Mean mu"),
    sig.level=slider(0.01,0.10,step=0.01,initial=0.05,label="Level of Significance alpha"),
    n=slider(2,50,initial=10,label="Sample Size n"),
    sim.reps=picker("One at a time","100 at a time","1000 at a time",label="Number of Repetitions"),
{
  conf <- (1-sig.level)/100
  current.sliders <- list(true,sig.level,n)
  prev.sliders <- list(tprev,lprev,nprev)
  if(!identical(current.sliders,prev.sliders)){
    tally <-c(0,0)
  }
  
  if(sim.reps=="One at a time"){
    samp=rnorm(n,mean=true,sd=20)
    color <- ifelse(beginning,"transparent","lightblue")
    border <- ifelse(beginning,"transparent","black")
    hist(samp,freq=FALSE,
         main="Histogram of Sample Male Weights\n(w/Density Curve for Population)",
          xlab="Weights",col=color,border=border,
         xlim=c(90,250),ylim=c(0,ymax))
    abline(v=true,col="red",lwd=2)
    
    curve(dnorm(x,mean=true,sd=20),from=110,to=230,lwd=2,col="red",type="l",add=TRUE)
    
    axis(1,at=170,labels=expression(mu[0]))
    text(expression(mu),x=true+4,y=ymax)
    
    if(!beginning) {
      xbar=mean(samp)
      t.input=sig.level/2
      t=qt(t.input,df=n-1)
      se=sd(samp)/sqrt(n)
      margin=t*se
      ci=c(xbar-margin,xbar+margin)
      segments(x0=ci[1],y0=0,x1=ci[2],y1=0,col="yellow",lwd=2)
      points(xbar,0,col="blue",pch=20)
      results <- t.test(samp,mu=170,alternative="two.sided")
      null.rejected <- (results$p.value<sig.level)   
      tally <- tally+c(null.rejected,!null.rejected)
    }
    
    
       
    lprev<<-sig.level
    nprev<<-n
    tprev <<- true
    if (!beginning){
    null.rejected <- ifelse(null.rejected,"yes","no")
    output <- data.frame(null.rejected,rejected.so.far=tally[1],
                         total.samples=sum(tally),
                         proportion.rejected=round(tally[1]/sum(tally),3))
     print(output) }
    beginning <<- FALSE
  }#end one at a time
  
  
  
  if(sim.reps=="100 at a time"){
    reps <- 100
    null.rejected <- logical(reps)
    for(i in 1:reps){
      samp=rnorm(n,mean=true,sd=20)
      results <- t.test(samp,mu=170,alternative="two.sided")
      null.rejected[i] <- (results$p.value<sig.level)
    }
    rejections <- sum(null.rejected)
    
    tally <- tally+c(rejections,reps-rejections)   
    
    lprev<<-sig.level
    nprev<<-n
    tprev <<- true
    
    output <- data.frame(rejected.so.far=tally[1],
                            total.samples=sum(tally),
                            proportion.rejected=round(tally[1]/sum(tally),3))
    print(output)
    
  }#end 100 at a time
  
  
  if(sim.reps=="1000 at a time"){
    reps <- 1000
    null.rejected <- logical(reps)
    for(i in 1:reps){
      samp=rnorm(n,mean=true,sd=20)
      results <- t.test(samp,mu=170,alternative="two.sided")
      null.rejected[i] <- (results$p.value<sig.level)
    }
    rejections <- sum(null.rejected)
    
    tally <- tally+c(rejections,reps-rejections)   
    
    lprev<<-sig.level
    nprev<<-n
    tprev <<- true
    
    output <- data.frame(rejected.so.far=tally[1],
                         total.samples=sum(tally),
                         proportion.rejected=round(tally[1]/sum(tally),3))
    print(output)
    
  }#end 1000 at a time
  
  
}#end manip body   
  )#end manip
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("picker", "sig.level","true","sim.reps"))
