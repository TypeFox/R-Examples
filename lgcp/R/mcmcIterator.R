###' iterator for MCMC loops
###'
###' control an MCMC loop with this iterator
###'
###' @param N number of iterations
###' @param burnin length of burn-in
###' @param thin frequency of thinning
###' @param trim whether to cut off iterations after the last retained iteration
###' @param progressor a function that returns a progress object
###' @export

mcmcLoop <- function(N,burnin,thin, trim=TRUE, progressor=mcmcProgressPrint){

  waste = (N-burnin) %% thin
  if(waste !=0){
    if(trim){
      N <- N - waste 
      warning(waste," iterations trimmed off")
    }else{
      warning(waste," iterations will be wasted after last retained iteration.")
    }
  }
  
  ii = 0
  nextEl <- function(){
    ## when nextEl is called in a while() loop, ii is the number of iterations completed so far.
    if(ii==0){
      p$start()
    }
    p$update()
    ii <<- ii + 1
    ## ii is now the number of the iteration we are about to do
    fin = !hasNx()
    if(fin){
      p$end()
    }
    return(!fin)
  }
  itrtr <- function(){
    return(ii)
  }
  isb <- function(){
    if(ii<=burnin){
      return(TRUE)
    }
    return(FALSE)
  }
  isr <- function(){
    if(isb()){
      return(FALSE)
    }
    if((ii-burnin) %% thin == 0){
      return(TRUE)
    }
    return(FALSE)
  }
  hasNx <- function(){
    if(ii <= N){
      return(TRUE)
    }
    return(FALSE)
  }

  rst <- function(){
    ii <<- 0
  }
  
  obj = list(nextElem=nextEl,
    hasNext=hasNx,
    is.retain=isr,
    is.burnin=isb,
    iteration=itrtr,
    restart=rst,
    N=N,
    burnin=burnin,
    thin=thin
    )
  class(obj) <- c("mcmc","abstractiter","iter")

  p = progressor(obj)
  obj$progress = p
  
  obj
}

###' print.mcmc function
###'
###' print method
###' print an mcmc iterator's details
###'
###' @method print mcmc
###' @param x a mcmc iterator
###' @param ... other args
###' @export
print.mcmc <- function(x,...){
  cat("mcmc loop controller\n\n")
  cat("   iterations:",x$N,"\n")
  cat("      burn-in:",x$burnin,"\n")
  cat("     thinning: 1/",x$thin,"\n",sep="")
  cat("      current:",iteration(x),"\n")
}

###' generic hasNext method
###'
###' test if an iterator has any more values to go
###' @param obj an iterator
###' @export
hasNext <- function(obj) {
  UseMethod("hasNext")
}

###' hasNext.iter function
###' 
###' method for iter objects
###' test if an iterator has any more values to go
###'
###' @method hasNext iter
###' @param obj an iterator
###' @export
hasNext.iter <- function(obj){
  obj$hasNext()
}

###' do we retain this iteration?
###'
###' if this mcmc iteration is one not thinned out, this is true
###' @param obj an mcmc iterator
###' @return TRUE or FALSE
###' @export
is.retain <- function(obj){
  obj$is.retain()
}

###' is this a burn-in iteration?
###'
###' if this mcmc iteration is in the burn-in period, return TRUE
###' @param obj an mcmc iterator
###' @return TRUE or FALSE
###' @export
is.burnin <- function(obj){
  obj$is.burnin()
}

###' iteration number
###'
###' within a loop, this is the iteration number we are currently doing.
###'
###' get the iteration number
###' @param obj an mcmc iterator
###' @return integer iteration number, starting from 1.
###' @export
iteration <- function(obj){
  obj$iteration()
}

###' reset iterator
###'
###' call this to reset an iterator's state to the initial
###' @param obj an mcmc iterator
###' @export
resetLoop <- function(obj){
  obj$restart()
}

###' summary.mcmc function
###' 
###' summary of an mcmc iterator
###' print out values of an iterator and reset it. DONT call this
###' in a loop that uses this iterator - it will reset it. And break.
###'
###' @method summary mcmc
###' @param object an mcmc iterator
###' @param ... other args
###' @export
summary.mcmc <- function(object,...){
  resetLoop(object)
   while(nextStep(object)){
    cat("Iteration ",iteration(object)," burnin: ",is.burnin(object)," retain: ",is.retain(object),"\n")
  }
   resetLoop(object)
 }

##' next step of an MCMC chain
##'
##' just a wrapper for nextElem really.
##'
##' @param object an mcmc loop object
##' @export
nextStep <- function(object){nextElem(object)}

##' loop over an iterator
##'
##' useful for testing progress bars
##' @param object an mcmc iterator
##' @param sleep pause between iterations in seconds
##' @export
loop.mcmc <- function(object,sleep=1){
  while(nextStep(object)){
    Sys.sleep(sleep)
  }
}
