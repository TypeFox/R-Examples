
###' null progress monitor
###'
###' a progress monitor that does nothing
###' @param mcmcloop an mcmc loop iterator
###' @return a progress monitor
###' @export

mcmcProgressNone <- function(mcmcloop){

  force(mcmcloop)
  
  start = function(){
  }
  update=function(){
  }
  end=function(){
  }

  obj = list(start=start,update=update,end=end)
  class(obj) <- "mcmcprogress"
  return(obj)
  
}

###' printing progress monitor
###'
###' a progress monitor that prints each iteration
###' @param mcmcloop an mcmc loop iterator
###' @return a progress monitor
###' @export

mcmcProgressPrint <- function(mcmcloop){
  force(mcmcloop)

  start = function(){
    cat("Initialising\n")
  }
  update = function(){

    cat("Iteration ",iteration(mcmcloop))
    if(is.burnin(mcmcloop)){
      cat(" - burn-in ")
    }
    if(is.retain(mcmcloop)){
      cat(" - retained")
    }
    cat("\n")
  }
      
  end = function(){
    cat("Finished\n")
  }
  obj = list(start=start,update=update,end=end)
  class(obj) <- "mcmcprogress"
  return(obj)
}

###' text bar progress monitor
###'
###' a progress monitor that uses a text progress bar
###' @param mcmcloop an mcmc loop iterator
###' @return a progress monitor
###' @export

mcmcProgressTextBar <- function(mcmcloop){
  if(options()$width<60){
    warning("progress bar will look bad in narrow window")
  }
  force(mcmcloop)
  e=environment()
  start = function(){
    w = options()$width
    options(width=max(c(10,w-25)))
    assign("burnbar",txtProgressBar2(1,mcmcloop$burnin,style=3,label="Burn-in"),envir=e)
    assign("iterbar",txtProgressBar2(1,mcmcloop$N-mcmcloop$burnin,style=3,label="Run"),envir=e)
    options(width=w)
  }
  update = function(){
    if(is.burnin(mcmcloop)){
      setTxtProgressBar2(with(e,burnbar),iteration(mcmcloop)) 
    }else{
      setTxtProgressBar2(with(e,iterbar),iteration(mcmcloop)-mcmcloop$burnin,
                         label=paste("Run [",iteration(mcmcloop),"/",mcmcloop$N,"]",sep=""))
    }
  }
      
  end = function(){
    close(with(e,burnbar))
    close(with(e,iterbar))
  }
  obj = list(start=start,update=update,end=end)
  class(obj) <- "mcmcprogress"
  return(obj)
}


### This function not used as it requires a dependency on tcltk
###
### graphical progress monitor
###
### a progress monitor that uses tcltk dialogs
### @param mcmcloop an mcmc loop iterator
### @return a progress monitor
### @export

#mcmcProgressTk <- function(mcmcloop){
#  force(mcmcloop)
#  e=environment()
#  start = function(){
#    assign("bar",tkProgressBar(title="Burn-in",label="Burn-in phase",min=1,max=mcmcloop$N),envir=e)
#  }
#  update = function(){
#    if(is.burnin(mcmcloop)){
#      title="Burn-in"
#      label=paste("Burn-in ",iteration(mcmcloop),"/",mcmcloop$burnin,sep="")
#    }else{
#      title="Running"
#      label = paste(iteration(mcmcloop),"/",mcmcloop$N,sep="")
#      if(is.retain(mcmcloop)){
#        label=paste(label," (retained)")
#      }
#    }
#    setTkProgressBar(with(e,bar),iteration(mcmcloop),label=label,title=title)
#  }
#  end = function(){
#    close(with(e,bar))
#  }
#  obj = list(start=start,update=update,end=end)
#  class(obj) <- "mcmcprogress"
#  return(obj)
#}
