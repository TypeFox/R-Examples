#' Set parameter bounds
#'
#' Set values of lower and upper bounds and check lengths of any user-specified
#' values
#'
#' @param lowerbounds vector of lower bounds
#' @param upperbounds vector of upper bounds
#' @param initialvalues vector of initial parameter estimates
#' @param ddfobj distance detection function object
#' @return \item{lower}{vector of lower bounds} \item{upper}{vector of upper
#'   bounds} \item{setlower}{logical indicating whether user set lower bounds}
#'   \item{setupper}{logical indicating whether user set upper bounds}
#' @author Jeff Laake
setbounds <- function(lowerbounds,upperbounds,initialvalues,ddfobj){

  # Set values of bounds and check lengths of any user-specified values
  if(!any(is.na(lowerbounds))){
    setlower<-TRUE
  }else{
    setlower<-FALSE
  }
  if(!any(is.na(upperbounds))){
    setupper<-TRUE
  }else{
    setupper<-FALSE
  }

  if(any(is.na(lowerbounds))){
    lowerbounds <- apply(matrix(c(initialvalues - .5*abs(initialvalues),
                                  initialvalues - 1),
                                byrow=FALSE,ncol=2),1,min)
    # fix jll 17-Aug-05 constrain power parameter in hazard rate to be >=1
    if(ddfobj$type=="hr"){
      lowerbounds[1] <- 0
    }else{
      if(ddfobj$type=="gamma") lowerbounds[1] <- -300
    }
  }else{
    if(length(lowerbounds)!=length(initialvalues))
      stop(paste("Error: incorrect number of values for lowerbounds given =",
                 length(lowerbounds),"need ",length(initialvalues),"\n"))
  }
  if(any(is.na(upperbounds))){
    upperbounds <- apply(matrix(c(initialvalues + .5*abs(initialvalues),
                                  initialvalues + 1),
                                  byrow=FALSE,ncol=2),1,max)
  }else{
    if(length(upperbounds)!=length(initialvalues)){
      stop(paste("Error: incorrect number of values for upperbounds given =",
                 length(upperbounds),"need ",length(initialvalues),"\n"))
    }
  }

  return(list(lower = lowerbounds,
              upper = upperbounds,
              setlower=setlower,
              setupper=setupper))
}

