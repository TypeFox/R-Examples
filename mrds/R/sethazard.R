sethazard <- function(ddfobj,dmat,width){
  # temporary function to evaluate the likelihood for hazard rate
  evallike<-function(fpar){
    ddfobj$intercept.only <- TRUE
    ddfobj$shape$parameters[1] <- fpar[1]
    ddfobj$scale$parameters[1] <- fpar[2]
    # set adjustment parameters to zero for now...
    if(!is.null(ddfobj$adjustment))
      ddfobj$adjustment$parameters <- rep(0,length(ddfobj$adjustment$order))

    fpar <- getpar(ddfobj)
    ddfobj$shape$parameters[1] <- NA
    ddfobj$scale$parameters[1] <- NA
    if(!is.null(ddfobj$adjustment)){
      ddfobj$adjustment$parameters <- rep(NA,length(ddfobj$adjustment$order))
    }

    flnl(fpar, ddfobj, FALSE, misc.options=list(width=width,
                                                int.range=c(0,width),
                                                showit=0,
         integral.numeric=FALSE,standardize=FALSE,fitting="none",point=FALSE))
  }

  # Using code from CDS in Distance.
  #  First find the scale parameter as you (CDS) would for half-normal,
  #  by finding:
  ss <- sqrt(sum(dmat$distance^2)/length(dmat$distance))

  # Then find the shape parameter by using the following algorithm:
  for(j in 1:10){
    scale <- log(ss/2+j*0.2*ss)
    for(i in seq(1,30,2)){
      shape <- log(i)
      # evaluate likelihood at this value
      newlnl <- evallike(c(shape,scale))
      if(is.nan(newlnl)){
        break
      }
      if(i==1 & j==1){
        lnl <- newlnl
        ival <- 2
        jval <- 1
      }else{
        if(newlnl<lnl){
          lnl <- newlnl
          ival <- i
          jval <- j
        }
      }
    }
  }
  if(is.nan(newlnl)){
    shape <- log(2)
  }else{
    scale <- log(ss/2+jval*0.2*ss)
    shape <- log(ival)
  }
  return(list(shape=shape,
              scale=scale))
}
