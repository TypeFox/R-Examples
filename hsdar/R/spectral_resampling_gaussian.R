
get.gaussian.response <- function(fwhm)
{
  if (is.null(attr(fwhm, "fwhm")))
  {
    if (any(toupper(names(fwhm))=="FWHM"))
    {
      fwhm_vec <- fwhm[, which(toupper(names(fwhm))=="FWHM")]
    } else {
      fwhm_vec <- if (pmatch("FWHM", toupper(names(fwhm)))==0) NULL else fwhm[, pmatch("FWHM", toupper(names(fwhm)))]
    }
    if (any(toupper(names(fwhm))=="CENTER"))
    {
      centerwl <- fwhm[, which(toupper(names(fwhm))=="CENTER")]
    } else {
      centerwl <- if (pmatch("CENTER", toupper(names(fwhm)))==0) NULL else fwhm[, pmatch("CENTER", toupper(names(fwhm)))]
    }
    if (any(c(is.null(fwhm_vec), is.null(centerwl))))
    {
      lb<-fwhm[,2]
      ub<-fwhm[,3]
      centerwl <- lb + (ub - lb)/2
      fwhm_vec <- (centerwl - lb) * 2
    }
    fwhm <- data.frame(channel=c(1:length(centerwl)),center=centerwl,fwhm=fwhm_vec)
  } else {
    if (!attr(fwhm, "fwhm"))
    {
      lb <- fwhm[,attr(fwhm, "50pass")[1]]
      ub <- fwhm[,attr(fwhm, "50pass")[2]]
      centerwl <- lb + (ub - lb)/2
      fwhm_vec <- (centerwl - lb) * 2
      fwhm <- data.frame(No=c(1:length(centerwl)), center=centerwl, fwhm=fwhm_vec)
    }
  }
  
  lb <- fwhm[,2]-fwhm[,3]/2
  ub <- fwhm[,2]+fwhm[,3]/2
  nch <- nrow(fwhm)
  responsedim <- c(min(lb)-(ub[which(lb==min(lb))]-lb[which(lb==min(lb))]),
                   max(ub)+(ub[which(ub==max(ub))]-lb[which(ub==max(ub))]),
                   1)
  response <- matrix(data = 0, ncol = nch, nrow =responsedim[2]-responsedim[1]+1)
  range_wl <- seq.int(responsedim[1],responsedim[2],responsedim[3])
  for (i in 1:ncol(response))
  {
    gauss <- dnorm(range_wl, mean = mean(c(lb[i],ub[i])), sd = (ub[i]-lb[i])/2)
    gauss <- (gauss-min(gauss))/(max(gauss)-min(gauss))
    response[,i] <- gauss
  }
  response <- as.data.frame(response)
  names(response) <- paste("Band",c(1:nch),sep="_")
  attr(response,"minwl") <- responsedim[1]
  attr(response,"maxwl") <- responsedim[2]
  attr(response,"stepsize") <- responsedim[3]
  return(response)
}