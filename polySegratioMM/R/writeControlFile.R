`writeControlFile` <-
function(jags.control,
                             file=paste(jags.control$stem,".cmd",sep=""))
{

  if (class(jags.control) != "jagsControl")
    stop("'jags.control' must be of class 'jagsControl'")

##  if (class(priors) != "priorsSegratioMM")
##    stop("'priors' must be of class 'priorsSegratioMM'")

  cat(file=file, jags.control$jags.code, sep="\n")

}

