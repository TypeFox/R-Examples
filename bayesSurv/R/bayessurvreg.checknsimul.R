###########################################################
#### AUTHOR:     Arnost Komarek                        ####
####             (06/02/2007)                          ####
####                                                   ####
#### FILE:       bayessurvreg.checknsimul.R            ####
####                                                   ####
#### FUNCTIONS:  bayessurvreg.checknsimul              ####
###########################################################

bayessurvreg.checknsimul <- function(nsimul)
{
  if(length(nsimul) == 0) innsimul <- "arnost"
  else                    innsimul <- names(nsimul)

  tmp <- match("niter", innsimul, nomatch=NA)
  if(is.na(tmp)) stop("nsimul$niter must be given")
  if (nsimul$niter <= 0) stop("nsimul$niter must be positive")
  
  tmp <- match("nburn", innsimul, nomatch=NA)
  if(is.na(tmp)) stop("nsimul$nburn must be given")
  if (nsimul$nburn < 0) stop("nsimul$nburn must be non-negative")
  if (nsimul$nburn > nsimul$niter) stop("nsimul$nburn must not be higher than nsimul$niter")

  tmp <- match("nwrite", innsimul, nomatch=NA)
  if(is.na(tmp)) nsimul$nwrite <- nsimul$niter
  
  tmp <- match("nthin", innsimul, nomatch=NA)
  if(is.na(tmp)) nsimul$nthin <- 1

  ## Only needed by bayessurvreg1:
  tmp <- match("nnoadapt", innsimul, nomatch=NA)
  if(is.na(tmp)) nsimul$nnoadapt <- 0
  
  return(nsimul)  
}
