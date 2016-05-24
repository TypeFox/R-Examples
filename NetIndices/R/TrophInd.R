
##==============================================================================
##
## Trophic analysis function
##
##==============================================================================

TrophInd <- function (Flow = NULL,             # from-to
                      Tij  = t(Flow),          # to-from
                      Import =NULL,            # flow from external (colNr Tij)
                      Export =NULL,            # flow to external (colNr Tij)
                      Dead=NULL)               # flow to dead matter


{
  if (is.character(Dead))
     dead <- which(rownames(Tij)%in%Dead) else
     dead <- Dead
  if (length(dead) != length(Dead)) stop("Dead not recognized")

## Check input and calculate internal network
  N          <- InternalNetwork (Tij,Import,Export)
  p          <- Diet(N$Tint,dead,N$iN)
  ncomp      <- ncol(N$Tint)    # the number of compartments (without externals)

##
## Trophic level TL  TL(i) = 1 + sumj(pij*TL(j))
##                   TL(i)-sumj(pij*TL(j))= 1
##
  A       <- -p
#  diag(A) <- diag(A)+ 1  # was this
# Karline: changed Januari 2012: organisms that "feed on themselves" are problematic!
  diag(A) <- 1
  B       <- rep(1,ncomp)
  TL      <- ginv(A) %*% B

##
## Omnivory index: variance of trophic levels of preys
##

  OI <- vector (length = ncomp)
  for (i in 1:ncomp)
    OI[i] <-sum((TL-(TL[i]-1))^2*p[i,])

  return(data.frame(TL, OI,row.names=rownames(N$Tint)))
}               ## END Trophic

