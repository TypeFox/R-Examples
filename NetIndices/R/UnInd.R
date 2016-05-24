
##==============================================================================
##
## THE AMI SUITE - AMI, STATISTICAL UNCERTAINTY (HR), CONDITIONAL UNCERTAINTY (DR)
## REALIZED UNCERTAINTY (AMI/HR)
##
##==============================================================================
UncInd <- function (Flow = NULL,                      # from-to
                    Tij  = t(Flow),                   # to-from
                    Import =NULL,                     # flow from external (colNr Tij)
                    Export =NULL)                     # flow to external (colNr Tij)
{

  N          <- InternalNetwork (Tij,Import,Export)

  # number of columns and rows for total flows (Tij)
  ncTij       <- ncol(Tij)
  nrTij       <- nrow(Tij)
  ncomp       <- ncol(N$Tint)          # the number of compartments (without externals)
  compNames   <- rownames(N$Tint)      # names of internal compartment

   # Sum of all flows, including externals
  Throughput  <- sum(Tij)

##
## AVERAGE MUTUAL INFORMATION
##

  # ascendency for the entire network (nrTij: number of rows in Tij)
  Ascend <- Ascendfun(Tij,1:nrow(Tij),1:ncol(Tij))

  AMI        <- Ascend[[1]] / Throughput

##
## STATISTICAL UNCERTAINTY
##

  Q          <- N$FlowFrom/Throughput
  Q          <- Q[Q>0]
  HR         <- -sum(Q*log2(Q))

##
## CONDITIONAL UNCERTAINTY INDEX
##

  DR         <- HR - AMI

##
## REALIZED UNCERTAINTY INDEX
##

  RU         <- AMI/HR

##
## MAXIMUM UNCERTAINTY (Hmax)
##

  Hmax <- ncomp*log2(nrTij)

##
## CONSTRAINT INFORMATION (Hc)
##

  blsum <- 0
  for (i in 1:nrTij) {     # all the rows (including externals)
    for (j in N$jN) {    # only internal columns
      if (Tij[i,j]>0) blsum <- blsum +
                      (Tij[i,j]/N$FlowFrom[j])*log2(Tij[i,j]/N$FlowFrom[j])
    }
  }

  Hc <- Hmax + blsum
  names(Hc)<-NULL
##
## CONSTRAINT EFFICIENCY
##

  CE <- Hc/Hmax
##
## NETWORK EFFICIENCY
##

  Hsys <- Hmax - Hc

  list(AMI=AMI, HR=HR, DR=DR, RU=RU,Hmax=Hmax,Hc=Hc,
    Hsys=Hsys,CE=CE)

}
