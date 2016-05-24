
##==============================================================================
## Ascendency function
##==============================================================================

Ascendfun <- function(Tij,irow,icol) {

  FlowFrom   <- colSums(Tij)
  FlowTo     <- rowSums(Tij)

    # Sum of all flows, including externals
  Throughput  <- sum(Tij)

  Asc   <- 0  # The ascendency of network
  Overh <- 0  # The overhead   of network
  Cap   <- 0  # The capacity   of network

  for (i in irow)   {
    for (j in icol){
      if (Tij[i,j]>0) {
        Asc   <- Asc   + Tij[i,j]*log2(Tij[i,j]*Throughput/FlowTo[i]/FlowFrom[j])
#       Overh <- Overh - Tij[i,j]*log2(Tij[i,j]*Tij[i,j]  /FlowTo[i]/FlowFrom[j])
        Cap   <- Cap   - Tij[i,j]*log2(Tij[i,j]/Throughput)
       }
     }
  }
  Overh <- Cap-Asc
  return (c(Asc,Overh,Cap))
}


##==============================================================================
##
##   THE ASCENDENCY SUITE
##
##   Calculates Ascendency, Overhead and Capacity
##     for list of rows (irow) and columns (icol) in Tij
##==============================================================================

AscInd <- function (Flow = NULL,               # from-to
                    Tij  = t(Flow),            # to-from
                    Import =NULL,              # flow from external (colNr Tij)
                    Export =NULL,              # flow to external (colNr Tij)
                    Dissipation=NULL)          # external dissipation flow
{

  N          <- InternalNetwork (Tij,Import,Export)

  # number of columns and rows for total flows (Tij)
  ncTij       <- ncol(Tij)
  nrTij       <- nrow(Tij)

  # THROUGHPUT = Sum of all flows, including externals
  Throughput  <- sum(Tij)

  Ascend           <- matrix (nrow=5,ncol=4)
  colnames(Ascend) <- c("Ascendency","Overhead","Capacity","ACratio")
  rownames(Ascend) <- c("Total","Internal","Import","Export","Dissipation")

  # ascendency for the entire network (nrTij: number of rows in Tij)
  Ascend[1,1:3] <- Ascendfun(Tij,1:nrTij,1:ncTij)

  # ascendency for the internal network (iN: indices to internal rows in Tij)
  Ascend[2,1:3] <- Ascendfun(Tij,N$iN,N$jN)

  # ascendency for the import (Import: index to external column in Tij)
  Ascend[3,1:3] <- Ascendfun(Tij,1:nrTij,N$import)

  # ascendency for the export
  Export1       <- setdiff(N$export,Dissipation)       # true export flow
  Ascend[4,1:3] <- Ascendfun(Tij,Export1,1:ncTij)

  # ascendency for the dissipation
  Ascend[5,1:3] <- Ascendfun(Tij,Dissipation,1:ncTij)

  # Ascendency - Capacity Ratio
  Ascend[,4]    <-  Ascend[,1]/  Ascend[,3]

  return(Ascend)
}
