
##==============================================================================
##
## THE EFFECTIVE MEASURES SUITE
##
##==============================================================================

EffInd<- function (Flow = NULL,               # from-to
                   Tij  = t(Flow),            # to-from
                   Import =NULL,              # flow from external (colNr Tij)
                   Export =NULL)              # flow to external (colNr Tij)
{

  N          <- InternalNetwork (Tij,Import,Export)

  # Sum of all flows, including externals

  Throughput  <- sum(Tij)
  ncomp       <- ncol(N$Tint)          # the number of compartments (without externals)

  CZ <-1    # EFFECTIVE CONNECTANCE
  FZ <-1    # EFFECTIVE FLOWS
  NZ <-1    # EFFECTIVE NODES
  RZ <-1    # EFFECTIVE ROLES

    for (i in 1:ncomp) {
      if (N$FlowToC[i]>0) {
         for (j in 1:ncomp) {
           if (N$FlowFromC[j]>0)
            NZ <-NZ*(Throughput^2        /
               N$FlowFromC[j]/N$FlowToC[i])^( 0.5*N$Tint[i,j]/Throughput)
         } # end j
      } # end if flowtoC
    } # end i
    for (i in 1:ncomp) {
      for (j in 1:ncomp) {
        if (N$Tint[i,j]>0) {
          CZ <-CZ*(N$Tint[i,j]^2         /N$FlowFromC[j]/
              N$FlowToC[i])^(-0.5*N$Tint[i,j]/Throughput)
          RZ <-RZ*(N$Tint[i,j]*Throughput/N$FlowFromC[j]/
              N$FlowToC[i])^(     N$Tint[i,j]/Throughput)
          FZ <-FZ*(N$Tint[i,j]           /Throughput)^
              (-N$Tint[i,j]/Throughput)
        }
      }
    }
  names(CZ) <- names(RZ) <- names(NZ) <- names(FZ) <- NULL
  list(CZ=CZ,FZ=FZ,NZ=NZ,RZ=RZ)
}
