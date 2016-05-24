
##==============================================================================
##
## FINN's SUITE   - Pathway analysis
##
##==============================================================================

PathInd <- function (Flow = NULL,            # from-to
                     Tij  = t(Flow),         # to-from
                     Import =NULL,           # flow from external (colNr Tij)
                     Export =NULL)           # flow to external (colNr Tij)
{
  N          <- InternalNetwork (Tij,Import,Export)

  # Rate of change of each compartment
  RateComp   <- N$FlowToC-N$FlowFromC

  ncTij       <- ncol(Tij)
  nrTij       <- nrow(Tij)
  ncomp       <- ncol(N$Tint)
  compNames   <- rownames(N$Tint)

  ExportSum   <- sum(N$FlowTo[N$export])
  ImportSum   <- sum(N$FlowFrom[N$import])

## THROUGHFLOW - Rates of change taken into consideration
## Throughflow based on rows and columns has to be the same

  Throughflow <- sum(N$Tint) + ImportSum - sum(RateComp[RateComp<0])
  Throughput  <- sum(Tij)

##
## THE PATHLENGTH
##

  Pathlength   <- Throughflow/ (ExportSum + sum(RateComp[RateComp>0]))

##
## TRANSITIVE CLOSURE MATRIX
##

  CompThroughflow <- pmax(N$FlowFromC,N$FlowToC)  # total flow, internal compartments
  Qij <- matrix(nrow=ncomp,ncol=ncomp,0)

  for (i in 1:ncomp)
    Qij[i,]  <- N$Tint[i,]/ CompThroughflow[i]

  diagnl    <- diag(nrow=ncomp,1)   # unity matrix
  IQ        <- diagnl-Qij
  M         <- ginv(IQ)           # The generalised inverse

  # Cycled throughflow
  diaM      <- diag(M)
  TSTC      <- sum((1-1/diaM)*N$FlowFromC)

  # Noncycled throughflow
  TSTS      <- Throughflow - TSTC

  # Finn's cycling index
  FCI       <- TSTC/Throughflow

  # Finn's cycling index revisited
  FCIb      <- TSTC/Throughput

  Finn<-list(TSTC=TSTC,TSTS=TSTS,FCI=FCI,FCIb=FCIb, APL=Pathlength)
  return(Finn)
}
