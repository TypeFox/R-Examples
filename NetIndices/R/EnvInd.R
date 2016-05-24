
##==============================================================================
##
##   THE ENVIRONS SUITE
##
##==============================================================================

EnvInd <- function (Flow = NULL,             # from-to
                        Tij  = t(Flow),      # to-from
                        Import =NULL,        # flow from external (colNr Tij)
                        Export =NULL,        # flow to external (colNr Tij)
                        full=FALSE)          # if True, also returns matrices
{

  N          <- InternalNetwork (Tij,Import,Export)

  ncTij       <- ncol(Tij)
  nrTij       <- nrow(Tij)
  ncomp       <- ncol(N$Tint)
  compNames   <- rownames(N$Tint)
  ImportSum   <- sum(N$FlowFrom[N$import])
  RateComp    <- N$FlowToC-N$FlowFromC
##
## Network Aggradation - ratio of internal disorder to generated disorder
##

  Throughflow  <- sum(N$Tint) + ImportSum - sum(RateComp[RateComp<0])
  TotalOutflow <- sum(Tij[N$export,N$jN])
  Naggr        <- Throughflow/TotalOutflow

##
## Homogenization Index - measure of evenness of flows
##
  # G1: transitive closure matrix
  G   <- matrix(nrow=ncomp,ncol=ncomp,0)

  for (j in 1:ncomp)
    G[,j] <- N$Tint[,j]/N$FlowFromC[j]

  # N1: Integral nondimensional matrix
  N1 <- NULL
  if (any(is.nan(G)))
    stop("Some of the flows are NANs")
  I   <- diag(nrow=ncomp,1)
  IG  <- I-G
  N1  <- ginv(IG)

  # coefficient of variation of G and N1
  # note: sd of a matrix takes standard deviation per column -
  # overrule by as.vector

  MG  <- mean(G)
  MN  <- mean(N1)
  CVG <- sd(as.vector(G))/MG
  CVN <- sd(as.vector(N1))/MN

  HP  <- CVG/CVN  # Homogenization index
##
## Dominance of Indirect effects index- ratio i/d
##

  D1  <- N$Tint-t(N$Tint)
  ID  <- (sum(N1-I-G))/sum(G)    #delta(ij)=1 diagonal elements!

##
## Utility Matrix (UP)
##

  D     <- D1/N$FlowToC
  U     <- ginv(I-D)
  # synergism index - was wrong in Latham 2006;
  # here is the correct(?) value:
  Gamma <- N$FlowToC*U
  bc    <- sum(Gamma[Gamma>0])/abs(sum(Gamma[Gamma<0]))


  if (full) {
    rownames(U) <- rownames(N1) <-rownames(G) <-compNames
    colnames(U) <- colnames(N1) <-colnames(G)<- compNames

    Res <- list(NAG=Naggr,HP=HP,BC=bc,ID=ID,CVN=CVN,CVG=CVG, U=U,N1=N1,G=G)

  } else
    Res<-list(NAG=Naggr,HP=HP,BC=bc,ID=ID,MN=MN,MG=MG,CVN=CVN,CVG=CVG)

  return(Res)
}
