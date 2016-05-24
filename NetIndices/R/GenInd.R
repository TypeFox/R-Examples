
##==============================================================================
## GENERAL SUITE - COMPARTMENTS, TOTAL SYSTEM THROUGHPUT (T..),
##           TOTAL SYSTEM THROUGHFLOW (TST), LINK DENSITY (LD),
##           NUMBER OF LINKS (L), AV. COMPARTMENT THROUGHFLOW (TSTbar),
##           CONNECTANCE (C), AVERAGE LINK WEIGHT (Tijbar),
##           COMPARTMENTALIZATION (Cbar)
##==============================================================================

GenInd <- function (Flow = NULL,            # from-to
                    Tij  = t(Flow),         # to-from
                    Import =NULL,           # flow from external (colNr Tij)
                    Export =NULL,           # flow to external (colNr Tij)
                    tol=0)                  # flow<=tol is assumed absent
{                      

##------------------------------------------------------------------------
## Flow is a matrix with Flow[i,j] the flow from i (row) to j (column)
## Tij[i,j] contains flow from j to i
## note: component position in rows and columns must be the same - not checked
##------------------------------------------------------------------------

  N          <- InternalNetwork (Tij,Import,Export)

  # Rate of change of each compartment
  RateComp   <- N$FlowToC-N$FlowFromC

  # number of columns and rows for total flows (Tij)
  ncTij       <- ncol(Tij)
  nrTij       <- nrow(Tij)
  ncomp       <- ncol(N$Tint)   # the number of compartments (without externals)
  compNames   <- rownames(N$Tint)

 
##
## NUMBER OF TOTAL AND INTERNAL LINKS AND LINK DENSITY
##
  intlinks    <- length(which(N$Tint>tol))
  links       <- length(which(Tij >tol))
  LD          <- links/ncomp

##
## THROUGHPUT AND THROUHFLOW
##
  ExportSum   <- sum(N$FlowTo[N$export])
  ImportSum   <- sum(N$FlowFrom[N$import])

  # THROUGHFLOW - Rates of change taken into consideration
  # Throughflow based on rows and columns has to be the same

  Throughflow <- sum(N$Tint) + ImportSum - sum(RateComp[RateComp<0])
  # or:  Throughflow <- sum(N$Tint) + ExportSum + sum(RateComp[RateComp>0])

  # THROUGHPUT = Sum of all flows, including externals

  Throughput  <- sum(Tij)

##
## AVERAGE COMPARTMENT THROUGHFLOW (TSTbar)
##

  Avthrflow <- Throughflow/ncomp

##
## CONNECTANCE
##
    
  Connectance <- intlinks/ncomp/(ncomp-1)

##
## AVERAGE LINK WEIGHT
##

  Avlinkweight <- Throughput/links

##
## COMPARTMENTALIZATION
##
  linkmat             <- N$Tint
  linkmat[linkmat>0]  <- 1   # 1 if there is a link, 0 elsewhere

## The number of components with which both i and j interact divided by
## the number of species by which either i or j interact
  Cij <- matrix(nrow=ncomp,ncol=ncomp,0)
  for (i in 1:ncomp) {
    int_i <- union(which(linkmat[i,]>0),which(linkmat[,i]>0))
    for (j in 1:ncomp) {
      int_j <- union(which(linkmat[j,]>0),which(linkmat[,j]>0))
      sect <- intersect(int_i,int_j)
      uni  <- union    (int_i,int_j)
      Cij[i,j]  <- length(sect)/length(uni)
    }
  }
  Compart <- (sum(Cij)-ncomp) / ncomp /(ncomp-1)


##    "No. of Compartments (n)","Total System Throughput (T..)",
##    "Total System Throughflow (TST)","Link Density (LD)",
##    "No. of Internal links, (Lint)","Total No. of links (Ltot)",
##    "Av. Compartment Throughflow (TSTbar)","Connectance (internal)(C)",
##    "Av. Link Weight (Tijbar)","Compartmentalization (Cbar)"
        
  list(N=ncomp,T..=Throughput,TST=Throughflow,
      Lint=intlinks,Ltot=links,LD=LD,C=Connectance,
      Tijbar=Avlinkweight,TSTbar=Avthrflow,Cbar=Compart)

}              #end generalIndices
