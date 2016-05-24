
#########################################################################
####
#### R-script used to generate the monte carlo sample
#### as used in the article:
####
#### Kones, J., Soetaert, K., van Oevelen D. and J.O. Owino, 2009.
#### Are network indices robust indicators of food web functioning?
#### a Monte Carlo approach
#### Ecological Modelling. In press.
####
#### Script written by Karline Soetaert
####
#### Note:
#### In the paper, the monte carlo procedure generates 3000 samples;
#### As this takes a while, here it is done based on 100 samples !
####
#########################################################################

# packages to be loaded...
require(LIM)            # contains the food web specifications
require(NetIndices)     # contains the network indices

# number of samples
nsamp <- 3000

# smaller number is faster...
# toggle off (by inserting "#" in front) to have it based on 3000 samples
nsamp <- 100

#==============================================================================
# part 1   - create networks: parsimoniuos solution and random sample
#==============================================================================

sample.lim<-function( lim, ... )  xsample(E=lim$A,F=lim$B,G=lim$G,H=lim$H,...)

# Parsimonious solutions
RigaAutumn <- Ldei.lim(LIMRigaAutumn)    # Gulf of riga, autumn
RigaSummer <- Ldei.lim(LIMRigaSummer)    # Gulf of riga, summer
RigaSpring <- Lsei.lim(LIMRigaSpring,parsimonious=TRUE)    # Gulf of riga, spring
Takapoto   <- Ldei.lim(LIMTakapoto )     # Takapoto atoll

# randomly sampled solutions of the four food webs - jump lenghts ~ convergence
RigaAutumn.X <- Xsample(LIMRigaAutumn,iter=nsamp,jmp=20)
RigaSpring.X <- Xsample(LIMRigaSpring,iter=nsamp,jmp=20)
RigaSummer.X <- Xsample(LIMRigaSummer,iter=nsamp,jmp=20)
Takapoto.X   <- Xsample(LIMTakapoto  ,iter=nsamp,jmp=500)

# Check and plot the output
panel.hist <- function(x, ...)
      {
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(usr[1:2], 0, 1.5) )
          h <- hist(x, plot = FALSE)
          breaks <- h$breaks; nB <- length(breaks)
          y <- h$counts; y <- y/max(y)
          rect(breaks[-nB], 0, breaks[-1], y, col="grey", ...)
      }

panel.trace <- function(x, ...)
      {
          usr <- par("usr"); on.exit(par(usr))
          nn  <- min(length(x),1000)
          par(usr = c(0, nn, usr[1:2]) )
          lines(x[1:nn],type="l")
      }
# Note: these figures are used to check convergence.
# Note: in the paper it is based on 3000 random solutions, which of course
# is more robust.

pairs(RigaAutumn.X, pch='.',diag.panel=panel.hist,gap=0,
      xaxt="n",yaxt="n",upper.panel=NULL,cex=2)
pairs(RigaAutumn.X, pch='.',diag.panel=panel.trace,gap=0,
      xaxt="n",yaxt="n",upper.panel=NULL,cex=2)
pairs(RigaSpring.X, pch='.',diag.panel=panel.trace,gap=0,
      xaxt="n",yaxt="n",upper.panel=NULL,cex=2)
pairs(RigaSummer.X, pch='.',diag.panel=panel.trace,gap=0,
      xaxt="n",yaxt="n",upper.panel=NULL,cex=2)
pairs(Takapoto.X, pch='.',diag.panel=panel.trace,gap=0,
      xaxt="n",yaxt="n",upper.panel=NULL,cex=2)

#==============================================================================
# part 2   - estimate indices
#==============================================================================

Indices <- function(LIM.web,    # web specifications
                    Pars.web,   # parsimonious solution
                    Sample.web, # random samples
                    ext,        # external compartments (import and export)
                    dead,       # dead compartments (TL=1)
                    specs)      # species for which trophic analysis is wanted
{
  # binary flow matrix
  flowmat    <- LIM.web$Flowmatrix
  flowmatrix <- flowmat
  ii         <- which(flowmat > 0, arr.ind = TRUE)

  # Add parsimonious solution to random sample
  X          <- rbind(t(Pars.web$X),Sample.web)

  Indices    <- NULL  # will contain all network indices
  indices    <- NULL  # a subset

  for (i in 1:nrow(X))
   {
    # generate required flowmatrix for this solution
     flowmatrix[ii] <- X[i,flowmat[ii]]

    # calculate all network indices
    UU<-UncInd    (flowmatrix,Import=ext,Export=ext)
    EF<-EffInd    (flowmatrix,Import=ext,Export=ext)
    AA<-AscInd    (flowmatrix,Import=ext,Export=ext)
    EE<-EnvInd    (flowmatrix,Import=ext,Export=ext)
    GG<-GenInd    (flowmatrix,Import=ext,Export=ext)
    PP<-PathInd   (flowmatrix,Import=ext,Export=ext)
    TT<-TrophInd  (flowmatrix,Import=ext,Export=ext,Dead=dead)

    # select the indices that will be investigated
    Ind <- unlist(c(UU,AA[1,],EE,GG,PP,TT,EF))

    ind <- c(Ind["TST"],Ind["Cbar"],Ind["Ascendency"],Ind["Overhead"],
             Ind["Capacity"],Ind["ACratio"],Ind["AMI"],Ind["HR"],
             Ind["Hsys"],Ind["CZ"],Ind["FZ"],Ind["NZ"],Ind["RZ"],
             Ind["FCI"],Ind["CVN"],Ind["CVG"],Ind["HP"])
    ind <- c(ind,TL=TT$TL[specs])
    ind <- c(ind,OI=TT$OI[specs])

 # and add them to the results matrices
    Indices <- rbind(Indices,Ind)
    indices <- rbind(indices,ind)
  }
  list(Full=Indices,   # all calculated indices
       Sub=indices)    # selected set
}

# Here they are estimated and saved to file; as this takes quite some time,
# it is done only once...
# riga:  P1, P2, Bact, Nano, Zoo, Det, DOC

IndRigaAutumn  <- Indices(LIMRigaAutumn,RigaAutumn,RigaAutumn.X,
               ext=c("CO2","Sedimentation"),dead=c(6,7),specs=c(4,5))

IndRigaSpring  <- Indices(LIMRigaSpring,RigaSpring,RigaSpring.X,
               ext=c("CO2","Sedimentation"),dead=c(6,7),specs=c(4,5))

IndRigaSummer  <- Indices(LIMRigaSummer,RigaSummer,RigaSummer.X,
               ext=c("CO2","Sedimentation"),dead=c(6,7),specs=c(4,5))

# Phytoplankton, Bacteria, Protozoa, Microzooplankton,
# Mesozooplankton, Detritus, DOC
IndTakapoto    <- Indices(LIMTakapoto,Takapoto,Takapoto.X,
               ext=c("CO2","Sedimentation","Grazing"),dead=c(6,7),specs=c(3,4,5))

#save(file="Indices",IndRigaAutumn,IndRigaSpring,IndRigaSummer,IndTakapoto,
#     RigaAutumn,RigaAutumn.X,RigaSpring,RigaSpring.X,RigaSummer,RigaSummer.X,
#     Takapoto,Takapoto.X)
