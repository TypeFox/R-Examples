

###  $Id: funcGenerators.R,v 1.3 2008/02/05 20:21:23 goswami Exp $
###  
###  File:    funcGenerators.R
###  Package: EMC
###  
###  Copyright (C) 2006-present Gopi Goswami
###
###  This program is free software; you can redistribute it and/or modify
###  it under the terms of the GNU General Public License as published by
###  the Free Software Foundation; either version 2 of the License, or
###  (at your option) any later version.
###
###  This program is distributed in the hope that it will be useful,
###  but WITHOUT ANY WARRANTY; without even the implied warranty of
###  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###  GNU General Public License for more details.
###
###  For a copy of the GNU General Public License please write to the
###  Free Software Foundation, Inc.
###  59 Temple Place, Suite 330.
###  Boston, MA  02111-1307 USA.
###
###  For bugs in the code please contact:
###  <goswami@stat.harvard.edu>
###
###  SYNOPSIS
###
###
###
###  DESCRIPTION
###
###
###


### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## The Cigar-shaped distribution with (symmetric) normal-proposal
CigarShapedFuncGenerator1 <-
    function (seed)
{
    set.seed(seed)    
    dd     <- 2
    ARDisp <-
        function (rho)
        {
            tmp <- rep(1, dd)
            diag((1 - rho) * tmp) + rho * tmp %*% t(tmp)
        }

    means          <- c(1, 1)
    disp           <- ARDisp(-0.95)
    logTarDensFunc <-
        function (draw, ...)
            dmvnorm(draw, means, disp, log = TRUE)

    proposalSD  <- c(1, 2)
    propNewFunc <-
        function (block, currentDraw, ...)
        {
            proposalDraw        <- currentDraw
            proposalDraw[block] <- rnorm(1, currentDraw[block], proposalSD[block])
            proposalDraw
        }            
    
    list(logTarDensFunc  = logTarDensFunc,
         propNewFunc     = propNewFunc)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## The Cigar-shaped distribution with independent t-proposal
CigarShapedFuncGenerator2 <-
    function (seed)
{
    set.seed(seed)    
    dd     <- 2
    ARDisp <-
        function (rho)
        {
            tmp <- rep(1, dd)
            diag((1 - rho) * tmp) + rho * tmp %*% t(tmp)
        }

    means          <- c(1, 1)
    disp           <- ARDisp(-0.95)
    logTarDensFunc <-
        function (draw, ...)
            dmvnorm(draw, means, disp, log = TRUE)

    tDF         <- 3
    propNewFunc <-
        function (block, currentDraw, ...)
        {
            proposalDraw        <- currentDraw
            proposalDraw[block] <- rt(1, tDF)
            proposalDraw
        }            
    
    logPropDensFunc <-
        function (block, currentDraw, proposalDraw, ...)
            dt(proposalDraw[block], tDF, log = TRUE)

    list(logTarDensFunc  = logTarDensFunc,
         propNewFunc     = propNewFunc,
         logPropDensFunc = logPropDensFunc)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
## The V-shaped distribution
VShapedFuncGenerator <-
    function (seed)
{
    set.seed(seed)    
    dd     <- 2
    ARDisp <-
        function (rho)
        {
            tmp <- rep(1, dd)
            diag((1 - rho) * tmp) + rho * tmp %*% t(tmp)
        }

    nMixComps  <- 2
    logWeights <- log(rep(1 / nMixComps, nMixComps))
    meanMat    <- matrix(c(1, 1, 15, 1), nMixComps, byrow = TRUE)
    dispParam  <- c(-0.95, 0.95)
    dispArr    <- array(dim = c(2, 2, nMixComps))
    for (ii in seq_len(nMixComps)) {
        dispArr[ , , ii] <- ARDisp(dispParam[ii])
    }
    logTarDensFunc <-
        function (draw, ...)
        {
            ld <- sapply(seq_len(nMixComps), FUN =
                         function (ii)
                     {
                         dmvnorm(draw, meanMat[ii, ], dispArr[ , , ii], log = TRUE)
                     })
            ww <- logWeights + ld
            mm <- max(ww)
            mm + log(sum(exp(ww - mm)))
        }
    
    MHProposalSD  <- c(1.0, 1.0)
    MHPropNewFunc <-
        function (temperature, block, currentDraw, ...)
        {
            proposalDraw        <- currentDraw
            proposalDraw[block] <- rnorm(1, currentDraw[block],
                                         sqrt(temperature) * MHProposalSD[block])
            proposalDraw
        }

    list(logTarDensFunc = logTarDensFunc,
         MHPropNewFunc  = MHPropNewFunc)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## The W-shaped distribution
WShapedFuncGenerator <-
    function (seed)
{
    set.seed(seed)    
    dd     <- 2
    ARDisp <-
        function (rho)
        {
            tmp <- rep(1, dd)
            diag((1 - rho) * tmp) + rho * tmp %*% t(tmp)
        }

    nMixComps  <- 4
    logWeights <- log(rep(1 / nMixComps, nMixComps))
    meanMat    <- matrix(c(-12, 1, -4, 1, 4, 1, 12, 1), nMixComps, byrow = TRUE)
    dispParam  <- c(-0.95, 0.95, -0.95, 0.95)
    dispArr    <- array(dim = c(2, 2, nMixComps))
    for (ii in seq_len(nMixComps)) {
        dispArr[ , , ii] <- ARDisp(dispParam[ii])
    }
    logTarDensFunc <-
        function (draw, ...)
        {
            ld <- sapply(seq_len(nMixComps), FUN =
                         function (ii)
                     {
                         dmvnorm(draw, meanMat[ii, ], dispArr[ , , ii], log = TRUE)
                     })
            ww <- logWeights + ld
            mm <- max(ww)
            mm + log(sum(exp(ww - mm)))
        }
    
    MHProposalSD  <- c(0.8, 0.8)
    MHPropNewFunc <-
        function (temperature, block, currentDraw, ...)
        {
            proposalDraw        <- currentDraw
            proposalDraw[block] <- rnorm(1, currentDraw[block],
                                         sqrt(temperature) * MHProposalSD[block])
            proposalDraw
        }

    list(logTarDensFunc = logTarDensFunc,
         MHPropNewFunc  = MHPropNewFunc)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### The unimodal normal distribution
uniModeFuncGenerator <-
    function (seed)
{
    set.seed(seed)
    dd             <- 2
    mean           <- rep(0, dd)
    disp           <- diag(dd)
    logTarDensFunc <-
        function (draw, ...)
            dmvnorm(draw, mean, disp, log = TRUE)
    
    MHProposalSD  <- c(2.5, 2.5)
    MHPropNewFunc <-
        function (temperature, block, currentDraw, ...)
        {
            proposalDraw        <- currentDraw
            proposalDraw[block] <- rnorm(1, currentDraw[block],
                                         sqrt(temperature) * MHProposalSD[block])
            proposalDraw
        }
    
    list(logTarDensFunc = logTarDensFunc,
         MHPropNewFunc  = MHPropNewFunc)    
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### The twenty-mode problem
twentyModeFuncGenerator <-
    function (seed)
{
    set.seed(seed)
    params     <- read.table(file = 'twentyModeParams.txt', header = F)
    nMixComps  <- 20
    logWeights <- log(params[ , 1])
    meanMat    <- as.matrix(params[ , c(2, 3)])    
    dispArr    <- array(dim = c(2, 2, nMixComps))
    for (ii in seq_len(nMixComps)) {
        dispArr[ , , ii] <- params[ii, 4]^2 * diag(2)
    }    
    logTarDensFunc <-
        function (draw, ...)
        {
            ld <- sapply(seq_len(nMixComps), FUN =
                         function (ii)
                     {
                         dmvnorm(draw, meanMat[ii, ], dispArr[ , , ii], log = TRUE)
                     })
            ww <- logWeights + ld
            mm <- max(ww)
            mm + log(sum(exp(ww - mm)))
        }
    
    MHProposalSD  <- c(0.25, 0.25)
    MHPropNewFunc <-
        function (temperature, block, currentDraw, ...)
        {
            proposalDraw        <- currentDraw
            proposalDraw[block] <- rnorm(1, currentDraw[block],
                                         sqrt(temperature) * MHProposalSD[block])
            proposalDraw
        }
    
    list(logTarDensFunc = logTarDensFunc,
         MHPropNewFunc  = MHPropNewFunc)    
}


### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### The three dimensional Normal example

threeDimNormalFuncGenerator <-
    function (seed)
{
    set.seed(seed)
    dd <- 3
    mu <- rep(0, dd)
    Sigma <- matrix(c(1,   0.9,   0,
                      0.9,   1,   0,
                      0,     0,   1), nrow = dd, ncol = dd)

    logTarDensFunc <-
        function (draw, ...)
        {
            dmvnorm(draw, mu, Sigma, log = TRUE)
        }

    proposalInfo  <-
        list(list(pos       =
                  c(1, 2),
                  propSigma =
                  matrix(c(4, 2,
                           2, 4), nrow = 2, ncol = 2)),
             
             list(pos    = 3,
                  propSD = 2))
    propNewFunc <-
        function (block, currentDraw, ...)
        {
            proposalDraw <- currentDraw
            info         <- proposalInfo[[block]]
            pos          <- info$pos
            if (block == 1) {
                proposalDraw[pos] <- mvrnorm(1, mu = currentDraw[pos],
                                             Sigma = info$propSigma)
            } else if (block == 2) {
                proposalDraw[pos] <- rnorm(1, currentDraw[pos], info$propSD)
            }
            return(proposalDraw)
        }
    
    list(logTarDensFunc = logTarDensFunc,
         propNewFunc    = propNewFunc)
}
