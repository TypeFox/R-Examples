### - nlogL method will return the negative logL value for given phi.
###   These are only called by the optim function to find the MLE.
###
### - These functions are based on per gene list in the following order
###     genes -> amino acids -> synonymous codons -> counts
###   for reu13.list.g, y.g, and n.g.


### Return negative log likelihood of multinomial distribution.
### For ROC + NSEf model.
my.objectivePhiOne.nlogL.rocnsef <- function(phi, fitlist, reu13.list.g, y.g,
    n.g){
  ret <- 0
  for(i.aa in 1:length(reu13.list.g)){ # i'th amino acid
    for(i.scodon in 1:length(reu13.list.g[[i.aa]])){ # i'th synonymous codon
      Pos <- reu13.list.g[[i.aa]][[i.scodon]]
      if(length(Pos) > 0){ # avoid gene has no synonymous codon
        xm <- matrix(cbind(1, phi, phi * Pos), ncol = 3)
        # lp.vec <- my.inverse.mlogit(xm %*% fitlist[[i.aa]]$coef.mat,
        #                             log = TRUE)

        ### Since 1 * log(mu) + phi * Delta.t + phi * Pos * omega is the
        ### exponent term, and it can be scaled by x to yield
        ### (1/x) * M + 1 * S_1 + Pos * S_2
        ### where 1/x -> 0 as x -> infinity.
        exponent <- xm %*% fitlist[[i.aa]]$coef.mat
        id.infinite <- rowSums(!is.finite(exponent)) > 0
        if(any(id.infinite)){
          xm.tmp <- matrix(cbind(1, Pos[id.infinite]), ncol = 2)
          coef.tmp <- matrix(fitlist[[i.aa]]$coef.mat[-1,], nrow = 2)
          exponent[id.infinite,] <- xm.tmp %*% coef.tmp
        }
        lp.vec <- my.inverse.mlogit(exponent, log = TRUE)

        ### Only add the observed codon to the log lilelihood.
        ret <- ret + sum(lp.vec[, i.scodon])
      }
    }
  }
  -ret
} # End of my.objectivePhiOne.nlogL.rocnsef().

### For ROC model.
my.objectivePhiOne.nlogL.roc <- function(phi, fitlist, reu13.list.g, y.g, n.g){
  ret <- 0
  for(i.aa in 1:length(y.g)){ # i'th amino acid
    if(n.g[[i.aa]] > 0){ # avoid gene has no synonymous codon.
      xm <- matrix(cbind(1, phi), ncol = 2)
      # lp.vec <- my.inverse.mlogit(xm %*% fitlist[[i.aa]]$coef.mat, log = TRUE)

      ### Since 1 * log(mu) + phi * Delta.t is the exponent term, and
      ### it can be scaled by x to yield (1/x) * M + 1 * S_1
      ### where 1/x -> 0 as x -> infinity.
      exponent <- xm %*% fitlist[[i.aa]]$coef.mat
      id.infinite <- rowSums(!is.finite(exponent)) > 0
      if(any(id.infinite)){
        xm.tmp <- matrix(1, nrow = sum(id.infinite), ncol = 1)
        coef.tmp <- matrix(fitlist[[i.aa]]$coef.mat[-1,], nrow = 1)
        exponent[id.infinite,] <- xm.tmp %*% coef.tmp
      }
      lp.vec <- my.inverse.mlogit(exponent, log = TRUE)

      ### Only add the observed codon to the log lilelihood.
      # ret <- ret + sum((y.g[[i.aa]] * lp.vec))
      ### 0 * (-Inf) produces NaN
      ret <- ret + sum((y.g[[i.aa]] * lp.vec)[y.g[[i.aa]] != 0])
    }
  }

  -ret
} # End of my.objectivePhiOne.nlogL.roc().

### For NSEf model.
my.objectivePhiOne.nlogL.nsef <- function(phi, fitlist, reu13.list.g, y.g, n.g){
  ret <- 0
  for(i.aa in 1:length(reu13.list.g)){ # i'th amino acid
    for(i.scodon in 1:length(reu13.list.g[[i.aa]])){ # i'th synonymous codon
      Pos <- reu13.list.g[[i.aa]][[i.scodon]]
      if(length(Pos) > 0){ # avoid gene has no synonymous codon
        xm <- matrix(cbind(1, phi * Pos), ncol = 2)
        # lp.vec <- my.inverse.mlogit(xm %*% fitlist[[i.aa]]$coef.mat,
        #                             log = TRUE)

        ### Since 1 * log(mu) + phi * Pos * omega is the exponent term, and
        ### it can be scaled by x to yield (1/x) * M + Pos * S_2
        ### where 1/x -> 0 as x -> infinity.
        exponent <- xm %*% fitlist[[i.aa]]$coef.mat
        id.infinite <- rowSums(!is.finite(exponent)) > 0
        if(any(id.infinite)){
          xm.tmp <- matrix(Pos[id.infinite], ncol = 1)
          coef.tmp <- matrix(fitlist[[i.aa]]$coef.mat[-1,], nrow = 1)
          exponent[id.infinite,] <- xm.tmp %*% coef.tmp
        }
        lp.vec <- my.inverse.mlogit(exponent, log = TRUE)

        ### Only add the observed codon to the log lilelihood.
        ret <- ret + sum(lp.vec[, i.scodon])
      }
    }
  }
  -ret
} # End of my.objectivePhiOne.nlogL.nsef().
