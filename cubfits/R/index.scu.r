### This file is to compute the strength of translational selection
### on codon usage (SCU).
### Ref: Wallace, et.al. (2013) p. 1451.
###
### \pi_c^M = exp(M_c) / \sum{c' | a} exp(M_{c'})
###
### S'_c = S_c - \sum{c' | a} \pi_{c'}^M S_{c'}
###
### SCU_g = x_g \frac{1}{L_g} \sum_i S'_{c_i}
###
### mSCU_g = \frac{1}{L_g} \sum_i S'_{c_i}

calc_scu_values <- function(b, y.list, phi.Obs = NULL){
  aa.names <- names(b)

  S.c.new <- list()
  for(i.aa in 1:length(aa.names)){
    M.c <- b[[i.aa]]$coef.mat[1,]
    tmp <- exp(c(M.c, 0))
    pi.M.c <- tmp / sum(tmp)
    S.c <- c(b[[i.aa]]$coef.mat[2,], 0)
    S.c.new[[i.aa]] <- S.c - sum(pi.M.c * S.c)
  }

  mSCU <- NULL
  for(i.y in 1:length(y.list)){
    mSCU.g <- 0.0
    L.g <- 0
    for(i.aa in 1:length(aa.names)){
      mSCU.g <- mSCU.g + sum(y.list[[i.y]][[i.aa]] * S.c.new[[i.aa]])
      L.g <- L.g + sum(y.list[[i.y]][[i.aa]])
    }
    mSCU <- c(mSCU, mSCU.g / L.g)
  }
  names(mSCU) <- names(y.list)

  SCU <- NULL
  if(!is.null(phi.Obs)){
    SCU <- mSCU * phi.Obs
  }

  list(SCU = SCU, mSCU = mSCU)
} # End of calc_scu_values().

