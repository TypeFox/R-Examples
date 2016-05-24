### Wrapper of find.prop.model.nse().
prop.model.nse <- function(b.Init, reu13.df, phi.Obs.lim = c(0.01, 10), phi.Obs.scale = 1,
                           nclass = 40, x.log10 = TRUE, delta_a12=0, a_2=1, positions=c(100, 200, 400, 800))
  {
  aa.names <- names(b.Init)
  if("Z" %in% aa.names){
    synonymous.codon <- .CF.GV$synonymous.codon.split[aa.names]
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon[aa.names]
  }

  ### Compute mutation and elong
  for(i.aa in aa.names){
    b.Init[[i.aa]]$u.codon <- synonymous.codon[[i.aa]]
  }

  ### For SCUO w/o phi.Obs simulations, the prior fixes the mean of phi.Obs at
  ### 1, but lets std drift only. Therefore, E[Phi] may or may not be
  ### distributed in the same range of phi.Obs.
  ### Convert phi.Obs.lim to SCUO's scale, and evaluate the proportional
  ### frequences of synonymous codons, then convert back to phi.Obs's
  ### scale (for plotting) inside the function find.prop.model.nse().
  phi.Obs.lim <- phi.Obs.lim / phi.Obs.scale

  ### For better plottling
  if(x.log10){
    phi.bin <- 10^seq(log10(phi.Obs.lim[1]), log10(phi.Obs.lim[2]), length = nclass)
  } else{
    phi.bin <- seq(phi.Obs.lim[1], phi.Obs.lim[2], length = nclass)
  }

  ### Call find.prop.model.nse().
  ret <- lapply(1:length(positions), function(ii) {
    find.prop.model.nse(b.Init, reu13.df, phi.bin, phi.Obs.scale, delta_a12=delta_a12, a_2=a_2, positions[ii])
  })
  ret
} # End of prop.model.nse().


### Summarize by amino acid for bInint could be from MCMC outputs.
find.prop.model.nse <- function(b.Init, reu13.df, phi.bin, phi.Obs.scale = 1,
                                delta_a12=0, a_2=1, pos=0)
  {
  u.aa <- unique(names(b.Init))
  x <- cbind(1, phi.bin, phi.bin)

  predict.nse <- list()
  for(aa in u.aa){
    #b.Init[[aa]]$coef.mat <- -b.Init[[aa]]$coef.mat
    b.Init[[aa]]$coef.mat <- matrix(rbind(b.Init[[aa]]$coef.mat, b.Init[[aa]]$coef.mat[2,]), nrow=3)
    b.Init[[aa]]$coef.mat[1,] <- b.Init[[aa]]$coef.mat[1,] * -1
    b.Init[[aa]]$coef.mat[2,] <- b.Init[[aa]]$coef.mat[2,] * delta_a12
    b.Init[[aa]]$coef.mat[3,] <- b.Init[[aa]]$coef.mat[3,] * pos * a_2
    #exponent <- x %*% b.Init[[aa]]$coef.mat
    exponent <- x %*% (-b.Init[[aa]]$coef.mat)
    scodon.prob <- my.inverse.mlogit(exponent)
    predict.nse[[aa]] <- cbind(scodon.prob, phi.bin * phi.Obs.scale)
    predict.nse[[aa]] <- as.data.frame(predict.nse[[aa]],
                                       stringsAsFactors = FALSE)

    u.codon <- b.Init[[aa]]$u.codon
    colnames(predict.nse[[aa]]) <- c(u.codon, "center")

    ### Add attr for u.codon.star
    if(all(b.Init[[aa]]$coef.mat[2,] < 0)){
      u.codon[length(u.codon)] <- paste(u.codon[length(u.codon)], "*", sep = "")
    } else{
      id.max <- which.max(b.Init[[aa]]$coef.mat[2,])
      u.codon.max <- colnames(b.Init[[aa]]$coef.mat)[id.max]
      u.codon[u.codon == u.codon.max] <- paste(u.codon.max, "*")
    }
    attr(predict.nse[[aa]], "u.codon.star") <- u.codon
  }

  names(predict.nse) <- u.aa
  predict.nse
} # End of find.prop.model.nse().
