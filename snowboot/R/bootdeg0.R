bootdeg0 <- function(sam.out, num.sam, n.boot) {
      # This function obtains the bootstrap samples for each sample from a network sam.out is the output of Oempdegreedistrib
      # num.sam is the number of different samples taken from the same network (Scalar or vector) n.boot is the number of
      # bootstrap samples taken from each sample
      n.seeds <- sam.out$n.seeds
      n.neigh <- sam.out$n.neigh
      if (length(num.sam) == 1)
            num.sam <- 1:num.sam
      empd <- as.list(rep(NA, length(num.sam)))
      i <- 1
      for (m in num.sam) {
            val.seeds <- sam.out$val.seeds[[m]]
            freq.deg.seeds <- sam.out$samples[[m]]$freq.deg.seeds
            bsam.seeds <- myBsample(val.seeds, n.seeds, n.boot, prob = freq.deg.seeds)
            values <- sam.out$values[[m]]  #all the possible degree values toresample
            #### Frequency ##### (Not the relative frequency)
            Fseed <- t(apply(bsam.seeds, 1, table.row, vect = values))  #freq (sorted according to values)
            # browser()
            empd.seeds <- Fseed/n.seeds
            empd[[i]] <- list(empd.seeds = empd.seeds)
            i <- i + 1
      }  # for(m in num.sam)
      list(values = sam.out$values[num.sam], empd = empd, num.sam = num.sam,
           n.boot = n.boot, n.neigh = n.neigh, seeds1 = sam.out$seeds1[num.sam,])
}
