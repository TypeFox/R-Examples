OparametersEst <- function(outOempd) {
      # This function obtains the parameters based on the empirical distritributions derived from the samples outOempd is
      # output of Oempdegreedistrib
      num.sam <- outOempd$num.sam
      mean <- rep(NA, num.sam)
      quartiles <- array(NA, c(num.sam, 3))  #1rst,2nd and 3rd quartile
      rfreq <- array(NA, c(num.sam, 5))  #freq of 0,1,2,3,4
      deciles <- array(NA, c(num.sam, 9))  #10,20,30,40,50,60,70,80,90
      for (m in 1:outOempd$num.sam) {
            distrib <- outOempd$Oempd[[m]]$Oempd
            # vals<-as.numeric(names(distrib)) #checar diversos formatos
            vals <- outOempd$values[[m]]
            mean[m] <- sum(vals * distrib)
            quartiles[m, ] <- vals[sapply(X = c(0.25, 0.5, 0.75), FUN = min.greater, v = cumsum(distrib), ge = TRUE)]
            rfreq[m, ] <- c(distribvals(distrib, vals, 0), distribvals(distrib, vals, 1), distribvals(distrib, vals, 2), distribvals(distrib,
                                                                                                                                     vals, 3), distribvals(distrib, vals, 4))
            deciles[m, ] <- vals[sapply(X = seq(0.1, 0.9, by = 0.1), FUN = min.greater, v = cumsum(distrib), ge = TRUE)]
      }
      # browser()
      list(mean = mean, quartiles = quartiles, rfreq = rfreq, deciles = deciles, m = m)
}
