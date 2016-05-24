"lcomCOPpv" <-
function(n, lcom, cop=NULL, para=NULL, repcoe=5E3, type="gno",
                  larsimn=1E4, larsimrep=15, uselarmu=FALSE, digits=5, ...) {

   type <- "gno"
   if(any("gno" == lmomco::dist.list())) {
      warning("the 'type' argument is providing a distribution not available ",
              "in the lmomco package, resetting to the 'gno'")
      type <- "gno"
   }
   if(uselarmu & larsimn == 0) {
       warning("argument 'uselarmu' is TRUE but 'larsimn == 0', incompatible, ",
               "resetting 'uselarmu' to FALSE")
       uselarmu <- FALSE
   }
   NAs <- rep(NA, 3)
   T2b <- T3b <- T4b <- data.frame(LcomType=NAs, N=NAs, Nrep=NAs,
                                   Mean=NAs, Lscale=NAs, Lskew=NAs, Lkurt=NAs)
   if(larsimn > 0) {
      if(larsimrep < 4) {
         warning("larsimrep is too small to compute first four L-moments, ",
                 "resetting to four")
         larsimrep <- 4
      }
      T2COPb <- vector(mode="numeric", larsimrep)
      T3COPb <- T4COPb <- T2COPb12 <- T3COPb12 <- T4COPb12 <- T2COPb
      T2COPb21 <- T3COPb21 <- T4COPb21 <- T2COPb
      message(" STATUS: Performing ", larsimrep,
              " replicate(s) of Tau_(234)[12:21] with ",
                                      larsimn,
              " sample size [large sam. sim.]")
      message("         Simulating replication ", appendLF=FALSE)
      for(i in 1:larsimrep) {
         message(i,"-", appendLF=FALSE)
         UV <- simCOP(n=larsimn, cop=cop, para=para, graphics=FALSE, ...)
         lmcr <- lmomco::lcomoms2(data.frame(U=UV$U, V=UV$V), nmom=4)
         T2COPb12[i] <- lmcr$T2[1,2];  T2COPb21[i] <- lmcr$T2[2,1]
         T3COPb12[i] <- lmcr$T3[1,2];  T3COPb21[i] <- lmcr$T3[2,1]
         T4COPb12[i] <- lmcr$T4[1,2];  T4COPb21[i] <- lmcr$T4[2,1]
         T2COPb[i]   <- (T2COPb12[i] + T2COPb21[i]) / 2
         T3COPb[i]   <- (T3COPb12[i] + T3COPb21[i]) / 2
         T4COPb[i]   <- (T4COPb12[i] + T4COPb21[i]) / 2
      }
      message("done")
      message("         Computing univariate L-moments-", appendLF=FALSE)
      lmrT2b12 <- lmomco::lmoms(T2COPb12, nmom=4)
      lmrT2b21 <- lmomco::lmoms(T2COPb21, nmom=4)
      lmrT3b12 <- lmomco::lmoms(T3COPb12, nmom=4)
      lmrT3b21 <- lmomco::lmoms(T3COPb21, nmom=4)
      lmrT4b12 <- lmomco::lmoms(T4COPb12, nmom=4)
      lmrT4b21 <- lmomco::lmoms(T4COPb21, nmom=4)
      lmrT2b   <- lmomco::lmoms(T2COPb,   nmom=4)
      lmrT3b   <- lmomco::lmoms(T3COPb,   nmom=4)
      lmrT4b   <- lmomco::lmoms(T3COPb,   nmom=4)
      T2b <- data.frame(LcomType=c("Tau2[12]", "Tau2[21]", "Tau2[12:21]"),
          N=rep(larsimn, 3), Nrep=rep(larsimrep, 3),
          Mean=  c(lmrT2b12$lambdas[1], lmrT2b21$lambdas[1], lmrT2b$lambdas[1]),
          Lscale=c(lmrT2b12$lambdas[2], lmrT2b21$lambdas[2], lmrT2b$lambdas[2]),
          Lskew= c(lmrT2b12$ratios[3],  lmrT2b21$ratios[3],  lmrT2b$ratios[3]),
          Lkurt= c(lmrT2b12$ratios[4],  lmrT2b21$ratios[4],  lmrT2b$ratios[4]))
      T3b <- data.frame(LcomType=c("Tau3[12]", "Tau3[21]", "Tau3[12|21]"),
          N=rep(larsimn, 3), Nrep=rep(larsimrep, 3),
          Mean=  c(lmrT3b12$lambdas[1], lmrT3b21$lambdas[1], lmrT3b$lambdas[1]),
          Lscale=c(lmrT3b12$lambdas[2], lmrT3b21$lambdas[2], lmrT3b$lambdas[2]),
          Lskew= c(lmrT3b12$ratios[3],  lmrT3b21$ratios[3],  lmrT3b$ratios[3]),
          Lkurt= c(lmrT3b12$ratios[4],  lmrT3b21$ratios[4],  lmrT3b$ratios[4]))
      T4b <- data.frame(LcomType=c("Tau4[12]", "Tau4[21]", "Tau4[12|21]"),
          N=rep(larsimn, 3), Nrep=rep(larsimrep, 3),
          Mean=  c(lmrT4b12$lambdas[1], lmrT4b21$lambdas[1], lmrT4b$lambdas[1]),
          Lscale=c(lmrT4b12$lambdas[2], lmrT4b21$lambdas[2], lmrT4b$lambdas[2]),
          Lskew= c(lmrT4b12$ratios[3],  lmrT4b21$ratios[3],  lmrT4b$ratios[3]),
          Lkurt= c(lmrT4b12$ratios[4],  lmrT4b21$ratios[4],  lmrT4b$ratios[4]))
      message("done")
   } else {
      message(" STATUS: Skipping large sample sim. because sample size is zero.")
   }

   T2table <- data.frame(LcomType=NAs, n=NAs, nrep=NAs, Mean=NAs, Lscale=NAs,
                                                        Lskew=NAs, Lkurt=NAs,
                                     sample.est=NAs, p.value=NAs, signif=NAs)
   T4table <- T3table <- T2table
   if(n > 0) {
      nrep <- as.integer(repcoe/sqrt(n))
      if(nrep == 0) nrep <- 1
      message(" STATUS: Performing ", nrep, " replicate(s) for ",n,
              " sample size [small sam. sim.]")
      message("           ", appendLF=FALSE)
      T2COP <- vector(mode="numeric", length=nrep)
      T3COP <- T4COP <- T2COP12 <- T3COP12 <- T4COP12 <- T2COP
      T2COP21 <- T3COP21 <- T4COP21 <- T2COP
      for(i in 1:nrep) {
         if(length(grep('00+$', i, perl=TRUE)) == 1) message(i,"-",appendLF=FALSE)
         UV <- simCOP(n=n, para=para, cop=cop, graphics=FALSE)
         lmcr <- lmomco::lcomoms2(data.frame(U=UV$U, V=UV$V), nmom=4)
         T2COP12[i] <- lmcr$T2[1,2]; T2COP21[i] <- lmcr$T2[2,1]
         T3COP12[i] <- lmcr$T3[1,2]; T3COP21[i] <- lmcr$T3[2,1]
         T4COP12[i] <- lmcr$T4[1,2]; T4COP21[i] <- lmcr$T4[2,1]
           T2COP[i] <- (T2COP12[i] + T2COP21[i]) / 2
           T3COP[i] <- (T3COP12[i] + T3COP21[i]) / 2
           T4COP[i] <- (T4COP12[i] + T4COP21[i]) / 2
      }
      message("done")

      message("         Computing univariate L-moments for L-comoments 'r[12]'")
      # Work on [12]
      lmrT2 <- lmomco::lmoms(T2COP12); lmrT2t <- lmrT2
      lmrT3 <- lmomco::lmoms(T3COP12); lmrT3t <- lmrT3
      lmrT4 <- lmomco::lmoms(T4COP12); lmrT4t <- lmrT4
      if(uselarmu) {
         message("         Using mean of large simulation and rescaling ",
                 "L-scale by small simulation LCV")
         lmrT2t$lambdas[1] <- T2b[1,4]
         lmrT3t$lambdas[1] <- T3b[1,4]
         lmrT4t$lambdas[1] <- T4b[1,4]
         lmrT2t$lambdas[2] <- lmrT2$ratios[2]*T2b[1,4]
         lmrT3t$lambdas[2] <- lmrT3$ratios[2]*T3b[1,4]
         lmrT4t$lambdas[2] <- lmrT4$ratios[2]*T4b[1,4]
         lmrT2t$ratios[2] <- lmrT2$ratios[2]
         lmrT3t$ratios[2] <- lmrT3$ratios[2]
         lmrT4t$ratios[2] <- lmrT4$ratios[2]
      }
      par2      <- lmomco::lmom2par(lmrT2t, type=type);# print(par2)
      par3      <- lmomco::lmom2par(lmrT3t, type=type);# print(par3)
      par4      <- lmomco::lmom2par(lmrT4t, type=type);# print(par4)
      T2p.value <- lmomco::plmomco(lcom$T2[1,2], par2);# print(T2p.value)
      T3p.value <- lmomco::plmomco(lcom$T3[1,2], par3);# message("PARA",T3p.value)
      T4p.value <- lmomco::plmomco(lcom$T4[1,2], par4);# print(T4p.value)
      if(T2p.value > 0.5) T2p.value <- 1 - T2p.value
      if(T3p.value > 0.5) T3p.value <- 1 - T3p.value
      if(T4p.value > 0.5) T4p.value <- 1 - T4p.value
      T2table[1,] <- c("Tau2[12]", n, nrep, lmrT2$lambdas[1], lmrT2$lambdas[2],
                                lmrT2$ratios[3:4], lcom$T2[1,2], T2p.value, NA)
      T3table[1,] <- c("Tau3[12]", n, nrep, lmrT3$lambdas[1], lmrT3$lambdas[2],
                                lmrT3$ratios[3:4], lcom$T3[1,2], T3p.value, NA)
      T4table[1,] <- c("Tau4[12]", n, nrep, lmrT4$lambdas[1], lmrT4$lambdas[2],
                                lmrT4$ratios[3:4], lcom$T4[1,2], T4p.value, NA)

      message("         Computing univariate L-moments for L-comoments 'r[21]'")
      # Work on [21]
      lmrT2 <- lmomco::lmoms(T2COP21); lmrT2t <- lmrT2
      lmrT3 <- lmomco::lmoms(T3COP21); lmrT3t <- lmrT3
      lmrT4 <- lmomco::lmoms(T4COP21); lmrT4t <- lmrT4
      if(uselarmu) {
         message("         Using mean of large simulation and rescaling ",
                 "L-scale by small simulation LCV")
         lmrT2t$lambdas[1] <- T2b[2,4]
         lmrT3t$lambdas[1] <- T3b[2,4]
         lmrT4t$lambdas[1] <- T4b[2,4]
         lmrT2t$lambdas[2] <- lmrT2$ratios[2]*T2b[2,4]
         lmrT3t$lambdas[2] <- lmrT3$ratios[2]*T3b[2,4]
         lmrT4t$lambdas[2] <- lmrT4$ratios[2]*T4b[2,4]
         lmrT2t$ratios[2] <- lmrT2$ratios[2]
         lmrT3t$ratios[2] <- lmrT3$ratios[2]
         lmrT4t$ratios[2] <- lmrT4$ratios[2]
      }
      par2        <- lmomco::lmom2par(lmrT2t, type=type)
      par3        <- lmomco::lmom2par(lmrT3t, type=type)
      par4        <- lmomco::lmom2par(lmrT4t, type=type)
      T2p.value   <- lmomco::plmomco(lcom$T2[2,1], par2)
      T3p.value   <- lmomco::plmomco(lcom$T3[2,1], par3)
      T4p.value   <- lmomco::plmomco(lcom$T4[2,1], par4)
      if(T2p.value > 0.5) T2p.value <- 1 - T2p.value
      if(T3p.value > 0.5) T3p.value <- 1 - T3p.value
      if(T4p.value > 0.5) T4p.value <- 1 - T4p.value
      T2table[2,] <- c("Tau2[21]", n, nrep, lmrT2$lambdas[1], lmrT2$lambdas[2],
                               lmrT2$ratios[3:4], lcom$T2[2,1], T2p.value, NA)
      T3table[2,] <- c("Tau3[21]", n, nrep, lmrT3$lambdas[1], lmrT3$lambdas[2],
                               lmrT3$ratios[3:4], lcom$T3[2,1], T3p.value, NA)
      T4table[2,] <- c("Tau4[21]", n, nrep, lmrT4$lambdas[1], lmrT4$lambdas[2],
                               lmrT4$ratios[3:4], lcom$T4[2,1], T4p.value, NA)

      message("         Computing univariate L-moments for L-comoments 'r[12|21]'")
      # Working on [12|21]
      T2 <- (lcom$T2[1,2] + lcom$T2[2,1]) / 2
      T3 <- (lcom$T3[1,2] + lcom$T3[2,1]) / 2
      T4 <- (lcom$T4[1,2] + lcom$T4[2,1]) / 2
      lmrT2 <- lmomco::lmoms(T2COP); lmrT2t <- lmrT2 # Create temporary copies that we can do
      lmrT3 <- lmomco::lmoms(T3COP); lmrT3t <- lmrT3 # substitution on the mean and L-scale dependng
      lmrT4 <- lmomco::lmoms(T4COP); lmrT4t <- lmrT4 # on setting of 'uselarmu' for the p-value computation.
      if(uselarmu) {
         message("         Using mean of large simulation and rescaling ",
                 "L-scale by small simulation LCV")
         lmrT2t$lambdas[1] <- T2b[3,4]
         lmrT3t$lambdas[1] <- T3b[3,4]
         lmrT4t$lambdas[1] <- T4b[3,4]
         lmrT2t$lambdas[2] <- lmrT2$ratios[2]*T2b[3,4]
         lmrT3t$lambdas[2] <- lmrT3$ratios[2]*T3b[3,4]
         lmrT4t$lambdas[2] <- lmrT4$ratios[2]*T4b[3,4]
         lmrT2t$ratios[2] <- lmrT2$ratios[2]
         lmrT3t$ratios[2] <- lmrT3$ratios[2]
         lmrT4t$ratios[2] <- lmrT4$ratios[2]
      }
      par2        <- lmomco::lmom2par(lmrT2t, type=type)
      par3        <- lmomco::lmom2par(lmrT3t, type=type)
      par4        <- lmomco::lmom2par(lmrT4t, type=type)
      T2p.value   <- lmomco::plmomco(T2, par2)
      T3p.value   <- lmomco::plmomco(T3, par3)
      T4p.value   <- lmomco::plmomco(T4, par4)
      if(T2p.value > 0.5) T2p.value <- 1 - T2p.value
      if(T3p.value > 0.5) T3p.value <- 1 - T3p.value
      if(T4p.value > 0.5) T4p.value <- 1 - T4p.value
      T2table[3,] <- c("Tau2[12:21]", n, nrep,
       lmrT2$lambdas[1], lmrT2$lambdas[2], lmrT2$ratios[3:4], T2, T2p.value, NA)
      T3table[3,] <- c("Tau3[12:21]", n, nrep,
       lmrT3$lambdas[1], lmrT3$lambdas[2], lmrT3$ratios[3:4], T3, T3p.value, NA)
      T4table[3,] <- c("Tau4[12:21]", n, nrep,
       lmrT4$lambdas[1], lmrT4$lambdas[2], lmrT4$ratios[3:4], T4, T4p.value, NA)
      for(i in 2:9) { # annoying part of my logic flow, having to recast into numerics
         T2table[,i] <- as.numeric(T2table[,i])
         T3table[,i] <- as.numeric(T3table[,i])
         T4table[,i] <- as.numeric(T4table[,i])
      }
      sig.codes <- c("-", ".", "*", "**", "***")
      sig.val  <- c(2, 0.2, 0.1, 0.02, 0.002) # Field [1] is purposelfully present and not used
      sig.val <- sig.val/2
      sig.valc <- 1 - sig.val
      for(i in 1:3) {
         T2table[i,10] <- ifelse(T2table[i,9] <=  sig.val[5]  |
                                 T2table[i,9] >= sig.valc[5],  sig.codes[5],
                          ifelse(T2table[i,9] <=  sig.val[4]  |
                                 T2table[i,9] >= sig.valc[4],  sig.codes[4],
                          ifelse(T2table[i,9] <=  sig.val[3]  |
                                 T2table[i,9] >= sig.valc[3],  sig.codes[3],
                          ifelse(T2table[i,9] <=  sig.val[2]  |
                                 T2table[i,9] >= sig.valc[2],  sig.codes[2],
                                                               sig.codes[1]))))
         T3table[i,10] <- ifelse(T3table[i,9] <=  sig.val[5]  |
                                 T3table[i,9] >= sig.valc[5],  sig.codes[5],
                          ifelse(T3table[i,9] <=  sig.val[4]  |
                                 T3table[i,9] >= sig.valc[4],  sig.codes[4],
                          ifelse(T3table[i,9] <=  sig.val[3]  |
                                 T3table[i,9] >= sig.valc[3],  sig.codes[3],
                          ifelse(T3table[i,9] <=  sig.val[2]  |
                                 T3table[i,9] >= sig.valc[2],  sig.codes[2],
                                                               sig.codes[1]))))
         T4table[i,10] <- ifelse(T4table[i,9] <=  sig.val[5]  |
                                 T4table[i,9] >= sig.valc[5],  sig.codes[5],
                          ifelse(T4table[i,9] <=  sig.val[4]  |
                                 T4table[i,9] >= sig.valc[4],  sig.codes[4],
                          ifelse(T4table[i,9] <=  sig.val[3]  |
                                 T4table[i,9] >= sig.valc[3],  sig.codes[3],
                          ifelse(T4table[i,9] <=  sig.val[2]  |
                                 T4table[i,9] >= sig.valc[2],  sig.codes[2],
                                                               sig.codes[1]))))
      }
   } else {
      message(" STATUS: Skipping large sample sim. because sample size is zero.")
   }
   message("STATUS: Rounding digits")
   if(! is.na(digits)) {
      if(larsimn > 0) {
          T2b[,4:7]    <- round(T2b[,4:7],     digits=digits)
          T3b[,4:7]    <- round(T3b[,4:7],     digits=digits)
          T4b[,4:7]    <- round(T4b[,4:7],     digits=digits)
      }
      if(n > 0) {
         T2table[,4:9] <- round(T2table[,4:9], digits=digits)
         T3table[,4:9] <- round(T3table[,4:9], digits=digits)
         T4table[,4:9] <- round(T4table[,4:9], digits=digits)
      }
   }
   txt <- ifelse(uselarmu,
    "Large sample size simulated means in 'Ntable' were used for p.value computation",
    "Small sample size simulated means in 'ntable' were used for p.value computation")
   zz <- list(text="Sample L-comoment Characterization of a Copula",
              p.value.how=txt,
              Ntable=list(Lcocorr=T2b,     Lcoskew=T3b,     Lcokurt=T4b),
              ntable=list(Lcocorr=T2table, Lcoskew=T3table, Lcokurt=T4table))
   return(zz)
}

