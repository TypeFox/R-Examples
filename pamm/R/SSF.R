`SSF` <-
  function (numsim, tss, nbstep = 10, randompart, fixed = c(0, 
                                                            1, 0), n.X = NA, autocorr.X = 0, X.dist="gaussian", intercept = 0, 
            exgr = NA, exrepl = NA, heteroscedasticity=c("null")) 
  {
    o.warn <- getOption("warn")

    if (is.na(exgr)[[1]]) {
      grmin = 2
      grmax = tss/2
    }
    if (is.na(exrepl)[[1]]) {
      remin = 2
      remax = tss/2
    }
    if (!is.na(exgr)[[1]]) {
      grmin = exgr[[1]]
      grmax = exgr[[2]]
    }
    if (!is.na(exrepl)[[1]]) {
      remin = exrepl[[1]]
      remax = exrepl[[2]]
    }
    group = round(seq(grmin, grmax, I((grmax - grmin)/nbstep)))
    repl = ceiling(tss/group)
    mg.r0 = unique(matrix(c(group, repl), ncol = 2))
    mg.r1 = subset(mg.r0, mg.r0[, 2] >= remin)
    mg.r = subset(mg.r1, mg.r1[, 2] <= remax)
    stepvec = c(1:length(mg.r[, 1]))
    VI <- as.numeric(randompart[[1]])
    VS <- as.numeric(randompart[[2]])
    VR <- as.numeric(randompart[[3]])
    if (length(randompart <= 4)) {
      CorIS <- as.numeric(randompart[[4]])
      CovIS <- CorIS * sqrt(VI) * sqrt(VS)
    }
    else {
      if (randompart[[5]] == "cor") {
        CorIS <- as.numeric(randompart[[4]])
        CovIS <- CorIS * sqrt(VI) * sqrt(VS)
      }
      if (randompart[[5]] == "cov") {
        CovIS <- as.numeric(randompart[[4]])
      }
    }
    sigma <- matrix(c(VI, CovIS, CovIS, VS), ncol = 2)
    
    Hetero <- heteroscedasticity[[1]]
    het <- as.numeric(heteroscedasticity[-1])
    
    if (X.dist=="gaussian") {
      FM <- fixed[[1]]
      FV <- fixed[[2]]
      FE <- fixed[[3]]
    }
    if (X.dist=="unif") {
      Xmin <- fixed[[1]]
      Xmax <- fixed[[2]]
      FE <- fixed[[3]]
    }
    
    iD <- numeric(length(mg.r[, 1]))
    rp <- numeric(length(mg.r[, 1]))
    ss <- numeric(length(mg.r[, 1]))
    powersl <- numeric(numsim)
    pvalsl <- numeric(numsim)
    slpowestimate <- numeric(length(mg.r[, 1]))
    slpowCIlower <- numeric(length(mg.r[, 1]))
    slpowCIupper <- numeric(length(mg.r[, 1]))
    slpvalestimate <- numeric(length(mg.r[, 1]))
    slpvalCIlower <- numeric(length(mg.r[, 1]))
    slpvalCIupper <- numeric(length(mg.r[, 1]))
    powerint <- numeric(numsim)
    pvalint <- numeric(numsim)
    intpowestimate <- numeric(length(mg.r[, 1]))
    intpowCIlower <- numeric(length(mg.r[, 1]))
    intpowCIupper <- numeric(length(mg.r[, 1]))
    intpvalestimate <- numeric(length(mg.r[, 1]))
    intpvalCIlower <- numeric(length(mg.r[, 1]))
    intpvalCIupper <- numeric(length(mg.r[, 1]))
    nsim.used.sl <- numeric(length(mg.r[, 1]))
    nsim.used.int <- numeric(length(mg.r[, 1]))
    
    kk <- 0
    for (k in stepvec) {
      N <- tss
      n.x <- ifelse( is.na(n.X)==TRUE,N, n.X)
      
      for (i in 1:numsim) {
        options(warn=2)
        if (X.dist=="gaussian"){
          if (autocorr.X==0) { ef <- rnorm(n.x, FM, sqrt(FV)) }
          else {
            y <- numeric(n.x)
            phi <- autocorr.X
            y[1] <- rnorm(1, 0, sd = sqrt(FV))
            for (t in 2:n.x) { y[t] <- rnorm(1, y[t-1]*phi, sd = sqrt(FV)) }
            ef <- y+FM
          }
        }
        
        if (X.dist=="unif"){
          if (autocorr.X==0) { ef <- runif(n.x, Xmin, Xmax) }
          else { stop("autocorrelation in fixed effects is not yet implemented for uniform distribution") }
        }
        
        if (n.x!=N) {
          if (n.x>=mg.r[k, 2]) {
            inief <- sample(1:(n.x-mg.r[k, 2]+1),mg.r[k, 1],replace=TRUE)
            EFrk <- rep(inief,mg.r[k, 2]) + rep (0:(mg.r[k, 2]-1),each=mg.r[k, 1])  #EFrk <- rep(inief,each=r) + rep (0:(r-1),k)
            EF <- ef[EFrk][1:N]
          }
          if (n.x<mg.r[k, 2]) {
            EF <- numeric(N)
            EF[1:(n.x*mg.r[k, 1])] <- rep(ef,each=mg.r[k, 1])
            EF[(n.x*mg.r[k, 1]+1):N] <- sample(ef,length((n.x*mg.r[k, 1]+1):N),replace=TRUE)
          }
        }
        else { EF <- ef }
        
        er <- numeric(length(N))
        if (Hetero=="null") (er <- rnorm(N, intercept, sqrt(VR)))
        if (Hetero=="power") (
          for (n in 1:N) {er[n] <- rnorm(1, intercept, sqrt(VR*(het[1]+abs(EF[n])^het[2])^2))} )
        if (Hetero=="exp")  (
          for (n in 1:N) {er[n] <- rnorm(1, intercept, sqrt(VR*exp(2*het[1]*EF[n])))} )
        
        db <- data.frame(ID = rep(1:mg.r[k, 1], mg.r[k, 2])[1:N], 
                         obs = 1:N, error = er, EF = EF)
        x <- rmvnorm(mg.r[k, 1], c(0, 0), sigma, method = "svd")
        db$rand.int <- rep(x[, 1], mg.r[k, 2])[1:N]
        db$rand.sl <- rep(x[, 2], mg.r[k, 2])[1:N]
        db$Y <- db$rand.int + (db$rand.sl + FE) * db$EF + 
          db$error
        
        m1.lmer <- try(lmer(Y ~ EF + (1 | ID), data = db),TRUE)
        if (class(m1.lmer)!="merModLmerTest")  {          	
          powerint[i] <- NA
          pvalint[i] <- NA
        }
        else{           
          lrt1 <- rand(m1.lmer)
          pvint <- lrt1[[1]][1,3]
          powerint[i] <- pvint <= 0.05
          pvalint[i] <- pvint
        }
        
        m2.lmer <- try(lmer(Y ~ EF + (EF | ID), data = db),TRUE)
        if (class(m2.lmer)!="merModLmerTest" || class(m1.lmer)!="merModLmerTest")  {          	
          powersl[i] <- NA
          pvalsl[i] <- NA
        }
        else{ 
          anosl <- anova(m2.lmer, m1.lmer, refit=FALSE)
          powersl[i] <- anosl[2, "Pr(>Chisq)"] <= 0.05
          pvalsl[i] <- anosl[2, "Pr(>Chisq)"]
        }
      }
      ## add number of models dropped for estimates
      options(warn=o.warn)
      kk <- kk + 1
      iD[kk] <- mg.r[k, 1]
      rp[kk] <- round(N/mg.r[k, 1], digits = 2)
      ss[kk] <- N
      slCIpow <- ci(powersl, na.rm=TRUE)
      slpowestimate[kk] <- slCIpow["Estimate"]
      slpowCIlower[kk] <- slCIpow["CI lower"]
      slpowCIupper[kk] <- slCIpow["CI upper"]
      slCIpval <- ci(pvalsl, na.rm=TRUE)
      slpvalestimate[kk] <- slCIpval["Estimate"]
      slpvalCIlower[kk] <- slCIpval["CI lower"]
      slpvalCIupper[kk] <- slCIpval["CI upper"]
      intCIpow <- ci(powerint, na.rm=TRUE)
      intpowestimate[kk] <- intCIpow["Estimate"]
      intpowCIlower[kk] <- intCIpow["CI lower"]
      intpowCIupper[kk] <- intCIpow["CI upper"]
      intCIpval <- ci(pvalint, na.rm=TRUE)
      intpvalestimate[kk] <- intCIpval["Estimate"]
      intpvalCIlower[kk] <- intCIpval["CI lower"]
      intpvalCIupper[kk] <- intCIpval["CI upper"]
      nsim.used.sl[kk] <- numsim - sum(is.na(pvalsl))
      nsim.used.int[kk] <- numsim - sum(is.na(pvalint))
    }
    sim.sum <- data.frame(nb.ID = iD, nb.repl = rp, N = ss, int.pval = intpvalestimate, 
                          CIlow.ipv = intpvalCIlower, CIup.ipv = intpvalCIupper, 
                          int.power = intpowestimate, CIlow.ipo = intpowCIlower, 
                          CIup.ipo = intpowCIupper, sl.pval = slpvalestimate, CIlow.slpv = slpvalCIlower, 
                          CIup.slpv = slpvalCIupper, sl.power = slpowestimate, 
                          CIlow.slpo = slpowCIlower, CIup.slpo = slpowCIupper,
                          nsim.int = nsim.used.int, nsim.sl = nsim.used.sl)

    class(sim.sum) <- c("SSF", "data.frame")
    sim.sum
  }
