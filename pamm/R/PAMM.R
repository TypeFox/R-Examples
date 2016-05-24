`PAMM` <-
  function (numsim, group, repl, randompart, fixed = c(0, 1, 0), n.X=NA, autocorr.X=0,
            X.dist="gaussian", intercept=0,heteroscedasticity=c("null"),ftype="lmer",
            mer.sim=FALSE)
  {
    o.warn <- getOption("warn")
    
    VI <- randompart[1]
    VS <- randompart[2]
    VR <- randompart[3]
    if (length(randompart) == 5) {
      if (randompart[5] == "cor") {
        CorIS <- randompart[4]
        CovIS <- CorIS * sqrt(VI) * sqrt(VS)
      }
      else if (randompart[5] == "cov") {
        CovIS <-randompart[4]
      }
      
    }
    else {
      CorIS <- randompart[4]
      CovIS <- CorIS * sqrt(VI) * sqrt(VS)
    }
    M <- matrix(c(VI, CovIS, CovIS, VS), ncol = 2)
    
    Hetero <- heteroscedasticity[1]
    het <- as.numeric(heteroscedasticity[-1])
    
    if (X.dist=="gaussian") {
      FM <- fixed[1]
      FV <- fixed[2]
      FE <- fixed[3]
    }
    if (X.dist=="unif") {
      Xmin <- fixed[1]
      Xmax <- fixed[2]
      FE <- fixed[3]
    }
    
    iD <- numeric(length(repl) * length(group))
    rp <- numeric(length(repl) * length(group))
    powersl <- numeric(numsim)
    pvalsl <- numeric(numsim)
    slpowestimate <- numeric(length(repl) * length(group))
    slpowCIlower <- numeric(length(repl) * length(group))
    slpowCIupper <- numeric(length(repl) * length(group))
    slpvalestimate <- numeric(length(repl) * length(group))
    slpvalCIlower <- numeric(length(repl) * length(group))
    slpvalCIupper <- numeric(length(repl) * length(group))
    powerint <- numeric(numsim)
    pvalint <- numeric(numsim)
    intpowestimate <- numeric(length(repl) * length(group))
    intpowCIlower <- numeric(length(repl) * length(group))
    intpowCIupper <- numeric(length(repl) * length(group))
    intpvalestimate <- numeric(length(repl) * length(group))
    intpvalCIlower <- numeric(length(repl) * length(group))
    intpvalCIupper <- numeric(length(repl) * length(group))
    nsim.used.sl <- numeric(length(repl) * length(group))
    nsim.used.int <- numeric(length(repl) * length(group))
    
    kk <- 0
    for (k in group) {
      for (r in repl) {
        N <- k * r
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
            if (n.x>=r) {
              inief <- sample(1:(n.x-r+1),k,replace=TRUE)
              EFrk <- rep(inief,r) + rep (0:(r-1),each=k)  #EFrk <- rep(inief,each=r) + rep (0:(r-1),k)
              EF <- ef[EFrk]
            }
            if (n.x<r) {
              EF <- numeric(N)
              EF[1:(n.x*k)] <- rep(ef,each=k)
              EF[(n.x*k+1):N] <- sample(ef,length((n.x*k+1):N),replace=TRUE)
            }
          }
          else { EF <- ef }
          
          db <- data.frame(ID = rep(1:k, r), obs = 1:N, EF = EF)
          if (mer.sim == TRUE) {
            family <- gaussian
            sigma <- sqrt(VR)
            beta <- c(intercept,fixed[3])
            names(beta) <- c("(Intercept)","EF")
            theta <- as.vector(chol(M)/sigma)[c(1,3,4)]
            names(theta) <- c("ID.(Intercept)","ID.EF.(Intercept)","ID.EF")
            params <- list(beta=beta, theta= theta, sigma=sigma)
            y <- simulate(formula(~EF + (EF | ID)),newdata=db,family=family, newparams=params)
            db$Y <- y[,1]
          }
          else {
            er <- numeric(length(N))
            if (Hetero=="null") (er <- rnorm(N, intercept, sqrt(VR)))
            if (Hetero=="power") (
              for (n in 1:N) {er[n] <- rnorm(1, intercept, sqrt(VR*(het[1]+abs(EF[n])^het[2])^2))} )
            if (Hetero=="exp")  (
              for (n in 1:N) {er[n] <- rnorm(1, intercept, sqrt(VR*exp(2*het[1]*EF[n])))} )
            db$error <- er
            x <- rmvnorm(k, c(0, 0), M, method = "svd")
            db$rand.int <- rep(x[, 1], r)
            db$rand.sl <- rep(x[, 2], r)
            db$Y <- db$rand.int + (db$rand.sl + FE) * db$EF + db$error		
          }
          #                if (ftype=="lme") {
          #                m.lm <- lm(Y ~ EF, data = db)
          #                     m1.lme <- lme(Y ~ EF,random= ~1 | ID,weights=varConstPower(form=~EF), data = db) 
          #                     pvint <- pchisq(-2 * (logLik(m.lm, REML = TRUE) - 
          #                     logLik(m1.lme, REML = TRUE))[[1]], 1, lower.tail = FALSE)
          #                     powerint[i] <- pvint <= 0.05
          #                     pvalint[i] <- pvint
          #                     m2.lme <- lme(Y ~ EF,random= ~EF | ID,weights=varConstPower(form=~EF), data = db) 
          #                     anosl <- anova(m2.lme, m1.lme)
          #                     powersl[i] <- anosl[2, "Pr(>Chisq)"] <= 0.05
          #                     pvalsl[i] <- anosl[2, "Pr(>Chisq)"]               
          #                }
          #                else {
          
          m1.lmer <- try(lmer(Y ~ EF + (1 | ID), data = db), silent = TRUE )
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
          if (VS > 0) {
            #                     	m2.lmer1 <- lmer(Y ~ EF + (1|ID) + (0 + EF|ID), data = db)                     
            m2.lmer2 <- try(lmer(Y ~ EF + (EF | ID), data = db), TRUE) 
            if (class(m2.lmer2)!="merModLmerTest" || class(m1.lmer)!="merModLmerTest")  {                       
              powerint[i] <- NA
              pvalint[i] <- NA
            }
            else{
              anosl <- anova(m2.lmer2, m1.lmer,refit=FALSE)
              powersl[i] <- anosl[2, "Pr(>Chisq)"] <= 0.05
              pvalsl[i] <- anosl[2, "Pr(>Chisq)"]
            }
          }
          #                 }
        }
        options(warn=o.warn)
        kk <- kk + 1
        iD[kk] <- k
        rp[kk] <- r
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
    }
    if (VS > 0) {
      sim.sum <- data.frame(nb.ID = iD, nb.repl = rp, int.pval = intpvalestimate, 
                            CIlow.ipv = intpvalCIlower, CIup.ipv = intpvalCIupper, 
                            int.power = intpowestimate, CIlow.ipo = intpowCIlower, 
                            CIup.ipo = intpowCIupper, sl.pval = slpvalestimate, CIlow.slpv = slpvalCIlower, 
                            CIup.slpv = slpvalCIupper, sl.power = slpowestimate, 
                            CIlow.slpo = slpowCIlower, CIup.slpo = slpowCIupper,
                            nsim.int = nsim.used.int, nsim.sl = nsim.used.sl)
    }
    else {
      sim.sum <- data.frame(nb.ID = iD, nb.repl = rp, int.pval = intpvalestimate, 
                            CIlow.ipv = intpvalCIlower, CIup.ipv = intpvalCIupper, 
                            int.power = intpowestimate, CIlow.ipo = intpowCIlower, 
                            CIup.ipo = intpowCIupper,
                            nsim.int = nsim.used.int)
    }
    
    class(sim.sum) = c("PAMM", "data.frame")
    sim.sum
  }


