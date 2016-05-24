`EAMM` <-
  function (numsim, group, repl, fixed = c(0, 1, 0), VI = seq(0.05, 
                                                              0.95, 0.05), VS = seq(0.05, 0.5, 0.05), CoIS = 0, relIS = "cor", n.X = NA, autocorr.X = 0,
            X.dist = "gaussian", intercept=0,heteroscedasticity=c("null"), mer.sim=TRUE, mer.model=NULL) 
  {
    o.warn <- getOption("warn")
    
    M <- NULL
    if (is.null(mer.model)) {
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
    }
    vgi <- numeric(length(VI) * length(VS))
    vgs <- numeric(length(VI) * length(VS))
    powersl <- numeric(numsim)
    pvalsl <- numeric(numsim)
    slpowestimate <- numeric(length(VI) * length(VS))
    slpowCIlower <- numeric(length(VI) * length(VS))
    slpowCIupper <- numeric(length(VI) * length(VS))
    slpvalestimate <- numeric(length(VI) * length(VS))
    slpvalCIlower <- numeric(length(VI) * length(VS))
    slpvalCIupper <- numeric(length(VI) * length(VS))
    powerint <- numeric(numsim)
    pvalint <- numeric(numsim)
    intpowestimate <- numeric(length(VI) * length(VS))
    intpowCIlower <- numeric(length(VI) * length(VS))
    intpowCIupper <- numeric(length(VI) * length(VS))
    intpvalestimate <- numeric(length(VI) * length(VS))
    intpvalCIlower <- numeric(length(VI) * length(VS))
    intpvalCIupper <- numeric(length(VI) * length(VS))
    nsim.used.sl <- numeric(length(VI) * length(VS))
    nsim.used.int <- numeric(length(VI) * length(VS))
    
    kk <- 0
    for (k in VI) {
      for (r in VS) {
        if (is.null(mer.model)){
          N <- group * repl
          n.x <- ifelse( is.na(n.X)==TRUE, N, n.X)
        }
        VR <- 1 - k
        if (VR >= 0) {
          for (i in 1:numsim) {
            options(warn=2)
            if (relIS == "cor") { CovIS <- CoIS * sqrt(k) * sqrt(r) }
            if (relIS == "cov") { CovIS <- CoIS }
            M <- matrix(c(k, CovIS, CovIS, r), ncol = 2)
            
            if (is.null(mer.model)){
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
                if (n.x>=repl) {
                  inief <- sample(1:(n.x-repl+1),group,replace=TRUE)
                  EFrk <- rep(inief,repl) + rep (0:(repl-1),each=group)
                  EF <- ef[EFrk]
                }
                if (n.x<repl) { 
                  EF <- numeric(N)
                  EF[1:(n.x*group)] <- rep(ef,each=group)
                  EF[(n.x*group+1):N] <- sample(ef,length((n.x*group+1):N),replace=TRUE)
                }
              }
              else { EF <- ef }
              X <- sort(rep(c(1:repl), group))
              db <- data.frame(ID = rep(1:group, repl), obs = 1:N, X = X, EF = EF)
              
              if (mer.sim == TRUE) {
                sigma <- sqrt(VR)
                beta <- c(intercept,fixed[3])
                names(beta) <- c("(Intercept)","EF")
                theta <- as.vector(chol(M)/sigma)[c(1,3,4)]
                names(theta) <- c("ID.(Intercept)","ID.EF.(Intercept)","ID.EF")
                params <- list(beta=beta, theta= theta, sigma=sigma)
                y <- simulate(formula(~EF + (EF | ID)),newdata=db,family=gaussian, newparams=params)
                db$Y <- y[,1]
              }
              else if (mer.sim == FALSE) {
                er <- numeric(length(N))
                if (Hetero=="null") (er <- rnorm(N, intercept, sqrt(VR)))
                if (Hetero=="power") (
                  for (n in 1:N) {er[n] <- rnorm(1, intercept, sqrt(VR*(het[1]+abs(EF[n])^het[2])^2))} )
                if (Hetero=="exp")  (
                  for (n in 1:N) {er[n] <- rnorm(1, intercept, sqrt(VR*exp(2*het[1]*EF[n])))} )
                db$error <- er
                x <- rmvnorm(group, c(0, 0), M, method = "svd")
                db$rand.int <- rep(x[, 1], repl)
                db$rand.sl <- rep(x[, 2], repl)
                db$Y <- db$rand.int + (db$rand.sl + FE) * db$EF + 
                  db$error		
              }
              else {}
              
              #models
              if (r > 0) {
                m.full <- try(lmer(Y ~ EF + (EF | ID), data = db), silent=TRUE)
                m.nocov <- try(lmer(Y ~ EF + (1 | ID) + (0 + EF | ID), data = db), silent=TRUE) 
                m.nosl <- try(lmer(Y ~ EF + (1 | ID), data = db), silent=TRUE)
                m.noint <- try(lmer(Y ~ EF + (0 + EF | ID), data = db), silent=TRUE) 
                
                #anosl <- anova(m.nocov, m.nosl, refit =FALSE)
                #powersl[i] <- anosl[2, "Pr(>Chisq)"] <= 0.05
                #pvalsl[i] <- anosl[2, "Pr(>Chisq)"]
                if (class(m.full)!="merModLmerTest" || class(m.nosl)!="merModLmerTest")  { 
                  powersl[i] <- NA
                  pvalsl[i] <- NA
                }
                else{
                  anoIxE <- anova(m.full, m.nosl, refit =FALSE)
                  powersl[i] <- anoIxE[2, "Pr(>Chisq)"] <= 0.05
                  pvalsl[i] <- anoIxE[2, "Pr(>Chisq)"]
                }
                
                if (class(m.nocov)!="merModLmerTest" || class(m.noint)!="merModLmerTest")  {                       
                  powerint[i] <- NA
                  pvalint[i] <- NA
                }
                else{
                  anoint <- anova(m.nocov, m.noint, refit =FALSE)
                  powerint[i] <- anoint[2, "Pr(>Chisq)"] <= 0.05
                  pvalint[i] <- anoint[2, "Pr(>Chisq)"]
                }                      
                
              }
              else {
                powersl[i] <- 0
                pvalsl[i] <- 1
                m1.lmer <- try(lmer(Y ~ EF + (1 | ID), data = db), silent =TRUE )
                if (class(m1.lmer)!="merModLmerTest") {
                  powerint[i] <- NA
                  pvalint[i] <- NA
                }
                else{
                  lrt1 <- rand(m1.lmer)
                  pvint <- lrt1[[1]][1,3]
                  powerint[i] <- pvint <= 0.05
                  pvalint[i] <- pvint
                }
              }
            }
            else if( !is.null(mer.model) ){
              if (length(mer.model)!=3)         stop("mer.model should be a list of a lmer model, an evironmental covariate and a random effect")
              group <- length(unique(mer.model[[1]]@frame[,mer.model[[3]]]))
              repl <- nrow(mer.model[[1]]@frame) / group
              sigma <- sqrt(VR)*getME(mer.model[[1]],"sigma")
              beta <- fixef(mer.model[[1]])
              #theta <- getME(mer.model[[1]],"theta") could be used as a based for model with multiple random effects
              theta <- as.vector(chol(M)/sqrt(VR))[c(1,3,4)]
              names(theta) <- c(paste(mer.model[[3]],"(Intercept)",sep="."), 
                                paste(mer.model[[3]],mer.model[[2]],"(Intercept)",sep="."), paste(mer.model[[3]],mer.model[[2]],sep="."))
              params <- list(beta=beta, theta= theta, sigma=sigma)
              form <- paste(" ~ ",paste(names(fixef(mer.model[[1]]))[-1],collapse=" + "), "+ (",mer.model[[2]],"|",mer.model[[3]],")")
              dat <- mer.model[[1]]@frame[,-1]
              y <- simulate(as.formula(form), newdata= dat, newparams=params, family=gaussian)
              dat$Y <- y[,1]
              formnc <- paste(" ~ ",paste(names(fixef(mer.model[[1]]))[-1],collapse=" + "), "+ ( 1 | ",mer.model[[3]],") + ( 0 +",mer.model[[2]],"|",mer.model[[3]],")")
              formns <- paste(" ~ ",paste(names(fixef(mer.model[[1]]))[-1],collapse=" + "), "+ ( 1 | ",mer.model[[3]],")")
              formni <- paste(" ~ ",paste(names(fixef(mer.model[[1]]))[-1],collapse=" + "), "+ ( 0 +",mer.model[[2]],"|",mer.model[[3]],")")
              
              if (r > 0) {
                m.full <- try(lmer(as.formula(paste("Y",form)), data = dat), silent=TRUE) 
                m.nocov <- try(lmer(as.formula(paste("Y",formnc)), data = dat), silent=TRUE) 
                m.nosl <- try(lmer(as.formula(paste("Y",formns)), data = dat), silent=TRUE) 
                m.noint <- try(lmer(as.formula(paste("Y",formni)), data = dat), silent=TRUE) 
                
                if (class(m.full)!="merModLmerTest" || class(m.nosl)!="merModLmerTest")  { 
                  powersl[i] <- NA
                  pvalsl[i] <- NA
                }
                else{
                  anoIxE <- anova(m.full, m.nosl, refit =FALSE)
                  powersl[i] <- anoIxE[2, "Pr(>Chisq)"] <= 0.05
                  pvalsl[i] <- anoIxE[2, "Pr(>Chisq)"] 
                }
                
                if (class(m.nocov)!="merModLmerTest" || class(m.noint)!="merModLmerTest")  {                       
                  powerint[i] <- NA
                  pvalint[i] <- NA
                }
                else{
                  anoint <- anova(m.nocov, m.noint, refit =FALSE)
                  powerint[i] <- anoint[2, "Pr(>Chisq)"] <= 0.05
                  pvalint[i] <- anoint[2, "Pr(>Chisq)"]
                }
              }
              else {
                powersl[i] <- 0
                pvalsl[i] <- 1
                m.nosl <- try(lmer(as.formula(paste("Y",formns)), data = dat), silent = TRUE)
                if (class(m.nosl)!="merModLmerTest" ) {
                  powerint[i] <- NA
                  pvalint[i] <- NA
                }
                else{
                  lrt1 <- rand(m.nosl)
                  pvint <- lrt1[[1]][1,3]
                  powerint[i] <- pvint <= 0.05
                  pvalint[i] <- pvint
                }
              }
            }
          }       	
          options(warn=o.warn)
          kk <- kk + 1
          vgi[kk] <- k
          vgs[kk] <- r
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
        else {
          options(warn=o.warn)
          kk <- kk + 1
          vgi[kk] <- k
          vgs[kk] <- r
          slpowestimate[kk] <- NA
          slpowCIlower[kk] <- NA
          slpowCIupper[kk] <- NA
          slpvalestimate[kk] <- NA
          slpvalCIlower[kk] <- NA
          slpvalCIupper[kk] <- NA
          intpowestimate[kk] <- NA
          intpowCIlower[kk] <- NA
          intpowCIupper[kk] <- NA
          intpvalestimate[kk] <- NA
          intpvalCIlower[kk] <- NA
          intpvalCIupper[kk] <- NA
          nsim.used.sl[kk] <- NA
          nsim.used.int[kk] <- NA
        }
      }
    }
    sim.sum <- data.frame(nb.ID = rep(group, length(VI) * length(VS)), 
                          nb.repl = rep(repl, length(VI) * length(VS)), VI = vgi, 
                          VS = vgs, int.pval = intpvalestimate, CIlow.ipv = intpvalCIlower, 
                          CIup.ipv = intpvalCIupper, int.power = intpowestimate, 
                          CIlow.ipo = intpowCIlower, CIup.ipo = intpowCIupper, 
                          sl.pval = slpvalestimate, CIlow.slpv = slpvalCIlower, 
                          CIup.slpv = slpvalCIupper, sl.power = slpowestimate, 
                          CIlow.slpo = slpowCIlower, CIup.slpo = slpowCIupper,
                          nsim.int = nsim.used.int, nsim.sl = nsim.used.sl)
    
    class(sim.sum) <- c("EAMM", "data.frame")
    sim.sum
  }
