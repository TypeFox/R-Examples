cox.aalenBase<-function (times, fdata, designX, designG, status,
    id, clusters, Nit = 5, beta = 0, weights=NULL, detail = 0, 
    sim = 1, antsim = 1000, weighted.test= 0, robust = 1, 
    ratesim = 0, residuals = 0, covariance = 1,
    resample.iid=0,namesZ=NULL,namesX=NULL,beta.fixed=0,strata=NULL,
    entry=NULL,offsets=0,exactderiv=1,max.timepoint.sim=100,silent=1,basesim=0)  ## {{{
{
  additive.resamp <-0; ridge <- 0; XligZ <- 0; 
  Ntimes <- length(times)
  sim <- rep(sim,3); 
  sim[2] <- basesim; ##  tmp tmp tmp  
###  if (ratesim==0 && (max(by(status,clusters,sum))>1)) mjump <- 1 else mjump <- 0; 
  mjump <- 0        ### not working
  sim[3] <- mjump;  ### multiple jumps within cluster

  if (ratesim==0 && mjump==1) cat("Multiple jumps in some clusters, use rate.sim=1\n"); 

  if (is.null(max.timepoint.sim)) max.timepoint.sim <- Ntimes; 
  if (max.timepoint.sim< Ntimes)  {   ## {{{
       qq <- quantile(times, probs = seq(0, 1, length=max.timepoint.sim))  
       qqc <- cut(times, breaks = qq, include.lowest = TRUE)    
       time.group <- as.integer(factor(qqc, labels = 1:(max.timepoint.sim-1)))
       time.resolution <- qq
       mts <- max.timepoint.sim; 
  } else { time.resolution <- qq <- times; time.group <- 0:(Ntimes-1); 
            max.timepoint.sim <- mts <- Ntimes;
  } ## }}}
  if (robust==1) antclust <- fdata$antclust else antclust <- 1

  designX <- as.matrix(designX); designG <- as.matrix(designG)
  if (is.matrix(designX) == TRUE) px <- as.integer(dim(designX)[2])
  if (is.matrix(designX) == TRUE) nx <- as.integer(dim(designX)[1])
  if (is.matrix(designG) == TRUE) pg <- as.integer(dim(designG)[2])
  if (is.matrix(designG) == TRUE) ng <- as.integer(dim(designG)[1])
  if (nx != ng) print(" X design and Z not same number of rows\n")

  if (is.null(weights)==FALSE) mw<-1 else { mw <- 0; weights <- rep(1, nx);}
  if (sum(offsets)==0) mof <- 0 else mof <- 1; nb <- 1; aalen <- 1
  if (covariance == 1) covs <- matrix(0, mts, px * px) else covs <- 0
  cumAi <- 0; dNit <- 0 ; dM.iid<-0; 
  gammaiid <- matrix(0, pg, fdata$antclust * 1)
  if (residuals == 1)  cumAi <- matrix(0, Ntimes , fdata$antpers * 1) 
###  if (residuals == 1)  dNit <- matrix(0, Ntimes , fdata$antpers * 1) 
  if (residuals == 2)  cumAi <- rep(0, fdata$antpers * 1) 
  cumint <- vcum <- matrix(0, Ntimes , px + 1); 
  Rvcu <- matrix(0, mts, px + 1); 
  if (sum(abs(beta)) == 0) betaS <- rep(0, pg) else betaS <- beta
  if (length(betaS)!=pg) betaS <- rep(betaS[1],pg); 
  score <- betaS; 
  loglike <- rep(0,2); 
  RVarbeta <- Iinv <- Varbeta <- matrix(0, pg, pg); 
  if (sim[1] == 1) Uit <- matrix(0, mts , 50 * pg) else Uit <- NULL
  if (additive.resamp == 1) baseproc <- matrix(0, mts, 50 * px) else baseproc <- NULL
  if (resample.iid == 1) biid <- matrix(0, mts, fdata$antclust * px) else biid <- NULL; 
  gamiid<- matrix(0,fdata$antclust,pg); 
  test <- matrix(0, antsim, 2 * px); testOBS <- rep(0, 2 * px)
  unifCI <- c(); testval <- c(); rani <- -round(runif(1) * 10000)
  Ut <- matrix(0, mts , pg + 1); simUt <- matrix(0, antsim, pg)
  var.score<- matrix(0, Ntimes , pg + 1); 

  if (is.diag(  t(designX) %*% designX  )==TRUE) stratum <- 1 else stratum <- 0
  if (!is.null(strata)) ostratum <- c(stratum[1],length(unique(strata)),strata) else if (stratum==1) ostratum <- c(stratum,px,designX %*% (0:(px-1))) else ostratum <- c(stratum,1,rep(0,nx)); 
  stratum <- ostratum

### print(stratum); ### print(weights); print(offsets)

  nparout <- .C("score", as.double(times), as.integer(Ntimes), 
                as.double(designX), as.integer(nx), as.integer(px), 
                as.double(designG), as.integer(ng), as.integer(pg), 
                as.integer(fdata$antpers), as.double(fdata$start), as.double(fdata$stop),
                as.double(betaS), as.integer(Nit), as.double(cumint), 
                as.double(vcum), as.double(weights), as.integer(mw), 
                as.double(loglike), as.double(Iinv), as.double(Varbeta), 
                as.integer(detail), as.double(offsets), as.integer(mof), 
                as.integer(sim), as.integer(antsim), as.integer(rani), 
                as.double(Rvcu), as.double(RVarbeta), as.double(test), 
                as.double(testOBS), as.double(Ut), as.double(simUt), 
                as.double(Uit), as.integer(XligZ), as.double(aalen), 
                as.integer(nb), as.integer(id), as.integer(status), 
                as.integer(weighted.test), as.double(dNit), as.integer(ratesim), 
                as.double(score), as.double(cumAi), as.double(gammaiid), 
                as.double(dM.iid), as.integer(residuals), as.integer(robust), 
                as.integer(covariance), as.double(covs), as.integer(additive.resamp),
                as.double(baseproc), as.integer(resample.iid), as.double(gamiid), 
                as.double(biid),as.integer(clusters),as.integer(antclust),
                as.double(var.score),as.integer(beta.fixed),
		as.double(weights),as.integer(entry) ,as.integer(exactderiv),
	        as.integer(time.group), as.integer(max.timepoint.sim),as.integer(stratum),
		as.integer(silent),PACKAGE = "timereg")

  Iinv <- matrix(nparout[[19]], pg, pg); 
  RVarbeta <- -matrix(nparout[[28]], pg, pg)
  rvcu <- matrix(nparout[[27]], mts , px + 1); ## convert to approx for times 
  Rvcu <- times; 
  for (i  in 2:(px+1)) Rvcu <- cbind(Rvcu,approx(rvcu[,1],rvcu[,i],times,f=0.5)$y)

  Varbeta <- -matrix(nparout[[20]], pg, pg); 
  cumint <- matrix(nparout[[14]], Ntimes, px + 1); 
  vcum <- matrix(nparout[[15]], Ntimes, px + 1)
  gamma <- matrix(nparout[[12]], pg, 1); 
  score <- matrix(nparout[[42]], pg, 1)

  Ut <- matrix(nparout[[31]], mts , pg + 1)
  if (beta.fixed==1) var.score<-matrix(nparout[[57]],Ntimes,pg+1)

  gammaiid <-t( matrix(nparout[[44]],pg,fdata$antclust * 1))
  gamiid<-matrix(nparout[[53]],fdata$antclust,pg);
  if (resample.iid==1)  {
    biid<-matrix(nparout[[54]],mts,fdata$antclust*px);
    B.iid<-list();
    for (i in (0:(fdata$antclust-1))*px) {
      B.iid[[i/px+1]]<-as.matrix(biid[,i+(1:px)]);
      colnames(B.iid[[i/px+1]])<-namesX; }
    colnames(gamiid)<-namesZ
  } else B.iid<-gamiid<-NULL; 

  if (covariance == 1) {
    covit <- matrix(nparout[[49]], mts, px * px)
    cov.list <- list()
    for (i in 1:mts) cov.list[[i]] <- matrix(covit[i,], px, px) 
  } else cov.list <- NULL
  if (residuals == 1) cumAi <- matrix(nparout[[43]],Ntimes,fdata$antpers * 1)
###  if (residuals == 1) dNit <- matrix(nparout[[40]],Ntimes,fdata$antpers * 1)
  if (residuals == 2) cumAi <- nparout[[43]]
  cumAi <- list(time = times, dM = cumAi)

    testUt <- test <- unifCI <- supUtOBS <- UIt <- testOBS <- testval <- pval.testBeq0 <- 
    pval.testBeqC <- obs.testBeq0 <- obs.testBeqC <- sim.testBeq0 <- sim.testBeqC <- testUt <- sim.supUt <- NULL 
	               
  if (sim[1] >= 1) {
    Uit <- matrix(nparout[[33]], mts, 50 * pg)
    UIt <- list()
    for (i in (0:49) * pg) UIt[[i/pg + 1]] <- as.matrix(Uit[, i + (1:pg)])
    simUt <- matrix(nparout[[32]], antsim, pg)
    supUtOBS <- apply(abs(as.matrix(Ut[, -1])), 2, max)
    testUt <- c()
    for (i in 1:pg) testUt <- c(testUt, pval(simUt[, i], supUtOBS[i]))
    sim.supUt <- as.matrix(simUt)

if (basesim==1) {
    test <- matrix(nparout[[29]], antsim, 2 * px)
	    testOBS <- nparout[[30]]
	    for (i in 1:(2 * px)) testval <- c(testval, pval(test[, i], testOBS[i]))
	    for (i in 1:px) unifCI <- c(unifCI, percen(test[, i], 0.95))
	    pval.testBeq0 <- as.vector(testval[1:px])
	    pval.testBeqC <- as.vector(testval[(px + 1):(2 * px)])
	    obs.testBeq0 <- as.vector(testOBS[1:px])
	    obs.testBeqC <- as.vector(testOBS[(px + 1):(2 * px)])
	    sim.testBeq0 <- as.matrix(test[, 1:px])
	    sim.testBeqC <- as.matrix(test[, (px + 1):(2 * px)])
    }
  }

  if (additive.resamp == 1) {
    baseproc <- matrix(nparout[[51]], mts, 50 * pg)
    additive.proc <- list()
    for (i in (0:49) * px)
      additive.proc[[i/px+1]]<-as.matrix(baseproc[,i+(1:px)])
  }
  else additive.proc <- NULL

  if (robust==0 & beta.fixed==0) var.score<-NULL;

  ud <- list(cum = cumint, var.cum = vcum, robvar.cum = Rvcu, 
             gamma = gamma, var.gamma = Varbeta, robvar.gamma = RVarbeta, 
             residuals = cumAi, loglike = nparout[[18]], D2linv = Iinv, 
             score = score, var.score=var.score, 
             pval.testBeq0 = pval.testBeq0, pval.testBeqC = pval.testBeqC, 
             obs.testBeq0 = obs.testBeq0, obs.testBeqC = obs.testBeqC, 
             sim.testBeq0 = sim.testBeq0, sim.testBeqC = sim.testBeqC, 
             conf.band = unifCI, test.procProp = Ut, sim.test.procProp = UIt,
             pval.Prop = testUt, sim.supProp = sim.supUt, covariance = cov.list, 
             B.iid=B.iid,gamma.iid=gammaiid,time.sim.resolution=qq,stratum=stratum)
  return(ud)
} ## }}}

