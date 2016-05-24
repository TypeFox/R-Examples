gsubpop.sim <-
function(z.early=NULL,z1=z1,z2=z2,sprev=sprev,
        corr=NULL,selim=NULL,nsim=nsim,seed=12345678,level=level,
        select="thresh",wt=NULL,method="CT-SD") {


 # validate inputs
 if (length(corr)==1) {
   if (abs(corr)> 1) {
   stop("Correlation must be between -1 and 1")
  }
  }
  if (length(corr)> 1) {
       stop("Correlation must be vector length=1")
  }
  if (is.null(corr)) {
   corr <- 1
  }
 if (length(z1)!=2){
    stop("Z-statistic stage1 for sub-pop and full-pop: input vector length=2")
 }
 if (is.null(z.early)) {
   z.early <- z1
  }
 if (length(z.early)!=2){
    stop("Z-statistic early for sub-pop and full-pop: input vector length=2")
 }
 if (length(z2)!=4){
    stop("Z-statistic stage1 for sub-pop and full-pop: input vector length=4")
 }

 if (length(sprev)>2) {
    stop("sprev must be vector length=2")
 }
 if (length(sprev)==2) {
  if (sprev[1]<=0 | sprev[1]>=1){
    stop("sprev must be >0 and <1")
  }
  if (sprev[2]<=0 | sprev[2]>=1){
    stop("sprev must be >0 and <1")
  }
 }

 sel.options <- c("thresh","futility")
 isel <- as.integer(match(select, sel.options, -1))
 if (isel < 1) {
   stop("Unknown method: current option thresh and futility")
  } # end if
  if(select=="thresh"){
   if(selim[1]>selim[2]){
     stop("Limits for threshold rule: selim[1]<selim[2]")
   } # end if
  } # end if
  if(select=="futility"){
   if(is.null(selim)==TRUE){
    selim[1] <- z.early[1]
    selim[2] <- z.early[2]
   }
  } # end if

  nsim <- abs(round(nsim, 0))
  if (nsim > 1e+07) {
    stop("Maximum of 10,000,000 simulations allowed")
  }
  if (level >= 1 | level <= 0) {
    stop("Level must be between 0 and 1")
  }

  meth.options <- c("CT-Simes","CT-Bonferroni","CT-SD","CEF")
  imeth <- as.integer(match(method, meth.options, -1))
  if (imeth < 1) {
     stop("unknown method: current options CT-Simes, CT-Bonferroni, CT-SD or CEF")
  }


 # simes test function
 simes.test <- function(Z=Z,select=rep(1,2)) {
     hyp.comb <- list(NULL)
     hyp.comb[[1]] <- matrix(c(1,2),nrow=1,ncol=2)
     hyp.comb[[2]] <- matrix(c(1,2),nrow=2,ncol=1)
     rownames(hyp.comb[[1]]) <- 1
     rownames(hyp.comb[[2]]) <- 1:2
     colnames(hyp.comb[[1]]) <- c("Hs","Hf")
     colnames(hyp.comb[[2]]) <- c("Hsf")
     psimes.test <- list(NULL)
     zscores <- list(NULL)
     psimes.test[[1]] <- matrix(pnorm(Z),nrow=1,ncol=2)
     if (sum(select)==2){
      psimes.test[[2]] <- matrix(min(2*min(pnorm(Z)),max(pnorm(Z))),nrow=1,ncol=1)
     } else {
      psimes.test[[2]] <- matrix(min(pnorm(Z)),nrow=1,ncol=1)
     } # end if
     for (j in 1:2){
      colnames(psimes.test[[j]]) <- colnames(hyp.comb[[j]])
      rownames(psimes.test[[j]]) <- 1
      zscores[[j]] <- qnorm(psimes.test[[j]])
     } # end j
      list(pvalues = psimes.test, zscores = zscores, hyp.comb = hyp.comb)
 } # end simes.test



 # bonferroni test function
 bonferroni.test <- function(Z=Z,select=rep(1,2)) {
     hyp.comb <- list(NULL)
     hyp.comb[[1]] <- matrix(c(1,2),nrow=1,ncol=2)
     hyp.comb[[2]] <- matrix(c(1,2),nrow=2,ncol=1)
     rownames(hyp.comb[[1]]) <- 1
     rownames(hyp.comb[[2]]) <- 1:2
     colnames(hyp.comb[[1]]) <- c("Hs","Hf")
     colnames(hyp.comb[[2]]) <- c("Hsf")
     pbonferroni.test <- list(NULL)
     zscores <- list(NULL)
     pbonferroni.test[[1]] <- matrix(pnorm(Z),nrow=1,ncol=2)
     if (sum(select)==2){
      pbonferroni.test[[2]] <- matrix(min(c(1,2*min(pnorm(Z)))),nrow=1,ncol=1)
     } else {
      pbonferroni.test[[2]] <- matrix(min(pnorm(Z)),nrow=1,ncol=1)
     } # end if
     for (j in 1:2){
      colnames(pbonferroni.test[[j]]) <- colnames(hyp.comb[[j]])
      rownames(pbonferroni.test[[j]]) <- 1
      zscores[[j]] <- qnorm(pbonferroni.test[[j]])
     } # end j
      list(pvalues = pbonferroni.test, zscores = zscores, hyp.comb = hyp.comb)
 } # end bonferroni.test


 # binormal test function
 binormal.test <- function (Z = Z, select = rep(1, 2), corr = corr) 
 {
     hyp.comb <- list(NULL)
     hyp.comb[[1]] <- matrix(c(1,2),nrow=1,ncol=2)
     hyp.comb[[2]] <- matrix(c(1,2),nrow=2,ncol=1)
     rownames(hyp.comb[[1]]) <- 1
     rownames(hyp.comb[[2]]) <- 1:2
     colnames(hyp.comb[[1]]) <- c("Hs","Hf")
     colnames(hyp.comb[[2]]) <- c("Hsf")
     pbinormal.test <- list(NULL)
     zscores <- list(NULL)
     int_binormal <- function(z,corr,zm) {
         pnorm((corr*z-zm)/sqrt(1-corr*corr))*dnorm(z)
     } 
     pbinormal.test[[1]] <- matrix(pnorm(Z),nrow=1,ncol=2)
     Zmin <- min(Z)
     if (sum(select)==2){
        binormal_integral <- integrate(int_binormal, lower = Zmin, 
                  upper = Inf, zm = Zmin, corr = corr)
        pbinormal.test[[2]] <- matrix(abs(1 - binormal_integral$value),nrow=1,ncol=1)
     } else {
        pbinormal.test[[2]] <- matrix(min(pnorm(Z)),nrow=1,ncol=1)
     }
     for (j in 1:2){
       colnames(pbinormal.test[[j]]) <- colnames(hyp.comb[[j]])
       rownames(pbinormal.test[[j]]) <- 1
       zscores[[j]] <- qnorm(pbinormal.test[[j]])
      }
    list(pvalues = pbinormal.test, zscores = zscores, hyp.comb = hyp.comb)
 }


  # set seed
  set.seed(abs(round(seed, 0)))

  # weights
  if (is.null(wt)){
        stop("wt must be between 0 and 1")
  } else {
  if (length(wt) > 1) {
        stop("wt must be vector length=1")
  }
  if (length(wt) == 1) {
        if (abs(wt) > 1) {
         stop("wt must be between 0 and 1")
        } else {
      weight <- wt
     }
   }
  }

  # stage 1
   z.mean1 <- c(z1,z.early)
   varcov.bdiag <- matrix(c(1,sqrt(sprev[1]),sqrt(sprev[1]),1),nrow=2,ncol=2)
   varcov.boffdiag <- matrix(c(corr,corr*sqrt(sprev[1]),corr*sqrt(sprev[1]),corr),nrow=2,ncol=2)
   varcov.mat1 <- matrix(NA,nrow=4,ncol=4)
   varcov.mat1[1:2,1:2] <- varcov.bdiag; varcov.mat1[3:4,3:4] <- varcov.bdiag
   varcov.mat1[1:2,3:4] <- varcov.boffdiag; varcov.mat1[3:4,1:2] <- varcov.boffdiag

  for (k in 1:nsim){
   if (k==1){
    isubpop <- 0; ifullpop <- 0; ibothpop <- 0; istop <- 0
    sub.cnt <- 0; full.cnt <- 0; both.cnt <- 0; stop.cnt <- 0
    reject_both <- 0; reject.int.s <- reject.int.f <- reject.int.b <- 0
   }
   if(corr==1){
     z1.randi <- rmvnorm(n = 1, mean = z1, sigma = varcov.bdiag)
     otest.stat1 <- c(z1.randi, z1.randi - (z1 - z.early))
   } else if(corr==-1){
     z1.randi <- rmvnorm(n = 1, mean = z1, sigma = varcov.bdiag)
     otest.stat1 <- c(z1.randi, z1 - (z1.randi - z.early))
   } else {
     otest.stat1 <- rmvnorm(n = 1, mean = z.mean1, sigma = varcov.mat1)
   }
   test.stat1 <- otest.stat1[1:2]
   etest.stat1 <- otest.stat1[3:4]
 
      # CEF
   c2 <- qmvnorm((1-level),corr=varcov.bdiag,tail="upper")$quantile
   cond1 <- (c2-sqrt(weight)*test.stat1[1])/sqrt(1-weight)
   cond2 <- (c2-sqrt(weight)*test.stat1[2])/sqrt(1-weight)
   qstat <- 1-pmvnorm(lower=c(cond1,cond2),upper=c(Inf,Inf),mean=c(0,0),
                 corr=varcov.bdiag)[1]

   # stage 1 p-values
   if (method=="CT-Simes"){
    pstage1 <- simes.test(test.stat1)
   } else if (method=="CT-SD") {
    pstage1 <- binormal.test(test.stat1,corr=sqrt(sprev[1])) 
   } else if (method=="CT-Bonferroni") {
    pstage1 <- bonferroni.test(test.stat1)
   } else if (method=="CEF"){
    pstage1 <- NULL
   }

  # threshold rule
  if(select=="thresh"){
   # difference between test statistics
   dtest.stat1 <- etest.stat1[1]-etest.stat1[2]
   if (dtest.stat1<=selim[1]){
    select.pop <- "sub"
   } else if (dtest.stat1>selim[2]){
    select.pop <- "full"
   } else{
    select.pop <- "both"
   } # end if
  } # end if thresh

 
  # futility rule
  if(select=="futility"){
   if (etest.stat1[1]<selim[1] & etest.stat1[2]<selim[2]){
    select.pop <- "both"
   } else if (etest.stat1[1]<selim[1] & etest.stat1[2]>=selim[2]){
    select.pop <- "sub"
   } else if (etest.stat1[1]>=selim[1] & etest.stat1[2]<selim[2]){
    select.pop <- "full"
   } else {
    select.pop <- "stop"
   } # end if
  } # end if futility

 # select procedure based on decision rule

   if (select.pop=="sub"){

    sub.cnt <- sub.cnt+1
    z.mean2 <- z2[1]
    test.stat2 <- c(rnorm(1,mean=z.mean2,sd=1),Inf)
    if (method=="CT-Simes"){
     pstage2 <- simes.test(test.stat2,select=c(1,0))
    } else if (method=="CT-SD") {
     pstage2 <- binormal.test(test.stat2,corr=sqrt(sprev[2]),select=c(1,0))
    } else if (method=="CT-Bonferroni") {
     pstage2 <- bonferroni.test(test.stat2,select=c(1,0))
    }

    if (method=="CEF"){
     fogs <- matrix(0,nrow=1,ncol=2)  
     rownames(fogs) <- 1; colnames(fogs) <- c("Hs","Hf")
     if (pnorm(test.stat2[1])<qstat){
      fogs[1,1] <- 1
      reject.int.s <- reject.int.s+1
     }
     reject_test <- list(reject=fogs)
    } else {
     comb.test <- combn.test(pstage1,pstage2,weight=weight,method="invnorm")
     reject_test <- hyp.test(comb.test,level=level,full.hyp=FALSE)
   
     # save p-values
     reject.int.s <- reject.int.s+as.numeric((1-pnorm(as.numeric(comb.test$zscores[[2]])))<level)

    }
    if (isubpop==0){
     subpop <- reject_test$reject
     isubpop <- 1
    } else {
     subpop <- subpop + reject_test$reject
    }

   } else if (select.pop=="full") {

    full.cnt <- full.cnt+1
    z.mean2 <- z2[2]
    test.stat2 <- c(Inf,rnorm(1,mean=z.mean2,sd=1))

    if (method=="CT-Simes"){
     pstage2 <- simes.test(test.stat2,select=c(0,1))
    } else if (method=="CT-SD") {
     pstage2 <- binormal.test(test.stat2,corr=sqrt(sprev[2]),select=c(0,1))
    } else if (method=="CT-Bonferroni") {
     pstage2 <- bonferroni.test(test.stat2,select=c(0,1))
    }

    if (method=="CEF"){
     fogs <- matrix(0,nrow=1,ncol=2)  
     rownames(fogs) <- 1; colnames(fogs) <- c("Hs","Hf")
     if (pnorm(test.stat2[2])<qstat){
      fogs[1,2] <- 1
      reject.int.f <- reject.int.f+1
     }
     reject_test <- list(reject=fogs)
    } else {
      comb.test <- combn.test(pstage1,pstage2,weight=weight,method="invnorm")
      reject_test <- hyp.test(comb.test,level=level,full.hyp=FALSE)

     # save p-values
     reject.int.f <- reject.int.f+as.numeric((1-pnorm(as.numeric(comb.test$zscores[[2]])))<level)

    }
    if (ifullpop==0){
     fullpop <- reject_test$reject
     ifullpop <- 1
    } else {
     fullpop <- fullpop + reject_test$reject
    }

   } else if (select.pop=="both"){

    both.cnt <- both.cnt+1
    z.mean2 <- z2[3:4]
    varcov.mat2 <- matrix(c(1,sqrt(sprev[2]),sqrt(sprev[2]),1),nrow=2,ncol=2)
    test.stat2 <- rmvnorm(n=1,mean=z.mean2,sigma=varcov.mat2)

    if (method=="CT-Simes"){
     pstage2 <- simes.test(test.stat2)
    } else if (method=="CT-SD") {
     pstage2 <- binormal.test(test.stat2,corr=sqrt(sprev[2]))
    } else if (method=="CT-Bonferroni") {
     pstage2 <- bonferroni.test(test.stat2)
    }
 
    if (method=="CEF"){
     zmin <- min(test.stat2)
     int_funct <- function(z,tau,zm) {
         pnorm((sqrt(tau)*z-zm)/sqrt(1-tau))*dnorm(z)
     }
     p_int <- 1-integrate(int_funct,lower=zmin,upper=Inf,zm=zmin,tau=sprev[2])$value
     fogs <- matrix(0,nrow=1,ncol=2)  
     rownames(fogs) <- 1; colnames(fogs) <- c("Hs","Hf")

     if(p_int<qstat){
      # reject p_int if less than qstat
      reject.int.b <- reject.int.b+1

     if(zmin==test.stat2[2]){
     # reject full
     fogs[1,2] <- 1
     } else if(test.stat2[2]<cond2){
      fogs[1,2] <- 1
     } # end if

     if(zmin==test.stat2[1]){
     # reject sub
      fogs[1,1] <- 1
     } else if (test.stat2[1]<cond1){
      fogs[1,1] <- 1
     } # end if

     }
     reject_test <- list(reject=fogs)
    } else {
     comb.test <- combn.test(pstage1,pstage2,weight=weight,method="invnorm")
     reject_test <- hyp.test(comb.test,level=level,full.hyp=FALSE)

     # save p-values
     reject.int.b <- reject.int.b+as.numeric((1-pnorm(as.numeric(comb.test$zscores[[2]])))<level)


    }
    ireject_both <- reject_test$reject[1]+reject_test$reject[2]
    if (ireject_both==2){
     reject_both <- reject_both+1
    }
    if (ibothpop==0){
     bothpop <- reject_test$reject
     ibothpop <- 1
    } else {
     bothpop <- bothpop + reject_test$reject
    }
    } else {
  stop.cnt <- stop.cnt+1
  if (method=="CEF"){
   pstage1 <- simes.test(test.stat1)
   comb.test <- combn.test(pstage1,pstage1,weight=weight,method="invnorm")
   reject_test <- hyp.test(comb.test,level=level,full.hyp=FALSE)
   fstop <- reject_test
   fstop <- c(0,0)
   istop <- 1
  } else {
   comb.test <- combn.test(pstage1,pstage1,weight=weight,method="invnorm")
   reject_test <- hyp.test(comb.test,level=level,full.hyp=FALSE)
   fstop <- reject_test
   fstop <- c(0,0)
   istop <- 1
  }
 }
 }

 # results
 if (isubpop==0){
  subpop <- reject_test$reject
  subpop[1,] <- c(0,0)
 } # end if
 if (ifullpop==0){
  fullpop <- reject_test$reject
  fullpop[1,] <- c(0,0)
 } # end if
 if (ibothpop==0){
  bothpop <- reject_test$reject
  bothpop[1,] <- c(0,0)
 } # end if
 sim.res <- matrix(c(subpop,0,reject.int.s,sub.cnt,
                     fullpop,0,reject.int.f,full.cnt,
                     bothpop,reject_both,reject.int.b,both.cnt),
                     nrow=3,ncol=5,byrow=T)
 rownames(sim.res) <- c("sub","full","both")
 colnames(sim.res) <- c(colnames(subpop),"Hs+Hf","Hs+f","n")

 # output
 list(results=sim.res)

}
