ipsHEGY <- function(data,itsd,Sel,pmax,CIPS=TRUE){

  csmean=apply(data,1,mean) 
if (CIPS==FALSE) csmean=0

STAT0=NULL
for (j in 1:ncol(data)){
hegy.out <- HEGY.test(wts=data[,j], itsd=itsd,regvar=csmean, selectlags=list(mode=Sel, Pmax=pmax))
STAT0=cbind(STAT0,hegy.out$stats[,1])
}
TEMP=data.frame(apply(STAT0,1,mean))
colnames(TEMP)=c("pHEGY")
colnames(STAT0)=colnames(data)
STAT=list(U_HEGY=STAT0,P_HEGY=TEMP)
return(STAT)
}


HEGY.test <- function(wts, itsd, regvar=0, selectlags=list(mode="signf", Pmax=NULL))
{
  s  <- frequency(wts); t0 <- start(wts); N <- length(wts)
  
  # Dependent variable.
  Deltay <- matrix(c(rep(NA, s), diff(wts, lag=s)), ncol=1)
  
  # Regressor variables.
  Intercept <- matrix(rep(1, N), ncol=1)[,itsd[1]]
  Trend     <- matrix(c(1:N), ncol=1)[,itsd[2]]
  SDummy <- data.frame(SeasDummy=SeasComponent(wts, "dummyCycle"))[,itsd[-c(1,2)]]
  
  if(!identical(regvar, 0) && length(names(regvar)) == 0)
    regvar <- data.frame(Regvar=regvar)
  
  if(identical(itsd, c(0,0,0)) && identical(regvar, 0))
    Mdetreg <- numeric(0)
  if(!identical(itsd, c(0,0,0)) && identical(regvar, 0))
    Mdetreg <- as.matrix(data.frame(Intercept, Trend, SDummy))
  if(!identical(itsd, c(0,0,0)) && !identical(regvar, 0))
    Mdetreg <- as.matrix(data.frame(Intercept, Trend, SDummy, regvar))
  if(identical(itsd, c(0,0,0)) && !identical(regvar, 0))
    Mdetreg <- as.matrix(data.frame(regvar))
  regvarnames <- dimnames(Mdetreg)[[2]]
  
  # HEGY regression without lags.
  Mhegyreg <- hegy.reg(wts)
  ifelse(length(Mdetreg) == 0,
         lmdf <- lmhegyp <- lm(Deltay[,1] ~ 0+Mhegyreg),
         lmdf <- lmhegyp <- lm(Deltay[,1] ~ 0+Mdetreg + Mhegyreg))
  
  # Lags selection.
  if(class(selectlags[[1]]) == "numeric"){
    selP <- selectlags[[1]]
  } else
    switch(selectlags[[1]],
           aic   = selP <- selPabic(lmdet=lmhegyp, type="aic", Pmax=selectlags[[2]]),
           bic   = selP <- selPabic(lmdet=lmhegyp, type="bic", Pmax=selectlags[[2]]),
           signf = selP <- selPsignf(lmdet=lmhegyp, cvref=NULL, Pmax=selectlags[[2]]),)
  
  # HEGY regression.
  # lmdetlag: regression with deterministic components and lags 
  # (without the hegy regressors).
  # lmhegy: lmdetlag including the hegy regressors.
  if(identical(selP, 0) || length(selP)==0){
    if(length(Mdetreg)==0){
      lmdetlag <- lm(Deltay[,1] ~ 0)
      lmhegy <- lmhegyout <- lm(Deltay[,1] ~ 0+Mhegyreg)
    } else{
      lmdetlag <- lm(Deltay[,1] ~ 0+Mdetreg)
      lmhegy <- lmhegyout <- lm(Deltay[,1] ~ 0+Mdetreg + Mhegyreg)
    }
  } else{
    Mlags <- ret(Deltay, max(selP)+2)[,-1]; aux <- dimnames(Mlags)[[2]]
    Mlags <- data.frame(Mlags[,selP]); lagnames <- aux[selP]
    Mlags <- as.matrix(Mlags)
    if(length(Mdetreg)==0){
      lmdetlag <- lm(Deltay[,1] ~ 0+ Mlags)
      lmhegy <- lmhegyout <- lm(Deltay[,1] ~ 0+Mhegyreg + Mlags)
    } else{
      lmdetlag <- lm(Deltay[,1] ~ 0+Mdetreg + Mlags)
      lmhegy <- lmhegyout <- lm(Deltay[,1] ~ 0+Mdetreg + Mhegyreg + Mlags)
    }
  }
  
  # lmhegy estimates.
  coefs <- coef(summary(lmhegy)); Ncoef <- length(coef(lmhegy))
  colnames <- dimnames(coefs)[[2]]
  ifelse(Ncoef==s, ref<-1, ref <- which(dimnames(coefs)[[1]] == "MhegyregYpi1"))
  
  if(ref > 1){
    regvarcoefs <- coefs[1:(ref-1),1:4]
    dim(regvarcoefs) <- c((ref-1), 4)
    dimnames(regvarcoefs) <- list(regvarnames, colnames)
  } else
    regvarcoefs <- NULL
  
  # hegyreg
  hegycoefs <- coefs[ref:(ref+s-1),1:4]
  dimnames(hegycoefs)[[1]] <- dimnames(Mhegyreg)[[2]]
  
  if((ref+s-1) < Ncoef){
    lagcoefs <- coefs[(ref+s):Ncoef,1:4]
    dim(lagcoefs) <- c(length((ref+s):Ncoef), 4); lagcoefs <- data.frame(lagcoefs)
    dimnames(lagcoefs) <- list(lagnames, colnames); lagcoefs <- as.matrix(lagcoefs)
  } else
    lagcoefs <- NULL
  
  # HEGY Statistics and p-values.
  # tpi
  if(s==4) c1 <- "HEGY";
  if(s==12) c1 <- "BM"
  c2 <- paste(itsd[1:2], sep="", collapse="")
  ifelse(itsd[3] != 0, c3 <-1, c3 <-0)
  
  for(i in 1:s){
    code <- paste(c(c1, c2, c3, "tpi", i), sep="", collapse="")
    hegycoefs[i,4] <- interpolpval(code=code, stat=hegycoefs[i,3], N=N)$pval
  }
  
  # tpi1-tpi2, Fpi_
  EtFst <- matrix(nrow=(s/2+3), ncol=2)
  EtFst[1:2,] <- hegycoefs[1:2,3:4]
  
  code <- paste(c(c1, c2, c3, "Foddeven"), sep="", collapse="")
  for(i in 1:(s/2-1)){
    lmpart <- update(lmdetlag, . ~ . + Mhegyreg[,c(1:s)[-c(2*i+1, 2*i+2)]])
    lmhegy <- update(lmpart, . ~ . + Mhegyreg[,c(2*i+1, 2*i+2)])
    EtFst[i+2,1] <- anova(lmpart, lmhegy, test="F")$F[2]
    EtFst[i+2,2] <- interpolpval(code=code, stat=EtFst[i+2,1], N=N)$pval
  }
  
  lmpart <- update(lmdetlag, . ~ . + Mhegyreg[,1])
  lmhegy <- update(lmpart, . ~ . + Mhegyreg[,2:s])
  EtFst[(s/2+2),1] <- anova(lmpart, lmhegy, test="F")$F[2]
  
  lmhegy <- update(lmdetlag, . ~ . + Mhegyreg[,1:s])
  EtFst[(s/2+3),1] <- anova(lmdetlag, lmhegy, test="F")$F[2]
  
  EtFsnames <- c("tpi_1", "tpi_2",
                 paste("Fpi", paste(seq(3,s,2),seq(4,s,2), sep=":"), sep="_"),
                 paste("Fpi_2:", s, sep=""), paste("Fpi_1:", s, sep=""))
  #dim(EtFst) <- c((s/2+3), 2);
  dimnames(EtFst) <- list(EtFsnames, c("Stat.", "p-value"))
  hegystat<-list(wts=wts, itsd=itsd, regvar=regvar, hegyreg=Mhegyreg, selectlags=selectlags,regvarcoefs=regvarcoefs, hegycoefs=hegycoefs, lagsorder=selP, lagcoefs=lagcoefs,res=residuals(lmhegy), lmhegy=lmhegyout, stats=EtFst)
  # new("hegystat", wts=wts, itsd=itsd, regvar=regvar, hegyreg=Mhegyreg, selectlags=selectlags,regvarcoefs=regvarcoefs, hegycoefs=hegycoefs, lagsorder=selP, lagcoefs=lagcoefs,res=residuals(lmhegy), lmhegy=lmhegyout, stats=EtFst)
  return(hegystat)
}



## Auxiliar functions

ret <- function(wts, k)
{
  N <- length(wts)
  wts.r <- matrix(NA, nrow=N, ncol=k); wts.r[,1] <- wts
  for(i in 1:(k-1))
    wts.r[(i+1):N, (i+1)] <- wts[1:(N-i)]
  wts.r <- data.frame(wts.r); dimnames(wts.r)[[2]] <- paste("Lag", 0:(k-1), sep=".")
  as.matrix(wts.r)
}

SeasComponent <- function(wts, type)
{
  s <- frequency(wts); t0 <- start(wts); N <- length(wts)
  
  if(type == "dummyCycle"){        
    # Barsky & Miron (1989)
    auxD <- matrix(0,nrow=N, ncol=s)
    sq   <- seq(1,N,s)
    k    <- 0
    
    for(j in 1:s){
      ifelse(sq[length(sq)] + k > N, n <- length(sq)-1, n <- N)
      for(i in 1:n)
        auxD[sq[i]+k,j] <- 1
      k <- k+1
    }
    VFE <- auxD
    if(t0[2] != 1){
      VFE <- matrix(nrow=N, ncol=s)
      VFE[,1:(t0[2]-1)] <- auxD[,(s-t0[2]+2):s]; VFE[,t0[2]:s] <- auxD[,1:(s-t0[2]+1)]
    }
    if(t0[2] == 1){ VFE <- auxD }
  }
  
  if(type == "trgCycle"){        
    # Granger & Newbold (1986)
    qq  <- s/2
    VFE <- matrix(nrow=N, ncol=(s-1))
    
    sq1 <- seq(1,qq*2,2)
    sq2 <- seq(2,qq*2,2)
    j   <- c(1:(qq-1))
    
    for(i in 1:N){
      for(k in 1:(s-qq-1)){
        VFE[i,sq1[k]] <- cos((j[k]*pi/qq)*i)
        VFE[i,sq2[k]] <- sin((j[k]*pi/qq)*i)
      }
      VFE[i,(s-1)] <- (-1)^i
    }
  }
  VFE
}

contts <- function(lm, a)
{
  var.coef <- vcov(lm)[a,a]
  se.coef  <- sqrt(var.coef)
  et       <- lm$coef[a]/se.coef
  list(se.coef=se.coef, t.stat=et)
}

# HEGY regressors

hegy.reg <- function(wts)
{
  s <- frequency(wts); ML <- ret(ret(wts, 2)[,2], s+1)

  if(s==4){
    y1 <- ML[,1] + ML[,2] + ML[,3] + ML[,4]
    y2 <- -(ML[,1] - ML[,2] + ML[,3] - ML[,4])
    y4 <- -(ML[,1] - ML[,3])
    y3 <- ret(y4, 2)[,2]
    Mypi <- data.frame(y1, y2, y3, y4)
  }

  if(s==12){
    y1 <- ML[,1] + ML[,2] + ML[,3] + ML[,4] + ML[,5] + ML[,6] +
          ML[,7] + ML[,8] + ML[,9] + ML[,10] + ML[,11] + ML[,12]
    y2 <- -(ML[,1] - ML[,2] + ML[,3] - ML[,4] + ML[,5] - ML[,6] +
            ML[,7] - ML[,8] + ML[,9] - ML[,10] + ML[,11] - ML[,12])
    y3 <- -(ML[,2] - ML[,4] + ML[,6] - ML[,8] + ML[,10] - ML[,12])
    y4 <- -(ML[,1] - ML[,3] + ML[,5] - ML[,7] + ML[,9] - ML[,11])
    y5 <- -0.5*(ML[,1] + ML[,2] - 2*ML[,3] + ML[,4] + ML[,5] - 2*ML[,6] +
                ML[,7] + ML[,8] - 2*ML[,9] + ML[,10] + ML[,11] - 2*ML[,12])
    y6 <- (sqrt(3)/2)*(ML[,1] - ML[,2] + ML[,4] - ML[,5] + ML[,7] - ML[,8] + ML[,10] - ML[,11])
    y7 <- 0.5*(ML[,1] - ML[,2] - 2*ML[,3] - ML[,4] + ML[,5] + 2*ML[,6] +
               ML[,7] - ML[,8] - 2*ML[,9] - ML[,10] + ML[,11] + 2*ML[,12])
    y8 <- -(sqrt(3)/2)*(ML[,1] + ML[,2] - ML[,4] - ML[,5] + ML[,7] + ML[,8] - ML[,10] - ML[,11])
    y9 <- -0.5*(sqrt(3)*ML[,1] - ML[,2] + ML[,4] - sqrt(3)*ML[,5] + 2*ML[,6] - sqrt(3)*ML[,7] +
                ML[,8] - ML[,10] + sqrt(3)*ML[,11] - 2*ML[,12])
    y10 <- 0.5*(ML[,1] - sqrt(3)*ML[,2] + 2*ML[,3] - sqrt(3)*ML[,4] + ML[,5] - ML[,7] +
                sqrt(3)*ML[,8] - 2*ML[,9] + sqrt(3)*ML[,10] - ML[,11])
    y11 <- 0.5*(sqrt(3)*ML[,1] + ML[,2] - ML[,4] - sqrt(3)*ML[,5] - 2*ML[,6] - sqrt(3)*ML[,7] -
                ML[,8] + ML[,10] + sqrt(3)*ML[,11] + 2*ML[,12])
    y12 <- -0.5*(ML[,1] + sqrt(3)*ML[,2] + 2*ML[,3] + sqrt(3)*ML[,4] + ML[,5] - ML[,7] -
                 sqrt(3)*ML[,8] - 2*ML[,9] - sqrt(3)*ML[,10] - ML[,11])
    Mypi <- data.frame(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12)
  }

  dimnames(Mypi)[[2]] <- paste("Ypi", 1:s, sep="")
  Mypi <- as.matrix(Mypi)
  Mypi
}

## selectP

selPabic <- function(lmdet, type, Pmax=NULL)
{
  if(mode(Pmax) != "numeric") Pmax <- round(10*log10(length(lmdet$model[,1])))
  switch(type, aic = k <- 2,
               bic = k <- log(length(lmdet$model[,1])),)

  ML <- ret(lmdet$model[,1], Pmax+1)[,-1]
  drop <- NULL
  for(i in Pmax:1){
    ifelse(length(drop) == 0,
      lm1 <- lm(lmdet$model[,1] ~ 0+as.matrix(lmdet$model[,-1]) + ML[,c(1:Pmax)]),
      lm1 <- lm(lmdet$model[,1] ~ 0+as.matrix(lmdet$model[,-1]) + ML[,c(1:Pmax)[-drop]]))
    aic1 <- AIC(lm1, k=k)
    #ifelse(length(drop) == (Pmax-1),
    ifelse(identical(drop, as.numeric(c(2:Pmax))),
      aic2 <- AIC(lm(lmdet$model[,1] ~ 0+as.matrix(lmdet$model[,-1])), k=k),
      aic2 <- AIC(lm(lmdet$model[,1] ~ 0+as.matrix(lmdet$model[,-1]) + ML[,c(1:Pmax)[-c(drop,i)]]), k=k))

    if(aic1 >= aic2)
      drop <- c(drop, i)
  }

  ifelse(length(drop) == 0, lags <- c(1:Pmax), lags <- c(1:Pmax)[-drop])
  lags
}

selPsignf <- function(lmdet, cvref=1.65, Pmax=NULL)
{
  if(mode(cvref) != "numeric") cvref <- 1.65  # 10% approx.
  if(mode(Pmax) != "numeric") Pmax <- round(10*log10(length(lmdet$model[,1])))

  ref <- ncol(model.matrix(lmdet))
  ML <- ret(lmdet$model[,1], Pmax+1)[,-1]
  drop <- NULL; cond <- TRUE

  while(cond == TRUE){
    lmref <- lm(lmdet$model[,1] ~ 0+as.matrix(lmdet$model[,-1]) + ML)
    #Nreg <- length(coef(lmref))
    #tstats <- coef(summary(lmref))[(ref+1):Nreg,3]; drop <- which(abs(tstats) < cvref)  **
    rcoefs <- na.omit(coef(summary(lmref))[,3]); Nreg <- length(coef(lmref))
    aux1 <- which(is.na(coef(lmref)[(ref+1):Nreg]))
    tstats <- rcoefs[(ref+1):length(rcoefs)]; drop <- c(aux1, which(abs(tstats) < cvref))

    if(length(aux1) == 0){
      cond <- FALSE
      aux <- names(data.frame(ML))[-drop]
      lags <- as.numeric(substr(aux, 5, nchar(aux)))
    }
    if(length(aux1) == ncol(ML)){
      cond <- FALSE
      lags <- numeric(0)
    }
    if(length(aux1) > 0 && length(drop) < ncol(ML)){
      aux <- dimnames(ML)[[2]][-drop]
      ML <- as.matrix(data.frame(ML[,-drop]))
      dimnames(ML)[[2]] <- aux
    }
  }
  lags
}

interpolpval <- function(code, stat, N, swarn=TRUE)
{
  table <- lookupCVtable(code)
  alphas <- dimnames(table)[[2]]
  Nref <- as.numeric(dimnames(table)[[1]]) 
  pvlref <- 0

  selrow <- rep(NA, ncol(table))
  for (i in 1:ncol(table))
    selrow[i] <- approx(Nref, table[,i], N, rule = 2)$y

  pval <- approx(selrow, alphas, stat, rule = 2)$y
  if(is.na(approx(selrow, alphas, stat, rule = 1)$y))
  {
    if(substr(code, 1,2) == "DF" || substr(code, 6,8) == "tpi" || substr(code, 8,10) == "tpi")
    {
      if(pval == min(as.numeric(alphas))){
        if(swarn==TRUE)
          warning("p-value is smaller than printed p-value")
        pvalw <- paste("<", format.pval(pval), sep=" ")
      } else{
          if(swarn==TRUE)
            warning("p-value is greater than printed p-value")
          pvalw <- paste(">", format.pval(pval), sep=" ")
      }
    }
    if(substr(code, 1,2) == "KP" || substr(code, 1,2) == "CH" ||
       substr(code, 6,9) == "Fodd" || substr(code, 8,11) == "Fodd")
    {
      if(pval == min(as.numeric(alphas))){
        if(swarn==TRUE)
          warning("p-value is smaller than printed p-value")
        pvalw <- paste("<", format.pval(pval), sep=" ")
      } else{
          if(swarn==TRUE)
            warning("p-value is greater than printed p-value")
          pvalw <- paste(">", format.pval(pval), sep=" "); pvlref <- -1
      }
    }
  } else
    pvalw <- format.pval(pval)

  ref1 <- c(1, 0.1, 0.05, 0.025, 0.01, 0.001)
  ref2 <- c(" ", ".", "*", "**", "***")
  pvl <- ref2[length(which(pval <= ref1))+pvlref]

  list(pval=pval, pvalw=as.character(pvalw), pvlab=pvl)
}

lookupCVtable <- function(code)
{

  if(code=="BM000tpi1")
  {
    table <- data.frame(CV1 = c(-2.51, -2.52, -2.57),
                        CV2 = c(-2.18, -2.21, -2.24),
                        CV3 = c(-1.89, -1.91, -1.95),
                        CV4 = c(-1.58, -1.59, -1.62))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM100tpi1")
  {
    table <- data.frame(CV1 = c(-3.35, -3.40, -3.41),
                        CV2 = c(-3.06, -3.11, -3.12),
                        CV3 = c(-2.80, -2.85, -2.86),
                        CV4 = c(-2.51, -2.55, -2.57))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM101tpi1")
  {
    table <- data.frame(CV1 = c(-3.32, -3.37, -3.41),
                        CV2 = c(-3.02, -3.06, -3.12),
                        CV3 = c(-2.76, -2.81, -2.86),
                        CV4 = c(-2.47, -2.53, -2.57))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM110tpi1")
  {
    table <- data.frame(CV1 = c(-3.87, -3.92, -3.97),
                        CV2 = c(-3.58, -3.63, -3.67),
                        CV3 = c(-3.32, -3.37, -3.40),
                        CV4 = c(-3.06, -3.09, -3.12))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM111tpi1")
  {
    table <- data.frame(CV1 = c(-3.83, -3.85, -3.97),
                        CV2 = c(-3.54, -3.57, -3.67),
                        CV3 = c(-3.28, -3.32, -3.40),
                        CV4 = c(-2.99, -3.04, -3.12))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }

  if(code=="BM000tpi2")
  {
    table <- data.frame(CV1 = c(-2.53, -2.52, -2.57),
                        CV2 = c(-2.16, -2.20, -2.24),
                        CV3 = c(-1.87, -1.91, -1.95),
                        CV4 = c(-1.57, -1.59, -1.62))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM100tpi2")
  {
    table <- data.frame(CV1 = c(-2.48, -2.54, -2.57),
                        CV2 = c(-2.15, -2.20, -2.24),
                        CV3 = c(-1.89, -1.91, -1.95),
                        CV4 = c(-1.57, -1.59, -1.62))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM101tpi2")
  {
    table <- data.frame(CV1 = c(-3.28, -3.37, -3.41),
                        CV2 = c(-3.01, -3.07, -3.12),
                        CV3 = c(-2.76, -2.81, -2.86),
                        CV4 = c(-2.48, -2.52, -2.57))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM110tpi2")
  {
    table <- data.frame(CV1 = c(-2.52, -2.55, -2.57),
                        CV2 = c(-2.18, -2.20, -2.24),
                        CV3 = c(-1.88, -1.93, -1.95),
                        CV4 = c(-1.55, -1.60, -1.62))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM111tpi2")
  {
    table <- data.frame(CV1 = c(-3.31, -3.40, -3.41),
                        CV2 = c(-3.02, -3.08, -3.12),
                        CV3 = c(-2.75, -2.84, -2.86),
                        CV4 = c(-2.47, -2.54, -2.57))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }

  if(code=="BM000tpi3" || code=="BM000tpi5" || code=="BM000tpi7" || code=="BM000tpi9" || code=="BM000tpi11")
  {
    table <- data.frame(CV1 = -c(2.50, 2.52, 2.56),
                        CV2 = -c(2.16, 2.18, 2.23),
                        CV3 = -c(1.88, 1.90, 1.95),
                        CV4 = -c(1.55, 1.57, 1.59))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM100tpi3" || code=="BM100tpi5" || code=="BM100tpi7" || code=="BM100tpi9" || code=="BM100tpi11")
  {
    table <- data.frame(CV1 = -c(2.51, 2.56, 2.56),
                        CV2 = -c(2.16, 2.20, 2.23),
                        CV3 = -c(1.87, 1.90, 1.95),
                        CV4 = -c(1.54, 1.57, 1.59))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM101tpi3" || code=="BM101tpi5" || code=="BM101tpi7" || code=="BM101tpi9" || code=="BM101tpi11")
  {
    table <- data.frame(CV1 = -c(3.83, 3.86, 3.91),
                        CV2 = -c(3.51, 3.55, 3.61),
                        CV3 = -c(3.25, 3.29, 3.35),
                        CV4 = -c(2.95, 2.99, 3.05))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM110tpi3" || code=="BM110tpi5" || code=="BM110tpi7" || code=="BM110tpi9" || code=="BM110tpi11")
  {
    table <- data.frame(CV1 = -c(2.49, 2.53, 2.56),
                        CV2 = -c(2.16, 2.20, 2.23),
                        CV3 = -c(1.88, 1.91, 1.95),
                        CV4 = -c(1.54, 1.57, 1.59))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="BM111tpi3" || code=="BM111tpi5" || code=="BM111tpi7" || code=="BM111tpi9" || code=="BM111tpi11")
  {
    table <- data.frame(CV1 = -c(3.79, 3.85, 3.91),
                        CV2 = -c(3.50, 3.55, 3.61),
                        CV3 = -c(3.24, 3.29, 3.35),
                        CV4 = -c(2.95, 3.00, 3.05))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }

  if(code=="BM000tpi4" || code=="BM000tpi6" || code=="BM000tpi8" || code=="BM000tpi10" || code=="BM000tpi12")
  {
    table <- data.frame(CV1 = c(-2.31, -2.33, -2.30),
                        CV2 = c(-1.95, -1.96, -1.94),
                        CV3 = c(-1.63, -1.65, -1.63),
                        CV4 = c(-1.27, -1.28, -1.28),
                        CV5 = c(1.25, 1.27, 1.27),
                        CV6 = c(1.61, 1.63, 1.63),
                        CV7 = c(1.93, 1.94, 1.94),
                        CV8 = c(2.29, 2.32, 2.32))
    dimnames(table) <- list(c(240, 480, 10000),
                            c("0.01", "0.025", "0.05", "0.10", "0.90", "0.95", "0.975", "0.99"))
  }
  if(code=="BM100tpi4" || code=="BM100tpi6" || code=="BM100tpi8" || code=="BM100tpi10" || code=="BM100tpi12")
  {
    table <- data.frame(CV1 = c(-2.30, -2.32, -2.30),
                        CV2 = c(-1.93, -1.95, -1.94),
                        CV3 = c(-1.62, -1.63, -1.63),
                        CV4 = c(-1.27, -1.27, -1.28),
                        CV5 = c(1.24, 1.27, 1.27),
                        CV6 = c(1.60, 1.62, 1.63),
                        CV7 = c(1.91, 1.93, 1.94),
                        CV8 = c(2.28, 2.30, 2.32))
    dimnames(table) <- list(c(240, 480, 10000),
                            c("0.01", "0.025", "0.05", "0.10", "0.90", "0.95", "0.975", "0.99"))
  }
  if(code=="BM101tpi4" || code=="BM101tpi6" || code=="BM101tpi8" || code=="BM101tpi10" || code=="BM101tpi12")
  {
    table <- data.frame(CV1 = c(-2.61, -2.65, -2.72),
                        CV2 = c(-2.21, -2.25, -2.31),
                        CV3 = c(-1.85, -1.90, -1.95),
                        CV4 = c(-1.45, -1.49, -1.54),
                        CV5 = c(1.46, 1.49, 1.53),
                        CV6 = c(1.86, 1.91, 1.95),
                        CV7 = c(2.20, 2.25, 2.30),
                        CV8 = c(2.60, 2.63, 2.72))
    dimnames(table) <- list(c(240, 480, 10000),
                            c("0.01", "0.025", "0.05", "0.10", "0.90", "0.95", "0.975", "0.99"))
  }
  if(code=="BM110tpi4" || code=="BM110tpi6" || code=="BM110tpi8" || code=="BM110tpi10" || code=="BM110tpi12")
  {
    table <- data.frame(CV1 = c(-2.28, -2.30, -2.30),
                        CV2 = c(-1.93, -1.94, -1.94),
                        CV3 = c(-1.61, -1.63, -1.63),
                        CV4 = c(-1.25, -1.27, -1.28),
                        CV5 = c(1.24, 1.25, 1.27),
                        CV6 = c(1.59, 1.61, 1.63),
                        CV7 = c(1.90, 1.92, 1.94),
                        CV8 = c(2.26, 2.28, 2.32))
    dimnames(table) <- list(c(240, 480, 10000),
                            c("0.01", "0.025", "0.05", "0.10", "0.90", "0.95", "0.975", "0.99"))
  }
  if(code=="BM111tpi4" || code=="BM111tpi6" || code=="BM111tpi8" || code=="BM111tpi10" || code=="BM111tpi12")
  {
    table <- data.frame(CV1 = c(-2.57, -2.66, -2.72),
                        CV2 = c(-2.18, -2.27, -2.31),
                        CV3 = c(-1.85, -1.91, -1.95),
                        CV4 = c(-1.45, -1.49, -1.54),
                        CV5 = c(1.45, 1.49, 1.53),
                        CV6 = c(1.86, 1.90, 1.95),
                        CV7 = c(2.19, 2.25, 2.30),
                        CV8 = c(2.60, 2.64, 2.72))
    dimnames(table) <- list(c(240, 480, 10000),
                            c("0.01", "0.025", "0.05", "0.10", "0.90", "0.95", "0.975", "0.99"))
  }

  if(code=="BM000Foddeven")
  {
    table <- data.frame(CV1 = c(2.34, 2.38, 2.40),
                        CV2 = c(3.03, 3.08, 3.10),
                        CV3 = c(3.71, 3.78, 3.79),
                        CV4 = c(4.60, 4.70, 4.68))
    dimnames(table) <- list(c(240, 480, 10000), c("0.10", "0.05", "0.025", "0.01"))
  }
  if(code=="BM100Foddeven")
  {
    table <- data.frame(CV1 = c(2.32, 2.36, 2.40),
                        CV2 = c(3.01, 3.06, 3.10),
                        CV3 = c(3.68, 3.76, 3.79),
                        CV4 = c(4.60, 4.66, 4.68))
    dimnames(table) <- list(c(240, 480, 10000), c("0.10", "0.05", "0.025", "0.01"))
  }
  if(code=="BM101Foddeven")
  {
    table <- data.frame(CV1 = c(5.27, 5.42, 5.64),
                        CV2 = c(6.26, 6.42, 6.67),
                        CV3 = c(7.19, 7.38, 7.63),
                        CV4 = c(8.35, 8.60, 8.79))
    dimnames(table) <- list(c(240, 480, 10000), c("0.10", "0.05", "0.025", "0.01"))
  }
  if(code=="BM110Foddeven")
  {
    table <- data.frame(CV1 = c(2.30, 2.36, 2.40),
                        CV2 = c(2.97, 3.05, 3.10),
                        CV3 = c(3.64, 3.72, 3.79),
                        CV4 = c(4.53, 4.62, 4.68))
    dimnames(table) <- list(c(240, 480, 10000), c("0.10", "0.05", "0.025", "0.01"))
  }
  if(code=="BM111Foddeven")
  {
    table <- data.frame(CV1 = c(5.25, 5.44, 5.64),
                        CV2 = c(6.23, 6.43, 6.67),
                        CV3 = c(7.14, 7.35, 7.63),
                        CV4 = c(8.33, 8.52, 8.79))
    dimnames(table) <- list(c(240, 480, 10000), c("0.01", "0.025", "0.05", "0.10"))
  }

  if(code=="HEGY000tpi1")
  {
    table <- data.frame(CV1 = c(-2.72, -2.60, -2.62, -2.62),
                        CV2 = c(-2.29, -2.26, -2.25, -2.23),
                        CV3 = c(-1.95, -1.97, -1.93, -1.94),
                        CV4 = c(-1.59, -1.61, -1.59, -1.62))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY100tpi1")
  {
    table <- data.frame(CV1 = c(-3.66, -3.47, -3.51, -3.48),
                        CV2 = c(-3.25, -3.14, -3.17, -3.13),
                        CV3 = c(-2.96, -2.88, -2.89, -2.87),
                        CV4 = c(-2.62, -2.58, -2.58, -2.57))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY101tpi1")
  {
    table <- data.frame(CV1 = c(-3.77, -3.55, -3.56, -3.51),
                        CV2 = c(-3.39, -3.22, -3.23, -3.18),
                        CV3 = c(-3.08, -2.95, -2.94, -2.91),
                        CV4 = c(-2.72, -2.63, -2.62, -2.59))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY110tpi1")
  {
    table <- data.frame(CV1 = c(-4.23, -4.07, -4.09, -4.05),
                        CV2 = c(-3.85, -3.73, -3.75, -3.70),
                        CV3 = c(-3.56, -3.47, -3.46, -3.44),
                        CV4 = c(-3.21, -3.16, -3.16, -3.15))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY111tpi1")
  {
    table <- data.frame(CV1 = c(-4.46, -4.09, -4.15, -4.05),
                        CV2 = c(-4.04, -3.80, -3.80, -3.74),
                        CV3 = c(-3.71, -3.53, -3.52, -3.49),
                        CV4 = c(-3.37, -3.22, -3.21, -3.18))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }

  if(code=="HEGY000tpi2")
  {
    table <- data.frame(CV1 = c(-2.67, -2.61, -2.60, -2.60),
                        CV2 = c(-2.27, -222, -2.23, -2.24),
                        CV3 = c(-1.95, -1.92, -1.94, -1.95),
                        CV4 = c(-1.60, -1.57, -1.61, -1.61))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY100tpi2")
  {
    table <- data.frame(CV1 = c(-2.68, -2.61, -2.60, -2.58),
                        CV2 = c(-2.27, -2.24, -2.21, -2.22),
                        CV3 = c(-1.95, -1.95, -1.91, -1.92),
                        CV4 = c(-1.60, -1.60, -1.58, -1.59))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY101tpi2")
  {
    table <- data.frame(CV1 = c(-3.75, -3.60, -3.49, -3.50),
                        CV2 = c(-3.37, -3.22, -3.15, -3.16),
                        CV3 = c(-3.04, -2.94, -2.90, -2.89),
                        CV4 = c(-2.69, -2.63, -2.59, -2.60))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY110tpi2")
  {
    table <- data.frame(CV1 = c(-2.65, -2.58, -2.65, -2.59),
                        CV2 = c(-2.24, -2.24, -2.25, -2.25),
                        CV3 = c(-1.91, -1.94, -1.96, -1.95),
                        CV4 = c(-1.57, -1.60, -1.63, -1.62))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY111tpi2")
  {
    table <- data.frame(CV1 = c(-3.80, -3.60, -3.57, -3.52),
                        CV2 = c(-3.41, -3.22, -3.18, -3.18),
                        CV3 = c(-3.08, -2.94, -2.93, -2.91),
                        CV4 = c(-2.73, -2.63, -2.61, -2.60))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }

  if(code=="HEGY000tpi3")
  {
    table <- data.frame(CV1 = -c(2.66, 2.55, 2.58, 2.58),
                        CV2 = -c(2.23, 2.18, 2.21, 2.24),
                        CV3 = -c(1.93, 1.90, 1.92, 1.92),
                        CV4 = -c(1.52, 1.53, 1.56, 1.55))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY100tpi3")
  {
    table <- data.frame(CV1 = -c(2.64, 2.23, 1.90, 1.52),
                        CV2 = -c(2.61, 2.23, 1.90, 1.54),
                        CV3 = -c(2.53, 2.18, 1.88, 1.53),
                        CV4 = -c(2.57, 2.21, 1.90, 1.53))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY101tpi3")
  {
    table <- data.frame(CV1 = -c(4.31, 4.06, 4.06, 4.00),
                        CV2 = -c(3.92, 3.72, 3.72, 3.67),
                        CV3 = -c(3.61, 3.44, 3.44, 3.38),
                        CV4 = -c(3.24, 3.14, 3.11, 3.07))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY110tpi3")
  {
    table <- data.frame(CV1 = -c(2.68, 2.56, 2.56, 2.58),
                        CV2 = -c(2.27, 2.19, 2.20, 2.21),
                        CV3 = -c(1.92, 1.89, 1.90, 1.92),
                        CV4 = -c(1.52, 1.54, 1.52, 1.56))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }
  if(code=="HEGY111tpi3")
  {
    table <- data.frame(CV1 = -c(4.46, 4.12, 4.05, 4.04),
                        CV2 = -c(4.02, 3.76, 3.72, 3.69),
                        CV3 = -c(3.66, 3.48, 3.44, 3.41),
                        CV4 = -c(3.28, 3.14, 3.12, 3.10))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.01", "0.025", "0.05", "0.10"))
  }

  if(code=="HEGY000tpi4")
  {
    table <- data.frame(CV1 = c(-2.51, -2.43, -2.44, -2.43),
                        CV2 = c(-2.11, -2.01, -1.99, -1.98),
                        CV3 = c(-1.76, -1.68, -1.68, -1.65),
                        CV4 = c(-1.35, -1.32, -1.31, -1.30),
                        CV5 = c(1.33, 1.31, 1.30, 1.29),
                        CV6 = c(1.72, 1.67, 1.66, 1.67),
                        CV7 = c(2.05, 2.00, 1.99, 1.97),
                        CV8 = c(2.49, 2.40, 2.38, 2.36))
    dimnames(table) <- list(c(48, 100, 136, 200),
                            c("0.01", "0.025", "0.05", "0.10", "0.90", "0.95", "0.975", "0.99"))
  }
  if(code=="HEGY100tpi4")
  {
    table <- data.frame(CV1 = c(-2.44, -2.38, -2.36, -2.36),
                        CV2 = c(-2.06, -1.99, -1.98, -1.98),
                        CV3 = c(-1.72, -1.68, -1.68, -1.66),
                        CV4 = c(-1.33, -1.30, -1.31, -1.29),
                        CV5 = c(1.30, 1.28, 1.27, 1.28),
                        CV6 = c(1.68, 1.65, 1.65, 1.65),
                        CV7 = c(2.04, 1.97, 1.97, 1.96),
                        CV8 = c(2.41, 2.32, 2.31, 2.30))
    dimnames(table) <- list(c(48, 100, 136, 200),
                            c("0.01", "0.025", "0.05", "0.10", "0.90", "0.95", "0.975", "0.99"))
  }
  if(code=="HEGY101tpi4")
  {
    table <- data.frame(CV1 = c(-2.86, -2.78, -2.72, -2.74),
                        CV2 = c(-2.37, -2.32, -2.31, -2.33),
                        CV3 = c(-1.98, -1.96, -1.96, -1.96),
                        CV4 = c(-1.53, -1.53, -1.52, -1.54),
                        CV5 = c(1.54, 1.52, 1.51, 1.53),
                        CV6 = c(1.96, 1.93, 1.92, 1.95),
                        CV7 = c(2.35, 2.29, 2.28, 2.32),
                        CV8 = c(2.81, 2.73, 2.71, 2.78))
    dimnames(table) <- list(c(48, 100, 136, 200),
                            c("0.01", "0.025", "0.05", "0.10", "0.90", "0.95", "0.975", "0.99"))
  }
  if(code=="HEGY110tpi4")
  {
    table <- data.frame(CV1 = c(-2.41, -2.38, -2.36, -2.35),
                        CV2 = c(-2.05, -1.97, -1.97, -1.97),
                        CV3 = c(-1.70, -1.65, -1.64, -1.66),
                        CV4 = c(-1.33, -1.28, -1.29, -1.29),
                        CV5 = c(1.26, 1.28, 1.26, 1.26),
                        CV6 = c(1.64, 1.65, 1.62, 1.64),
                        CV7 = c(1.96, 1.98, 1.92, 1.96),
                        CV8 = c(2.37, 2.32, 2.31, 2.30))
    dimnames(table) <- list(c(48, 100, 136, 200),
                            c("0.01", "0.025", "0.05", "0.10", "0.90", "0.95", "0.975", "0.99"))
  }
  if(code=="HEGY111tpi4")
  {
    table <- data.frame(CV1 = c(-2.75, -2.76, -2.71, -2.65),
                        CV2 = c(-2.26, -2.32, -2.78, -2.27),
                        CV3 = c(-1.91, -1.94, -1.94, -1.92),
                        CV4 = c(-1.48, -1.51, -1.51, -1.48),
                        CV5 = c(1.51, 1.51, 1.53, 1.55),
                        CV6 = c(1.97, 1.92, 1.96, 1.97),
                        CV7 = c(2.34, 2.28, 2.31, 2.31),
                        CV8 = c(2.78, 2.69, 2.78, 2.71))
    dimnames(table) <- list(c(48, 100, 136, 200),
                            c("0.01", "0.025", "0.05", "0.10", "0.90", "0.95", "0.975", "0.99"))
  }

  if(code=="HEGY000Foddeven")
  {
    table <- data.frame(CV1 = c(2.45, 2.39, 2.41, 2.42),
                        CV2 = c(3.26, 3.12, 3.14, 3.16),
                        CV3 = c(4.04, 3.89, 3.86, 3.92),
                        CV4 = c(5.02, 4.89, 4.81, 4.81))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.10", "0.05", "0.025", "0.01"))
  }
  if(code=="HEGY100Foddeven")
  {
    table <- data.frame(CV1 = c(2.32, 2.35, 2.36, 2.37),
                        CV2 = c(3.04, 3.08, 3.00, 3.12),
                        CV3 = c(3.78, 3.81, 3.70, 3.86),
                        CV4 = c(4.78, 4.77, 4.73, 4.76))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.10", "0.05", "0.025", "0.01"))
  }
  if(code=="HEGY101Foddeven")
  {
    table <- data.frame(CV1 = c(5.50, 5.56, 5.56, 5.56),
                        CV2 = c(6.60, 6.57, 6.63, 6.61),
                        CV3 = c(7.68, 7.72, 7.66, 7.53),
                        CV4 = c(9.22, 8.74, 8.92, 8.93))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.10", "0.05", "0.025", "0.01"))
  }
  if(code=="HEGY110Foddeven")
  {
    table <- data.frame(CV1 = c(2.23, 2.31, 2.33, 2.34),
                        CV2 = c(2.95, 2.98, 3.04, 3.07),
                        CV3 = c(3.70, 3.71, 3.69, 3.76),
                        CV4 = c(4.64, 4.70, 4.57, 4.66))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.10", "0.05", "0.025", "0.01"))
  }
  if(code=="HEGY111Foddeven")
  {
    table <- data.frame(CV1 = c(5.37, 5.52, 5.55, 5.56),
                        CV2 = c(6.55, 6.60, 6.62, 6.57),
                        CV3 = c(7.70, 7.52, 7.59, 7.56),
                        CV4 = c(9.27, 8.79, 8.77, 8.96))
    dimnames(table) <- list(c(48, 100, 136, 200), c("0.10", "0.05", "0.025", "0.01"))
  }

  table
}
