
## Auxiliar functions.

ret <- function(vari, k)
{
  N <- length(vari)
  vret <- matrix(c(1:k*N),N,k)
  vret[,1] <- vari
  i <- 1
  while(i < k){
   vret[1:i,(i+1)] <- NA
   vret[(i+1):N, (i+1)] <- vari[1:(N-i)]
   i <- i+1
              }
  vret
}

ysooys <- function(yso, t0, N, s)
{
  index <- matrix(-9999, ncol=3, nrow=N)
  index[,3] <- c(1:N)
  index[,2] <- c(c(t0[2]:s), rep(1:s,N/s))[1:N]
  index[1:length(t0[2]:s),1] <- rep(t0[1], length(t0[2]:s))
  i <- 1
  while(index[N,1] == -9999)
  {
     iaux <- which(index[,1] == -9999)
     reps <- ifelse(length(iaux) >= s, reps <- s, reps <- length(iaux))
     index[iaux[1]:(iaux[1]+(reps-1)),1] <- rep(t0[1]+i, reps)
     i <- i+1
  }

# year and season to observation
  if(length(yso)==2)
  {
    quest1 <- which(c(index[,1] == yso[1]) == TRUE)
    quest2 <- which(c(index[,2] == yso[2]) == TRUE)
     i <- 1; out <- quest1[1]
    while(length(which(quest2 == quest1[i])) != 1){
      i <- i+1
      out <- quest1[i]
    }
  }

# observation to year and season
  if(length(yso)==1)
    out <- index[which(index[,3]==yso),1:2]

  list(out, index)
}

Dummy <- function(wts, y0, s0, yN, sN)
{
   aux   <- ysooys(1, start(wts), length(wts), frequency(wts))[[2]][,1]
   years <- seq(aux[1], aux[length(aux)])
   Mdum  <- rep(0, wts$N)
   obs0  <- (which(years==y0)-1) * frequency(wts) + s0
   obsN  <- (which(years==yN)-1) * frequency(wts) + sN
   Mdum[obs0:obsN] <- 1
   Mdum <- matrix(Mdum, ncol=1)
   Mdum
}

SeasDummy <- function(wts, s, t0, type)
{
  N <- length(wts)

  if(type == "alg"){        # Empleada en Barsky & Miron (1989)
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

  if(type == "trg"){        # Empleada en Granger & Newbold (1986)
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

## partsm functions.

fit.ar.par <- function(wts, type, detcomp, p)
{
  s <- frequency(wts)
  t0 <- start(wts)
  N <- length(wts)
  MLag <- ret(wts, (p+1))

  Interc <- rep(detcomp$regular[1], N)
  Trnd <- c(1:N)
  SDum <- SeasDummy(wts, s, t0, "alg")
  STrnd <- Trnd*SDum

  Mdetreg <- data.frame(Interc, Trnd, SeasDum=SDum)
  ref1 <- which(detcomp$regular != 0)
  aux1 <- Mdetreg[,ref1]

  Mdetreg <- data.frame(SeasDum=SDum, SeasTrnd=STrnd)
  ref2 <- c(detcomp$seasonal[1] * c(1:s), detcomp$seasonal[2] * c((s+1):(2*s)))
  aux2 <- Mdetreg[,ref2]

  if(length(detcomp$regvar) == 1)
    MDT <- as.matrix(data.frame(.=aux1, .=aux2))
  if(length(detcomp$regvar) > 1)
    MDT <- as.matrix(data.frame(.=aux1, .=aux2, RegVar=detcomp$regvar))

  switch(type,
    "AR" = {
            if(p == 0){
              lm.ar <- lm(MLag[,1] ~ 0+MDT)
              lm.par <- par.coeffs <- ar.coeffs <- NULL
              #out <- list(lm=lm1, sumlm=summary(lm1))
            }
            if(p > 0){
              lm.ar <- lm(MLag[,1] ~ 0+MLag[,2:(p+1)] + MDT)
              ar.coeffs <- matrix(lm.ar$coef[1:p], nrow=1, byrow=TRUE)
              lm.par <- par.coeffs <- NULL
              #out <- list(lm=lm1, sumlm=summary(lm1), ar.coeffs=lm.ar$coef[1:p])
            } },

    "PAR" = {
             if(p == 0){
               lm.par <- lm(MLag[,1] ~ 0+MDT); par.coeffs <- NULL
               lm.ar <- par.coeffs <- ar.coeffs <- NULL
               #out <- list(lm=lm1, sumlm=summary(lm1))
             }
             if(p > 0){
               Yperlag <- matrix(nrow=length(wts), ncol=(p*s))
               j <- 2; for(i in seq(1,(p*s),s)){
                 Yperlag[,i:(i+s-1)] <- MLag[,j]*SDum
                 j <- j+1
               }
               lm.par <- lm(MLag[,1] ~ 0+Yperlag + MDT)
               par.coeffs <- matrix(lm.par$coef[1:(p*s)], nrow=p, byrow=TRUE)
               lm.ar <- ar.coeffs <- NULL
               #out <- list(lm=lm1, sumlm=summary(lm1),
               #            ar.coeffs=matrix(lm.par$coef[1:(p*s)], nrow=p, byrow=TRUE))
             } }
  )
  ##~ out
  new("fit.partsm", type=type, p=p, lm.ar=lm.ar, lm.par=lm.par, ar.coeffs=ar.coeffs, par.coeffs=par.coeffs)
}

Fnextp.test <- function(wts, detcomp, p, type)
{
  if(type=="AR"){
    test.name <- "Test for the significance of the p+1 autoregressive parameters"
    lm1 <- fit.ar.par(wts, detcomp=detcomp, type=type, p=p)@lm.ar
    lm2 <- fit.ar.par(wts, detcomp=detcomp, type=type, p=p+1)@lm.ar
  }
  if(type=="PAR"){
    test.name <- "Test for the significance of the p+1 periodic autoregressive parameters"
    lm1 <- fit.ar.par(wts, detcomp=detcomp, type=type, p=p)@lm.par
    lm2 <- fit.ar.par(wts, detcomp=detcomp, type=type, p=p+1)@lm.par
  }

  RSS0 <- sum(lm1$residuals^2)
  RSS1 <- sum(lm2$residuals^2)

  df   <- c(frequency(wts), length(wts) - length(coef(lm2)) - p)
  Fnextp <- ((RSS0-RSS1)/df[1]) / (RSS1/(df[2]-df[1]))
  pval <- 1 - pf(q=Fnextp, df1=df[1], df2=df[2], log.p=FALSE)

  ref1 <- c(1, 0.1, 0.05, 0.01, 0.001)
  ref2 <- c(" ", ".", "*", "**", "***")
  pvl <- ref2[length(which(pval <= ref1) == TRUE)]

  new("Ftest.partsm", test.label="Fnextp", test.name=test.name, p=p,
      Fstat=Fnextp, df=df, pval=pval, pvl=pvl, h0md=lm1, hamd=lm2)
  ##~list(Fnextp=Fnextp, df=df, pvalue=pval, lm1=summary(lm1), lm2=summary(lm2))
}

Fpar.test <- function(wts, detcomp, p)
{
  urmd <- fit.ar.par(wts, type="PAR", detcomp, p)@lm.par
  rmd  <- fit.ar.par(wts, type="AR", detcomp, p)@lm.ar

  RSS0 <- sum(rmd$residuals^2)
  RSS1 <- sum(urmd$residuals^2)

  df   <- c((frequency(wts)-1)*p, length(wts) - length(coef(rmd)) - p)
  Fpar <- ((RSS0-RSS1)/df[1]) / (RSS1/(df[2]-df[1]))
  pval <- 1 - pf(q=Fpar, df1=df[1], df2=df[2], log.p=FALSE)

  ref1 <- c(1, 0.1, 0.05, 0.01, 0.001)
  ref2 <- c(" ", ".", "*", "**", "***")
  pvl <- ref2[length(which(pval <= ref1) == TRUE)]

  new("Ftest.partsm", test.label="Fpar", test.name="Test for periodicity in the autoregressive parameters",
      p=p, Fstat=Fpar, df=df, pval=pval, pvl=pvl, h0md=rmd, hamd=urmd)
  ##~list(Fpar=Fpar, df=df, pvalue=pval, urmd=summary(urmd), rmd=summary(rmd))
}

Fsh.test <- function(res, s)
{
  res2 <- res^2
  SDum <- SeasDummy(res2, s, c(1,1), "alg")

  lm1 <- lm(res2 ~ 1)
  lm2 <- update(lm1, . ~ . + SDum[,1:(s-1)])

  RSS0 <- sum(residuals(lm1)^2)
  RSS1 <- sum(residuals(lm2)^2)

  df <- c((s-1), length(res) - length(coef(lm1)))
  Fsh <- ((RSS0-RSS1)/df[1]) / (RSS1/(df[2]-df[1]))
  pval <- 1 - pf(q=Fsh, df1=df[1], df2=df[2], log.p=FALSE)

  ref1 <- c(1, 0.1, 0.05, 0.01, 0.001)
  ref2 <- c(" ", ".", "*", "**", "***")
  pvl <- ref2[length(which(pval <= ref1) == TRUE)]

  new("Ftest.partsm", test.label="Fsh", test.name="Test for seasonal heteroskedasticity",
      Fstat=Fsh, df=df, pval=pval, pvl=pvl, h0md=lm1, hamd=lm2)
}

##~ Representaci?n posible si las ra?ces del modelo en forma VQ son reales (ver).
fit.piar <- function(wts, detcomp, p, initvalues=NULL)
{
  if(p > 2)
    stop("This function is implemented only for PAR models up to order 2.\n")

# Create the matrix with all the possible regressors.

  s <- frequency(wts)
  N <- length(wts)
  ML <- ret(wts, (p+1))

  Interc <- rep(detcomp$regular[1], N)
  Trnd <- c(1:N)
  SDum <- SeasDummy(wts, s, start(wts), "alg")
  STrnd <- Trnd*SDum

  Mdetreg <- data.frame(Interc, Trnd, SeasDum=SDum)
  ref1 <- which(detcomp$regular != 0)
  aux1 <- Mdetreg[,ref1]

  Mdetreg <- data.frame(SeasDum=SDum, SeasTrnd=STrnd)
  ref2 <- c(detcomp$seasonal[1] * c(1:s), detcomp$seasonal[2] * c((s+1):(2*s)))
  aux2 <- Mdetreg[,ref2]

  if(length(detcomp$regvar) == 1)
    Mreg <- as.matrix(data.frame(.=aux1, .=aux2))
  if(length(detcomp$regvar) > 1)
    Mreg <- as.matrix(data.frame(.=aux1, .=aux2, RegVar=detcomp$regvar))

  Yperlag <- matrix(nrow=length(wts), ncol=(p*s))
  j <- 2
  for(i in seq(1,(p*s),s)){
    Yperlag[,i:(i+s-1)] <- ML[,j]*SDum
    j <- j+1
  }

  Mreg <- as.matrix(data.frame(Yperlag=Yperlag, Mreg=Mreg))

# Define the nls model for p=1 and p=2.

  aux.anlr <- paste("coef", 1:s, "*Mreg[,", 1:s, "]", sep="")
  nlr <- paste("(1/(", paste(paste("coef", 1:(s-1), sep=""), collapse="*"), "))", sep="")
  aux.anlr[s] <- paste(nlr, "*Mreg[,", s, "]", sep="")

  anlr <- paste(aux.anlr, collapse=" + ")

  if(p==2){
    aux.b <- paste("coef", s:(p*s-1), "*Mreg[,", 1:s, "]", sep="")
    b <- paste(aux.b, collapse=" + ")

    ref <- paste("coef", s:(p*s-1), "*coef", 0:(s-1), sep="")
    ref[1] <- paste("coef", s, "*", nlr, sep="")
    aux.ab <- paste(ref[1:s], "*Mreg[,", (s+1):(p*s), "]", sep="")
    ab <- paste(aux.ab, collapse=" - ")
  }

  if(ncol(Mreg) > (p*s))
    aux.det <- paste("coef", (s*p):(ncol(Mreg)-1), "*Mreg[,", (s*p+1):ncol(Mreg), "]", sep="")
  det <- paste(aux.det, collapse=" + ")

  if(p==1)
    form.rpar <- as.formula(paste("ML[,1] ~ ", paste(anlr, det, sep=" + ")))
  if(p==2)
    form.rpar <- as.formula(paste("ML[,1] ~ ", paste(anlr, b, sep=" + "), "-", paste(ab, det, sep=" + ")))

# Compute initial values for the nls model.

  if(class(initvalues) == "NULL")
  {
  if(p==1)
    init <- lm(ML[,1] ~ 0+Mreg)$coef[-s]

  #nlsinit <- list(init[1])
  #for(i in 2:length(init))
  #  nlsinit <- c(nlsinit, list(init[i]), recursive=FALSE)
  #names(nlsinit) <- paste("coef", 1:length(init), sep="")

  if(p==2){
    lmaux <- lm(ML[,1] ~ 0+Mreg[,-c((s+1):(s*2))])
    initalpha <- coef(lmaux)[1:s]

    filc <- rowSums(initalpha*SDum)
    regaux <- SDum*ML[,2] - filc*SDum*ML[,3]

    lmbeta <- lm((ML[2:N,1]-residuals(lmaux)) ~ 0+regaux[2:N,])
    initbeta <- coef(lmbeta)

    init <- c(initalpha[1:(s-1)], initbeta, coef(lmaux)[(s+1):length(coef(lmaux))])
  }
  }

  nlsinit <- list(init[1])
  for(i in 2:length(init))
    nlsinit <- c(nlsinit, list(init[i]), recursive=FALSE)
  names(nlsinit) <- paste("coef", 1:length(nlsinit), sep="")

# nls.

  nlsrpar <- summary(nls(form.rpar, start=nlsinit, trace=FALSE))

# Periodic differences of the original series.

  pdiff.coeffs <- c(coef(nlsrpar)[1:(s-1),1], (1/prod(coef(nlsrpar)[1:(s-1),1])))
  pdfilc <- rowSums(pdiff.coeffs*SDum)
  pdiff.data <- ts(ML[,1] - pdfilc * ML[,2], frequency=s, start=start(wts))

  if(p==1)
    par.coeffs <- matrix(pdiff.coeffs, nrow=p, byrow=TRUE)
  if(p>1 && p<=4)
    par.coeffs <- matrix(c(pdiff.coeffs, coef(nlsrpar)[((p-1)*s):(p*s-1),1]), nrow=p, byrow=TRUE)

  new("fit.piartsm", p=p, nls.parameters=nlsrpar$parameters, nls.res=nlsrpar$residuals, par.coeffs=par.coeffs,
      pdiff.data=pdiff.data)
  ## pdiff.coeffs=matrix(pdiff.coeffs, nrow=1)
  ##~list(nls.rpar=nlsrpar, parcoeffs=parcoeffs, pdiff.coeffs=matrix(pdiff.coeffs, nrow=1), pddata=pddata)
}

Fpari.piar.test <- function(wts, detcomp, p, type)
{
  ML  <- ret(wts, 3)
  SDum <- SeasDummy(wts, frequency(wts), start(wts), "alg")
  Y1  <- ML[,2]*SDum
  Y2  <- ML[,3]*SDum

  switch(type,
    "PARI1"  = MLh0d <- ts(ML[,1]-ML[,2], frequency=frequency(wts), start=start(wts)),
    "PARI-1" = MLh0d <- ts(ML[,1]+ML[,2], frequency=frequency(wts), start=start(wts))
  )

  if(p==1){
    h0md <- fit.ar.par(wts=MLh0d, detcomp=detcomp, type="AR", p=0)@lm.ar
    hamd <- fit.piar(wts=wts, detcomp=detcomp, p=1)
  }
  if(p==2){
    switch(type,
      "PARI1"  = MLr <- (ML[,2]-ML[,3])*SDum,
      "PARI-1" = MLr <- (ML[,2]+ML[,3])*SDum)

    ifelse(detcomp$regvar == 0,
      dch0 <- list(regular=detcomp$regular, seasonal=detcomp$seasonal, regvar=MLr),
      dch0 <- list(regular=detcomp$regular, seasonal=detcomp$seasonal,
                   regvar=as.matrix(data.frame(detcomp$regvar, MLr)))
    )
    h0md <- fit.ar.par(wts=MLh0d, detcomp=dch0, type="AR", p=0)@lm.ar
    hamd <- fit.piar(wts=wts, detcomp=detcomp, p=2)
  }

  RSS0 <- sum(residuals(h0md)^2)
  RSS1 <- sum(hamd@nls.res^2)

  df <- c((frequency(wts)-1), length(wts) - length(coef(h0md)) - p)
  Fstat <- ((RSS0-RSS1)/df[1]) / (RSS1/(df[2]-df[1]))
  pval <- 1 - pf(q=Fstat, df1=df[1], df2=df[2], log.p=FALSE)

  ref1 <- c(1, 0.1, 0.05, 0.01, 0.001)
  ref2 <- c(" ", ".", "*", "**", "***")
  pvl <- ref2[length(which(pval <= ref1) == TRUE)]

  new("Ftest.partsm", test.label=paste("F", type, sep="-"),
      test.name="Test for a parameter restriction in a PAR model", p=p,
      Fstat=Fstat, df=df, pval=pval, pvl=pvl, h0md=h0md, hamd=hamd@nls.parameters)
  ##~list(Fstat=Fstat, df=df, pvalue=pval, h0md=summary(h0md), hamd=summary(hamd))
}

LRurpar.test <- function(wts, detcomp, p)
{
  if(p > 2)
    stop("This function is implemented only for PAR models up to order 2.\n")

  lmpar <- fit.ar.par(wts, type="PAR", detcomp=detcomp, p=p)
  par.phis <- lmpar@par.coeffs
  RSS1 <- sum(residuals(lmpar@lm.par)^2)

  nlsrpar <- fit.piar(wts=wts, detcomp=detcomp, p=p)
  RSS0 <- sum(nlsrpar@nls.res^2)

  n <- length(residuals(lmpar@lm.par))
  LR <- n * log(RSS0/RSS1, base = exp(1))
  LRtau <- sign(prod(par.phis[1,]) - 1) * sqrt(LR)
  ##~ poner output de fit.ar.par como lista en la que se nombre a alpha y beta, en lugar de coger la primera fila de una matriz para hacer prod(par.phis[1,]).

  new("LRur.partsm", test.label="LRurpar",
      test.name="Likelihood ratio test for a single unit root in a PAR model", p=p,
      LR=LR, LRtau=LRtau, h0nls=nlsrpar@nls.parameters, halm=lmpar@lm.par)
  ##~ list(LR=LR, LRtau=LRtau)
}

##~ Hacer predictpiar usando m?todos para cada parte, PAR.MVrepr, Omegas,...
predictpiar <- function(wts, p, hpred)
{
  t0 <- start(wts)
  N <- length(wts)
  s <- frequency(wts)

  sdetcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
  nlspiar <- fit.piar(wts=wts, detcomp=sdetcomp, p=p)
  nlscoef <- nlspiar@nls.parameters[,1]

  pred.se <- matrix(nrow=as.integer(hpred/s), ncol=s)
  sigma <- sd(nlspiar@nls.res)

  #phi <- matrix(nlspiar@par.coeffs[1,], nrow=1)
  #MPhi <- PAR.MVrepr(phi=phi, s=s)
  MV.out <- PAR.MVrepr(nlspiar)
  Phi0 <- MV.out@Phi0
  Phi1 <- MV.out@Phi1

  if(p==1)
  {
  # Forecast standard error for h=1.

    pred.se[1,] <- sqrt( sigma^2 * diag(solve(Phi0) %*% t(solve(Phi0))) )

  # Forecast standard error for h=2,3,...,hpred/s.

    Gamma <- solve(Phi0) %*% Phi1
    for(h in 2:(hpred/s))
      pred.se[h,] <- sqrt( sigma^2 * diag(
                       (solve(Phi0) %*% t(solve(Phi0)) +
                        (h-1) * (Gamma %*% solve(Phi0)) %*% t(Gamma %*% solve(Phi0)))) )

  # Forecast.

    mus <- matrix(nlscoef[(p*s):((p+1)*s-1)], ncol=1)
    pred <- matrix(nrow=as.integer(hpred/s), ncol=s)
    Yret1 <- matrix(wts[(N-s+1):N], ncol=1)

    for(i in 1:as.integer(hpred/s)){
      #pred[i,] <- solve(Phi0) %*% mus + solve(Phi0) %*% Phi1 %*% Yret1
      pred[i,] <- solve(Phi0) %*% mus + Gamma %*% Yret1
      Yret1 <- matrix(pred[i,], ncol=1)
    }
  }

  if(p==2)
  {
    albe <- nlspiar@par.coeffs
    Omega0 <- MV.out@Omega0
    Omega1 <- MV.out@Omega1
    Pi0 <- MV.out@Pi0
    Pi1 <- MV.out@Pi1

    Gamma <- solve(Phi0) %*% Phi1

    Xi <- list()
    Xi[[1]] <- Pi0
    for(i in 2:(hpred/s)){
      Xi[[i]] <- Gamma %*% Xi[[i-1]] + prod(albe[2,])^(i-2) * Pi1
    }

  # Forecast standard error for h=1.

    pred.se[1,] <- sqrt(sigma^2 * diag(Xi[[1]] %*% t(Xi[[1]])))

  # Forecast standard error for h=2,3,...,hpred/s.

    aux1 <- Xi[[1]] %*% t(Xi[[1]])
    for(h in 2:(hpred/s)){
      aux2 <- 0
      for(i in 2:h)
        aux2 <- aux2 + Xi[[i]] %*% t(Xi[[i]])
      pred.se[h,] <- sqrt(sigma^2 * diag(aux1 + aux2))
    }

  # Forecast.

    mus <- matrix(nlscoef[(p*s):((p+1)*s-1)], ncol=1)
    pred <- matrix(nrow=as.integer(hpred/s), ncol=s)
    Yret1 <- matrix(wts[(N-s+1):N], ncol=1)
    musx <- (1-prod(albe[2,]))^(-1) * (Omega0 + Omega1) %*% mus
    for(i in 1:as.integer(hpred/s)){
      #pred[i,] <- solve(Phi0) %*% mus + solve(Phi0) %*% Phi1 %*% Yret1
      pred[i,] <- solve(Phi0) %*% musx + Gamma %*% Yret1
      Yret1 <- matrix(pred[i,], ncol=1)
    }
  }

# Confidence intervals.

  #tl <- tu <- rep(NA, hpred); k <- 1
  tl <- tu <- matrix(NA, nrow=(hpred/s), ncol=s)
  for(i in 1:nrow(pred.se)){
    for(j in 1:s)
    {
      #tl[k] <- pred[k] - 1.96 * pred.se[i,j]
      #tu[k] <- pred[k] + 1.96 * pred.se[i,j]
      #k <- k+1
      tl[i,j] <- pred[i,j] - 1.96 * pred.se[i,j]
      tu[i,j] <- pred[i,j] + 1.96 * pred.se[i,j]
    }
  }

  t0aux <- ysooys(N, t0, N+1, s)
  t0p <- t0aux[[2]][(N+1),1:2]
  fcast <- ts(matrix(t(pred), ncol=1), frequency=s, start=t0p)
  fse <- ts(matrix(t(pred.se), ncol=1), frequency=s, start=t0p)
  ucb <- ts(matrix(t(tu)), frequency=s, start=t0p)
  lcb <- ts(matrix(t(tl)), frequency=s, start=t0p)

  new("pred.piartsm", wts=wts, p=p, hpred=hpred, fcast=fcast, fse=fse, ucb=ucb, lcb=lcb)
  ##~list(output, sigma=sigma, pred.se=pred.se)
}

plotpredpiar <- function(x){
  if (!(class(x) == "pred.piartsm"))
    stop("\n Object is not of class 'pred.piartsm'.\n")

  opar <- par(las=1)
  ts.plot(x@wts, x@lcb, x@ucb, x@fcast, lty = c(1,2,2,1), xlab="",
          col=c("black","red","red","blue"), main="Forecast and confidence intervals")
  par(opar)
}

## Set classes.

setClass("fit.partsm", representation(type="character", p="numeric", lm.ar="ANY", lm.par="ANY",
         ar.coeffs="ANY", par.coeffs="ANY"))

setClass("fit.piartsm", representation(p="numeric", nls.parameters="matrix", nls.res="numeric",
         par.coeffs="matrix", pdiff.data="ts"))

setClass("Ftest.partsm", representation(test.label="character", test.name="character", p="numeric",
         Fstat="numeric", df="numeric", pval="numeric", pvl="character", h0md="lm", hamd="ANY"))

setClass("LRur.partsm", representation(test.label="character", test.name="character", p="numeric",
         LR="numeric", LRtau="numeric", h0nls="matrix", halm="lm"))

setClass("pred.piartsm", representation(wts="ts", hpred="numeric", p="numeric", fcast="ts", fse="ts",
         ucb="ts", lcb="ts"))

setClass("MVPAR", representation(Phi0="matrix", Phi1="matrix", Gamma.eigenvalues="numeric",
         tvias="matrix"))

setClass("MVPIAR", representation(Phi0="matrix", Phi1="matrix", Gamma.eigenvalues="numeric",
         tvias="matrix", Omega0="ANY", Omega1="ANY", Pi0="ANY", Pi1="ANY"))

## Set methods for classes.

setMethod("show", "fit.partsm",
  function(object)
  {
    if(object@type == "AR"){
      ar.coeffs <- object@ar.coeffs
      model <- paste("  y_t = phi_1*y_{t-1} + phi_2*y_{t-2} + ... + phi_p*y_{t-p} + coeffs*detcomp + epsilon_t \n")
      dimnames(ar.coeffs) <- list("phi_p", paste("p=", 1:ncol(ar.coeffs), sep=""))
    }
    if(object@type == "PAR"){
      par.coeffs <- object@par.coeffs
      model <- paste("  y_t = alpha_{1,s}*y_{t-1} + alpha_{2,s}*y_{t-2} + ... + alpha_{p,s}*y_{t-p} + coeffs*detcomp + epsilon_t,  for s=1,2,...,", ncol(par.coeffs), "\n", sep="")

      dimnames(par.coeffs) <- list(paste("alpha_", 1:nrow(par.coeffs), "s", sep=""),
                                   paste("s=", 1:ncol(par.coeffs), sep=""))
    }

    cat("----\n ", paste(object@type, "model of order", object@p, ".\n\n"))
    cat(model)

    if(object@type == "AR"){
      if(object@p == 0)
        cat("----\n  None autoregressive lags (p=0). \n\n")
      if(object@p > 0){
        cat("----\n  Autoregressive coefficients. \n\n")
        print(round(ar.coeffs, 2)); cat("\n")
      }
    }
    if(object@type == "PAR"){
      if(object@p == 0)
        cat("----\n  None periodic autoregressive lags (p=0). \n\n")
      if(object@p > 0){
        cat("----\n  Autoregressive coefficients. \n\n")
        print(round(par.coeffs, 2)); cat("\n")
      }
    }
  }
)

setMethod("summary", "fit.partsm",
  function(object)
  {
    show(object)
    if(object@type == "AR")
      print(summary(object@lm.ar))
    if(object@type == "PAR")
      print(summary(object@lm.par))
  }
)

##~ A?adir slot con detcomp (.Rd), y poner en la primera l?nea de cat.
setMethod("show", "fit.piartsm",
  function(object)
  {
    model <- c(paste("  y_t - alpha_s*y_{t-1} = beta_s*(y_{t-1} - alpha_{s-1}*y_{t-2}) + \n"),
               paste("                         coeffs*detcomp + epsilon_t,  with prod(alpha_s=1) for s=1,2,...,",  ncol(object@par.coeffs), ".\n", sep=""))

    if(nrow(object@par.coeffs) == 1)
      dimnames(object@par.coeffs) <- list(c("  alpha_s"),
                                          paste("s=", 1:ncol(object@par.coeffs), sep=""))
    if(nrow(object@par.coeffs) == 2)
      dimnames(object@par.coeffs) <- list(c("  alpha_s", "  beta_s"),
                                          paste("s=", 1:ncol(object@par.coeffs), sep=""))

    cat("----\n ", paste("PIAR model of order", object@p, ".\n\n"))
    cat(model)

    cat("\n  Periodic autoregressive coefficients: \n\n")
    print(round(object@par.coeffs, 3)); cat("\n")
  }
)

setMethod("summary", "fit.piartsm",
  function(object, digits = max(3L, getOption("digits") - 3L),
           signif.stars = getOption("show.signif.stars"))
  {
    show(object)
    cat("----\n  Estimates of the non-linear model.\n")
    printCoefmat(object@nls.parameters, digits=digits, signif.stars=signif.stars)
    cat("\n----\n  Periodically differenced data.\n")
    print(round(object@pdiff.data, digits=digits)); cat("\n")
  }
)

#setMethod("plot", signature(x="fit.piartsm", y="missing"),
#  function(x){
#    opar <- par(las=1)
#    layout(matrix(c(1, 1, 2, 3), 2 , 2, byrow=TRUE))
#    plot(x@pdiff.data, main="Periodically differenced data", ylab="", xlab="")
#    bbplot(x@pdiff.data)
#    monthplot(x@pdiff.data, ylab="")
#    ##acf(x@pdiff.data, main="Autocorrelations", ylab="", na.action=na.pass)
#    ##pacf(x@pdiff.data, main="Partial autocorrelations", ylab="", na.action=na.pass)
#    par(opar)
#  }
#)
##~ monthplot, main="Seasonal subseries"

setMethod("show", "Ftest.partsm",
  function(object)
  {
    cat("----\n ", object@test.name, ".\n\n")

    if(object@test.label=="Fpar"){
      cat("  Null hypothesis: AR(", object@p, ") with the selected deterministic components.\n")
      cat("  Alternative hypothesis: PAR(", object@p, ") with the selected deterministic components.\n\n")
    }

    if(object@test.label=="Fsh"){
      ##~ HACER.
    }

    if(object@test.label=="Fnextp"){
      ifelse(substr(object@test.name, 38, 45) == "periodic", type <- "PAR", type <- "AR")
      cat("  Null hypothesis:", type, "(", object@p, ") with the selected deterministic components.\n")
      cat("  Alternative hypothesis:", type, "(", object@p+1, ") with the selected deterministic components.\n\n")
    }

    if(object@test.label=="F-PARI1" && object@test.label=="F-PARI-1"){
      switch(object@test.label,
        "F-PARI1"  = H0cat <- "long run unit root 1",
        "F-PARI-1" = H0cat <- "seasonal unit root -1")
      cat("  Null hypothesis: PAR(", object@p, ") for a series with the", H0cat, "\n")
      cat("  Alternative hypothesis: Periodically integrated AR(", object@p, "). \n\n")
    }

    cat("  F-statistic:", round(object@Fstat, 2), "on", object@df[1], "and", object@df[2], "DF,", "p-value:", object@pval, object@pvl, "\n\n")
    cat("  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n\n")
  }
)

setMethod("summary", "Ftest.partsm",
  function(object)
  {
    show(object)
    cat("----\n----\n## Fitted model for the null hypothesis.\n")
    print(summary(object@h0md))
    cat("----\n----\n## Fitted model for the alternative hypothesis.\n")
    print(summary(object@hamd))
  }
)

setMethod("show", "LRur.partsm",
  function(object)
  {
    cat("----\n ", c(object@test.name, "of order", object@p, ".\n\n"))

    cat("  Null hypothesis: PAR(", object@p, ") restricted to a unit root. \n")
    cat("  Alternative hypothesis: PAR(", object@p, "). \n\n")

    cat("  LR-statistic:", round(object@LR, 2), "\n  ---\n")
    cat("  5 and 10 per cent asymptotic critical values:\n")
    cat("    when seasonal intercepts are included: 9.24, 7.52. \n")
    cat("    when seasonal intercepts and trends are included: 12.96, 10.50. \n\n")

    cat("  LRtau-statistic:", round(object@LRtau, 2), "\n  ---\n")
    cat("  5 and 10 per cent asymptotic critical values:\n")
    cat("    when seasonal intercepts are included: -2.86, -2.57. \n")
    cat("    when seasonal intercepts and trends are included: -3.41, -3.12. \n\n")
  }
)

setMethod("summary", "LRur.partsm",
  function(object)
  {
    show(object)
    cat("----\n----\n## Fitted model for the null hypothesis.\n")
    print(object@h0nls)
    cat("----\n----\n## Fitted model for the alternative hypothesis.\n")
    print(object@halm)
  }
)

setMethod("show", "pred.piartsm",
  function(object)
  {
    yd <- as.integer(time(object@fcast))
    sd <- cycle(object@fcast)

    ysd <- rep(NA, object@hpred)
    for(i in 1:object@hpred){
      ifelse(nchar(sd[i]) == 1, sdi <- paste("0", sd[i], sep=""), sdi <- sd[i])
      ysd[i] <- paste(yd[i], sdi, sep=".")
    }

    out <- data.frame(object@fcast, object@fse, object@ucb, object@lcb)
    dimnames(out) <- list(ysd, c("fcast", "fse", "ucb", "lcb"))

    cat("----\n ", paste("  Forecasts for a PIAR model of order", object@p, ".\n\n"))
    print(out)
    cat("\n  'fcast': Forecast; 'fse': Forecast standard error; \n")
    cat("  'ucb': Upper confidence bound; 'lcb': Lower condidence bound. \n\n")
  }
)

PAR.MVrepr <- function(object)
{
   if(class(object) != "fit.partsm" || class(object) != "fit.piartms")
     stop("Object is not of class 'fit.partsm' or 'fit.piartsm'.\n")
}

setMethod("PAR.MVrepr", "fit.partsm",
  function(object)
  {
    if(is.null(object@par.coeffs))
      stop("Object is related to an AR model. A PAR model must be provided.\n")

    phi <- object@par.coeffs; s <- ncol(phi)
    Phi01 <- .PAR.MV(phi=phi, s=s)

    new("MVPAR", Phi0=Phi01$Phi0, Phi1=Phi01$Phi1, Gamma.eigenvalues=Phi01$Phi01ev, tvias=Phi01$tvias)
  }
)

setMethod("PAR.MVrepr", "fit.piartsm",
  function(object)
  {
    phi <- matrix(object@par.coeffs[1,], nrow=1)
    s <- ncol(phi); p <- nrow(object@par.coeffs)

    Phi01 <- .PAR.MV(phi=phi, s=s)
    Phi0 <- Phi01$Phi0
    Phi1 <- Phi01$Phi1

    if(p==1)
      Omega0 <- Omega1 <- Pi0 <- Pi1 <- NULL
    if(p==2)
    {
      albe <- object@par.coeffs

      Omega0 <- diag(1, s)
      for(i in 2:s)
        Omega0[i,(i-1)] <- albe[2,i]
      Omega0[3,1] <- albe[2,2]*albe[2,3]
      Omega0[4,1] <- albe[2,2]*albe[2,3]*albe[2,4]
      Omega0[4,2] <- albe[2,3]*albe[2,4]

      Omega1 <- diag(0, s)
      Omega1[1,4] <- albe[2,1]
      Omega1[1,2] <- albe[2,1]*albe[2,3]*albe[2,4]
      Omega1[1,3] <- albe[2,1]*albe[2,4]
      Omega1[2,3] <- albe[2,1]*albe[2,2]*albe[2,4]
      Omega1[2,4] <- albe[2,1]*albe[2,2]
      Omega1[3,4] <- albe[2,1]*albe[2,2]*albe[2,3]

      Pi0 <- solve(Phi0) %*% Omega0
      Pi1 <- solve(Phi0) %*% Omega1 + prod(albe[2,]) * Pi0
    }

    new("MVPIAR", Phi0=Phi01$Phi0, Phi1=Phi01$Phi1, Gamma.eigenvalues=Phi01$Phi01ev, tvias=Phi01$tvias,
        Omega0=Omega0, Omega1=Omega1, Pi0=Pi0, Pi1=Pi1)
  }
)

setMethod("show", "MVPAR",
  function(object)
  {
    Phi0 <- round(object@Phi0, 3)
    Phi1 <- round(object@Phi1, 3)
    tvias <- round(object@tvias, 3)
    Gev <- round(object@Gamma.eigenvalues, 3)

    dimnames(Phi0) <- list(rep("", nrow(Phi0)), rep("", nrow(Phi0)))
    dimnames(Phi1) <- list(rep("", nrow(Phi1)), rep("", nrow(Phi1)))
    dimnames(tvias) <- list(rep("", nrow(tvias)), rep("", nrow(tvias)))

    cat("----\n ", paste("  Multivariate representation of a PAR model.\n\n"))

    cat("\n  Phi0:\n")
    print(Phi0); cat("\n")
    cat("\n  Phi1:\n")
    print(Phi1); cat("\n")
    cat("\n  Eigen values of Gamma = Phi0^{-1} %*% Phi1:\n")
    cat(Gev, "\n")
    cat("\n  Time varing accumulation of shocks:\n")
    print(tvias); cat("\n")
  }
)

setMethod("show", "MVPIAR",
  function(object)
  {
    Phi0 <- round(object@Phi0, 3)
    Phi1 <- round(object@Phi1, 3)
    tvias <- round(object@tvias, 3)
    Gev <- round(object@Gamma.eigenvalues, 3)

    dimnames(Phi0) <- list(rep("", nrow(Phi0)), rep("", nrow(Phi0)))
    dimnames(Phi1) <- list(rep("", nrow(Phi1)), rep("", nrow(Phi1)))
    dimnames(tvias) <- list(rep("", nrow(tvias)), rep("", nrow(tvias)))

    cat("----\n ", paste("  Multivariate representation of a PIAR model.\n\n"))

    cat("\n  Phi0:\n")
    print(Phi0); cat("\n")
    cat("\n  Phi1:\n")
    print(Phi1); cat("\n")
    cat("\n  Eigen values of Gamma = Phi0^{-1} %*% Phi1:\n")
    cat(Gev, "\n")
    cat("\n  Time varing accumulation of shocks:\n")
    print(tvias); cat("\n")
  }
)

##~ Internal functions.

.PAR.MV <- function(phi, s)
{
  paux <- length(phi)/s
  P <- 1 + as.integer((paux-1)/s)

  Phi0 <- matrix(0, nrow=s, ncol=s) # ver elemento (s,1)
  for(i in 1:s){
    for(j in 1:s){
      if(i==j){ Phi0[i,j] <- 1 }
      if(j>i) { Phi0[i,j] <- 0 }
      if(j<i && (i-j) <= paux) { Phi0[i,j] <- -phi[i-j,i] }
    }
  }

  for(k in 1:P){
    #Phik <- matrix(0, nrow=s, ncol=s*P)
    Phik <- matrix(0, nrow=s, ncol=s)
    for(i in 1:s){
      for(j in 1:s)
        if(i+s*k-j <= paux) Phik[i,j] <- phi[i+s*k-j,i]
    }
    MPhi <- as.matrix(data.frame(Phi0, Phik))
  }

  Phi0 <- MPhi[,1:s]
  Phi1 <- MPhi[,(s+1):(2*s)]

  Gamma <- solve(Phi0) %*% Phi1
  Phi01ev <- eigen(Gamma, only.values = TRUE)$values
  tvias <- Gamma %*% solve(Phi0)  # Time-varing impact of accumulation of shocks

  list(Phi0=Phi0, Phi1=Phi1, Phi01ev=Phi01ev, tvias=tvias)
}



##~ monthplot, main="Seasonal subseries"

acf.ext1 <- function(wts, transf.type, perdiff.coeffs, type, lag.max, showcat, plot)
{

  switch(transf.type,
    orig    = out <- acf(wts, type=type, lag.max=lag.max,
                         plot=plot, na.action=na.exclude),
    fdiff   = out <- acf(diff(wts, lag=1), type=type, lag.max=lag.max,
                         plot=plot, na.action=na.exclude),
    sdiff   = out <- acf(diff(wts, lag=frequency(wts)), type=type, lag.max=lag.max,
                         plot=plot, na.action=na.exclude),
    fsdiff  = out <- acf(diff(diff(wts, lag=1), lag=frequency(wts)), type=type,
                         lag.max=lag.max, plot=plot, na.action=na.exclude),
    fdiffsd = { MLd <- c(NA, diff(wts, lag=1))
                SDum <- SeasDummy(wts, frequency(wts), start(wts), "alg")
                res <- lm(MLd ~ 0+SDum[,1:frequency(wts)])$residuals
                out <- acf(res, type=type, lag.max=lag.max, plot=plot, na.action=na.exclude) },

    perdiff = { ML <- ret(wts, 2)
                SDum <- SeasDummy(wts, frequency(wts), start(wts), "alg")
                filc <- rowSums(perdiff.coeffs*SDum)
                MLd <- ML[,1] - filc * ML[,2]
                out <- acf(MLd, type=type, lag.max=lag.max, plot=plot, na.action=na.exclude) },

    perdiffsd = { ML <- ret(wts, 2)
                  SDum <- SeasDummy(wts, frequency(wts), start(wts), "alg")
                  filc <- rowSums(perdiff.coeffs*SDum)
                  MLd <- ML[,1] - filc * ML[,2]
                  res <- lm(MLd ~ 0+SDum[,1:frequency(wts)])$residuals
                  out <- acf(res, type=type, lag.max=lag.max, plot=plot, na.action=na.exclude) }

  )

  se  <- 1/sqrt(out$n.used)
  tst <- out$acf/se

  ref1 <- c(1, 0.1, 0.05, 0.01, 0.001)
  ref2 <- c(" ", ".", "*", "**", "***")

  pval <- pvl <- rep(NA, length(tst))
  for(i in 1:length(tst)){
    pval[i] <- (1-pnorm(q=abs(tst[i]), mean=0, sd=1, lower.tail=TRUE, log.p=FALSE)) * 2
    pvl[i] <- ref2[length(which(ref1 >= pval[i]) == TRUE)]
  }

  if(showcat == TRUE){
    sout <- data.frame(Lag=c(0:(length(tst)-1)), acf=round(out$acf,3), pvalue=round(pval,3), pvl)
    cat("----\n  Estimated autocorrelation function for the\n")
    switch(transf.type,
      orig    = cat("  original series.\n\n"),
      fdiff   = cat("  first differences of the original series.\n\n"),
      sdiff   = cat("  seasonal differences of the original series.\n\n"),
      fsdiff  = cat("  fisrt and seasonal differences of the original series.\n\n"),
      fdiffsd = cat("  residuals of the first differences on four seasonal dummy variables.\n\n"),
      perdiff = cat("  periodic differences of the original series.\n\n"),
      perdiffsd = cat("  residuals of the periodic differences on four seasonal dummy variables.\n\n")
    )

    print(sout)
    cat("\n  s.e.=", round(se, 2), "\n")
    cat("  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
  }

  out <- data.frame(Lag=c(0:(length(tst)-1)), acf=out$acf, pvalue=pval, pvl)
  list(acf=out, se=se)
}
