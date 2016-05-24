#  ============== ConfApply  ==============================
ConfA <- function(confid)
{ 
  n <- confid[1]
  r <- confid[2]
  SL <- (r/n)*100    

  if(r > 0 & r < n)
  {  flow <- qf(1-0.05/2,(2*(n-r+1)),2*r)
     fup <- qf(1-0.05/2,(2*(r+1)),2*(n-r))
     CIlow <- 100*r/(r+(n-r+1)*flow)             
     CIup <- 100*(r+1)*fup/(n-r+(r+1)*fup)
  }

  if(r == n)
  { flow <- qf(1-0.05,2,2*n)
    CIlow <- 100*n/(n+flow)
    CIup <- 100
  }

  if(r == 0 )
  { CIlow <- 0
    fup <- qf(1-0.05,2,(2*n))
    CIup <- 100*fup/(n+fup)
  }

  if(r <0 | r > n)
  {  cat("====== out of range ==========\n")
  }

  CIwidth <- CIup-CIlow
  CIso    <- CIup-SL
  CIsu    <- SL-CIlow
   
  Confres <- rbind(n,r,SL,CIlow,CIup,CIwidth,CIso,CIsu)
  rownames(Confres) <- c("Gesamtzahl (n)","Anzahl Tote","Spontanrate",
                         "CIlow","CIup","CIwidth","CIso","CIsu")  
 return(Confres)
}

##  ============== TV.model ============================
TV.model <- function(targetSL,ntarget,sigma.bio.wahr,plot=FALSE)
{
 p.est <- targetSL/100
 CIu.TV <- 999999
 TVmax  <- 1
 while(CIu.TV > 0.5)
 { 
   TVmax <- TVmax + 1
   nTV <- ceiling( ntarget/TVmax )
   kTV <- round(p.est*nTV)
   Confres <- ConfA(c(nTV,kTV))  
   CIu.TV <- nTV * Confres["CIlow",1] 
 }

 TVmax <- TVmax - 1         
 TVmatrix <- matrix( nrow = TVmax, ncol = 5)
 colnames(TVmatrix) <- c("TV","nTV","E.TV","Var.TV","delta")
 TVmatrix[1,1] <- 1
 TVmatrix[1,2] <- ntarget
 TVmatrix[1,3] <- p.est
 TVmatrix[1,4] <- sigma.bio.wahr^2

 if(TVmax<2)
 { 
   return(c(TVmatrix[1, ])) 
 } else
 { 
   for (TV in 2:TVmax)
   {  
     nTV <- ceiling( ntarget/TV )
     kTV <- round(p.est*nTV)
     var.nTV <-  p.est*nTV*(1-p.est)
     var.p.est <-  p.est*(1-p.est)/nTV
     var.TV <- TV/(TV-1) * sigma.bio.wahr^2     
     E.TV   <- kTV/nTV       

     TVmatrix[TV,1] <- TV
     TVmatrix[TV,2] <- nTV
     TVmatrix[TV,3] <- E.TV
     TVmatrix[TV,4] <- var.TV
   }
   
    TVmatrix[ ,5] <- (TVmatrix[ ,3]-p.est)^2 + 
                     (sqrt(TVmatrix[ ,4])-sigma.bio.wahr)^2
    iopt  <- which(TVmatrix[ ,5] == min(TVmatrix[-1,5],na.rm=TRUE)) 
    TVm.opt <- TVmatrix[iopt, ]
   


  #  ---- plot ------------------------------
     if(plot)
    { xlimmin <- min(TVmatrix[ ,"E.TV"])
      xlimmax <- max(TVmatrix[ ,"E.TV"])
      deltalimx <- xlimmax - xlimmin
      ylimmin <- sqrt(min(TVmatrix[ ,"Var.TV"]))
      ylimmax <- sqrt(max(TVmatrix[ ,"Var.TV"]))

      deltalimy <- ylimmax - ylimmin
      deltalim <- max(deltalimx,deltalimy)

      xlimmax <- xlimmin + deltalim 
      ylimmax <- ylimmin + deltalim     
       
      dev.new()
      plot(TVmatrix[ ,3],sqrt(TVmatrix[ ,4]),type="p",
           xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax),xlab="E.TV",
           ylab="sqrt(Var.TV)",col="black",pch=21)
      points(TVmatrix[1,3],sqrt(TVmatrix[1,4]),type="p",col="red",pch=13)
      points(TVmatrix[iopt,3],sqrt(TVmatrix[iopt,4]),type="p",bg="green",pch=21)

     }
   return(c(TVm.opt))
  }
}


#  ============== OptNBinoApply  ==============================
OptNBinoA <- function(param,p0,alpha=0.05,beta=0.2,nmin=nmin)
{
    pA   <- param[1]
    nmax <- param[2]
    n <- nmin-1
    betaok <- FALSE

    while( (n < nmax) & !betaok)
    {
      n <- n + 1
      Test <- betaBino(n,p0,pA,alpha=alpha)
      betaok <- (Test$beta <= beta)
    }

    kkrit <- NA
    index.kkrit <- (1:(n+1))[Test$dist[ ,"k"]==Test$kkrit]
    alpha.eff <- 1-Test$dist[index.kkrit,"cdf|H0"]
    beta.eff <- Test$beta

    if (betaok)
    { 
      kkrit <- Test$kkrit
    }
    data.s <- c(betaok,n,kkrit,alpha.eff,beta.eff,pA,p0)
    names(data.s) <- c("beta ok?","n","k krit","alpha eff.","beta eff.","pA","p0")
    return(data.s)

}

#  ============== betaBino  ==============================
betaBino <- function(n,p0,pA,alpha)
{
    k <- 0:n
    cdfk.0 <- pbinom(k,n,p0)
    cdfk.A <- pbinom(k,n,pA)
    dist <- cbind(k,cdfk.0,cdfk.A)
    dimnames(dist) <- list(NULL,c("k","cdf|H0","cdf|HA"))
    kkrit <- k[cdfk.0 >= 1-alpha][1]
    beta.kkrit <- cdfk.A[k==kkrit]
  
    result.t <- list(dist,kkrit,beta.kkrit)
    names(result.t) <- c("dist","kkrit","beta")  
    return(result.t)
}


## ============ DW2s.R ===========================================
DW2s <- function(theta,d,gamma,rho,pNull=0)
{
  x    <- rep(NA,times=length(d))
  deq0 <- (d==0)
  dgt0 <- (d> 0)
  x[dgt0]  <- log(d[dgt0])

  p <- rep(NA,times=length(d))
  p[deq0] <- gamma
  p[dgt0] <- gamma + (1-gamma-rho)/(1+exp(-(theta[1]+theta[2]*x[dgt0])))

  zu.klein <-  p < pNull
  zu.gross <-  p > (1-pNull)  

  p[zu.klein] <-   pNull
  p[zu.gross] <- 1-pNull

  return(p)
}


# ============== EstM2Sn_apply.R ====================
EstM2Sn <- function(theta.ini,n,r,d,dNull,dpred,logxpred,tau,
                    SL=0,immunity=0,pNull=1.e-8,alpha=0.05)
{
  colnames(dpred) <- c("(Intercept)","d")
  npar            <- length(theta.ini)
  npred           <- nrow(dpred)
  lt              <- length(tau)
  ltau            <- matrix(tau,ncol=1)  # für apply-Funktionen
  M2Sn.theta.hut  <- rep(NA,times=npar) 
  M2Sn.Cov        <- matrix(NA,nrow=npar,ncol=npar)
  M2Sn.predp.hut  <- rep(NA,times=npred) 
  M2Sn.predp.ciu  <- rep(NA,times=npred) 
  M2Sn.predp.cio  <- rep(NA,times=npred) 
  M2Sn.logxi.hut     <- rep(NA,times=lt)
  M2Sn.logxi.hut.ciu <- rep(NA,times=lt)
  M2Sn.logxi.hut.cio <- rep(NA,times=lt)

  M2Sn.pgrid.hut  <- rep(NA,times=npred) 
  M2Sn.pgrid.ciu  <- rep(NA,times=npred) 
  M2Sn.pgrid.cio  <- rep(NA,times=npred) 

  F.M2Sn <- function (tau,immunity,SL,beta0,beta1)
            { hut <-  -(log((1-immunity-tau)/(tau-SL))+ beta0)/beta1
              return(hut)
            }

  empp      <- r/n
  existiert <- sum((0 < empp) & (empp < 1)) >= 2
  if (existiert)
  {
    deq0  <- dpred[ ,"d"] == 0
    dgt0  <- dpred[ ,"d"]  > 0
    igt0  <- (1:npred)[dgt0]
    M2Sn.SL.hut <- SL                
    M2Sn.immunity.hut   <- immunity                  
    gradtol <- 1.e-10

    M2Sn <- nlm(loglM2S,theta.ini,hessian=T,print.level=0,gradtol=gradtol,
                r=r,d=d,N=n,gamma=M2Sn.SL.hut,rho=M2Sn.immunity.hut,pNull=pNull)

    if (M2Sn$code==1) coco <- "OK, grad=0"
    if (M2Sn$code==2) coco <- "OK, delta=0"
    if (M2Sn$code==3) coco <- "Abbruch 3"
    if (M2Sn$code==4) coco <- "Zuviele Iter."
    if (M2Sn$code==5) coco <- "Abbruch 5"

    M2Sn.theta.hut <- M2Sn$estimate
    M2Sn.gradnorm  <- sqrt(M2Sn$gradient %*% M2Sn$gradient)
    M2Sn.beta0.hut <- M2Sn$estimate[1]
    M2Sn.beta1.hut <- M2Sn$estimate[2]

    NAinit       <- rep(NA,times=length(d))
    M2Sn.p.hut   <- NAinit
    M2Sn.D1beta0 <- NAinit
    M2Sn.D1beta1 <- NAinit
    M2Sn.p.hut[d==0]   <- 0
    M2Sn.D1beta0[d==0] <- 0
    M2Sn.D1beta1[d==0] <- 0
    logu <- -(M2Sn.beta0.hut + M2Sn.beta1.hut*log(d[d>0]))
   
    if (all(logu<100)) 
    {
      u <- exp(logu)
      M2Sn.p.hut[d>0]  <- DW2s(M2Sn.theta.hut,d[d>0],SL,immunity,pNull=pNull)
      M2Sn.D1beta0[d>0]    <- M2Sn.p.hut[d>0] * u / (1+u)^2
      M2Sn.D1beta1[d>0]    <- M2Sn.D1beta0[d>0] * log(d[d>0])  
      M2Sn.Inf <- matrix(NA,nrow=2,ncol=2)
      M2Sn.Inf[1,1] <-  -sum(u/(1+u)^2 * n[d>0])
      M2Sn.Inf[1,2] <-  -sum(u/(1+u)^2 * n[d>0] * log(d[d>0]) )
      M2Sn.Inf[2,1] <-   M2Sn.Inf[1,2]
      M2Sn.Inf[2,2] <-  -sum(u/(1+u)^2 * n[d>0] * (log(d[d>0])^2))

      deter <- M2Sn.Inf[1,1]*M2Sn.Inf[2,2] - M2Sn.Inf[1,2]^2
      M2Sn.Cov2 <- matrix(NA,nrow=2,ncol=2)
      if (abs(deter) > 1.e-12)
      {
        M2Sn.Cov2[1,1] <-  M2Sn.Inf[2,2]
        M2Sn.Cov2[1,2] <- -M2Sn.Inf[1,2]
        M2Sn.Cov2[2,1] <-  M2Sn.Cov2[1,2]
        M2Sn.Cov2[2,2] <-  M2Sn.Inf[1,1]
        M2Sn.Cov2      <- -M2Sn.Cov2 / deter
        cov.ok         <- TRUE
      } else
      {
        cov.ok         <- FALSE
      }
      M2Sn.Cov      <- M2Sn.Cov2

      npgrid <- npred
      pgridmax <- -(log((1-M2Sn.SL.hut-M2Sn.immunity.hut)/
                   (1-M2Sn.immunity.hut-0.01-M2Sn.SL.hut)  -1) +
                       M2Sn.beta0.hut
                   )  / M2Sn.beta1.hut 

      logdNull <- log(dNull)
      logdgrid <- data.frame("A" = rep(1,times=npgrid),
                             "B" = seq(logdNull,pgridmax,length.out=npgrid))
      colnames(logdgrid) <- c("(Intercept)","logx")
      dgrid      <- logdgrid
      dgrid[ ,2] <- exp(logdgrid[ ,2])
      names(dgrid) <- c("(Intercept)","d")
      deq0  <- dgrid[ ,"d"] == 0
      dgt0  <- dgrid[ ,"d"]  > 0
      igt0  <- (1:npgrid)[dgt0]

      M2Sn.logxbeta.hut <- rep(NA,times=npgrid)
      M2Sn.logxbeta.var <- rep(NA,times=npgrid)
      M2Sn.logxbeta.hut[deq0] <- M2Sn.SL.hut
      M2Sn.logxbeta.var[deq0] <- 0 
      M2Sn.pgrid.hut[deq0] <- 0
      M2Sn.logxbeta.hut[dgt0] <- as.matrix(logdgrid[dgt0, ]) %*% 
                                 matrix(M2Sn$estimate,ncol=1)
      M2Sn.pgrid.hut[dgt0] <- M2Sn.SL.hut +
                             (1-M2Sn.SL.hut-M2Sn.immunity.hut)/
                             (1+exp(-M2Sn.logxbeta.hut[dgt0]))

      M2Sn.logxi.hut <- apply(ltau,MARGIN=1,FUN=F.M2Sn,immunity=M2Sn.immunity.hut,
                              SL=M2Sn.SL.hut,beta0=M2Sn.beta0.hut,
                              beta1=M2Sn.beta1.hut)
      zalpha <- qnorm(1-alpha/2)

      if (cov.ok)
      {
        for (i in 1:npgrid)
        { M2Sn.logxbeta.var[i] <- as.matrix(logdgrid[i, ]) %*% 
                                  M2Sn.Cov2 %*%
                                  t(as.matrix(logdgrid[i, ]))
        }
        M2Sn.logxbeta.se <- sqrt(M2Sn.logxbeta.var)
        M2Sn.pgrid.ciu <- M2Sn.logxbeta.hut - zalpha*M2Sn.logxbeta.se
        M2Sn.pgrid.cio <- M2Sn.logxbeta.hut + zalpha*M2Sn.logxbeta.se
        M2Sn.pgrid.ciu <- M2Sn.SL.hut +
                         (1-M2Sn.SL.hut-M2Sn.immunity.hut)/
                         (1+exp(-M2Sn.pgrid.ciu))
        M2Sn.pgrid.cio <- M2Sn.SL.hut +
                         (1-M2Sn.SL.hut-M2Sn.immunity.hut)/
                         (1+exp(-M2Sn.pgrid.cio))
        M2Sn.logxi.hut.ciu <- apply(ltau,MARGIN=1,FUN=InterpolxA,
                                    x=logdgrid$logx,y=M2Sn.pgrid.cio)
        M2Sn.logxi.hut.cio <- apply(ltau,MARGIN=1,FUN=InterpolxA,
                                    x=logdgrid$logx,y=M2Sn.pgrid.ciu)
      }                                
    }  else  coco <- "Divergenz" 
  }   else  coco <- "Nichtexistent"

  if (coco == "OK, grad=0")    coconum = 1
  if (coco == "OK, delta=0")   coconum = 2
  if (coco == "Nichtexistent") coconum = 3
  if (coco == "Divergenz")     coconum = 4
  if (coco == "Zuviele Iter.") coconum = 5
  if (coco == "Abbruch 3")     coconum = 6
  if (coco == "Abbruch 5")     coconum = 7  

  Result <- list(coco,coconum,M2Sn.theta.hut,M2Sn.gradnorm,M2Sn.Cov,
                 M2Sn.pgrid.hut,M2Sn.pgrid.ciu,M2Sn.pgrid.cio,
                 M2Sn.logxi.hut,M2Sn.logxi.hut.ciu,M2Sn.logxi.hut.cio,
                 dgrid,logdgrid)
  names(Result) <- c("CondCode","CondCodeNum","theta.hut","GradNorm","Cov",
                     "predp.hut","predp.ciu","predp.cio",
                     "logxi.hut","logxi.hut.ciu","logxi.hut.cio",
                     "dgrid","logdgrid")
  return(Result)
}


# ============== loglM2S.R ==================== 
loglM2S <- function(theta,N,r,d,gamma,rho,pNull)
{
  p <- DW2s(theta,d,gamma,rho,pNull=pNull)
  logli <- r * log(p) + (N-r)* log(1-p)
  logl <- sum(logli)
  return(-logl)
}


# ============== Interpolxapply.R ====================  
InterpolxA <- function(ytarget,x,y)
{  
  xx <- x[!is.na(x) & !is.na(y)]
  yy <- y[!is.na(x) & !is.na(y)]
  xtarget <- NA

  if (any(is.na(xx)))
  {arg.ok <- FALSE  }
  if (any(is.na(yy)))
  {arg.ok <- FALSE  }
  arg.ok <- TRUE
  if (is.na(ytarget)) 
  {arg.ok <- FALSE  }
  n <- length(yy)
  if (n <= 1) 
  {arg.ok <- FALSE  }
  if (length(xx) != n) 
  {arg.ok <- FALSE  }

  if (arg.ok)
  {
    if ( !(yy[1] <= ytarget & ytarget <= yy[n]) &
         !(yy[n] <= ytarget & ytarget <= yy[1]) )
    { arg.ok <- FALSE    }
  }

  if (arg.ok)
  {
    for (i in 2:n)
    {
      if ( (yy[i-1] <= ytarget & ytarget < yy[i  ]) |
           (yy[i  ] <= ytarget & ytarget < yy[i-1]) )
      {
        xtarget <- xx[i-1] + (xx[i]-xx[i-1])*(ytarget-yy[i-1])/(yy[i]-yy[i-1])
        break
      }
    }
    if (ytarget == yy[n]) xtarget <- xx[n]
  }
  return(xtarget)
}


# ============== dosedistrib.R ==================
plan <- function(D,ntarget)
{  V <- data.frame(matrix(NA,nrow=nrow(D),ncol=max(ntarget)))
   colnames(V) <- paste("Conc",1:max(ntarget))
   nt <- length(ntarget)

   if (nt != nrow(D)) cat("******  plan: internal error 01  ******\n")
   for (i in 1:nt)
   {
     if (ntarget[i] > 0)
     {  nrest <- ntarget[i]
        if (nrest %% 2 == 0)      
        { if (ntarget[i] == 2)
          { j <- 1
            V[i,j] <- D[i,"EC.target"] + D[i,"FCIU95"]
            j <- 2
            V[i,j] <- D[i,"EC.target"] + D[i,"FCIO95"]
          } else
          { j <- 0
            faktor <- 1
            border <- (ntarget[i]/2)
            while(nrest > 0)
            { j <- j+1
              V[i,j] <- D[i,"EC.target"] + (faktor/border) * D[i,"FCIU99"]
              j <- j+1
              V[i,j] <- D[i,"EC.target"] + (faktor/border) * D[i,"FCIO99"]
              faktor <- faktor + 1
              nrest  <- nrest - 2
            }
          }
        }  else
        { j <- 1 
          V[i,j] <- D[i,"EC.target"]
          nrest <- nrest - 1
          if (nrest ==2)
          { j <- j+1
            V[i,j] <- D[i,"EC.target"] + D[i,"FCIU95"]
            j <- j+1
            V[i,j] <- D[i,"EC.target"] + D[i,"FCIO95"]
          } else
          { faktor <- 1
            border <- (ntarget[i]-1)/2
            while(nrest > 0)
            { j <- j+1
              V[i,j] <- D[i,"EC.target"] + (faktor/border) * D[i,"FCIU99"]
              j <- j+1
              V[i,j] <- D[i,"EC.target"] + (faktor/border) * D[i,"FCIO99"]
              faktor <- faktor + 1
              nrest  <- nrest - 2
            }
          }
        }
     }
   }
   return(V)
}




