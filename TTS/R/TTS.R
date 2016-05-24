##' Time-Temperature Superposition (TTS) analysis
##' 
##' The Master Curve at a specific temperature is estimated using
##' Time-Temperature Superposition (TTS) procedures.  The Master Curve means
##' the variation of a specific viscoelastic property of the selected material
##' depending on time or frequency. TTS procedures provide the viscoelastic
##' property variation at the selected temperature in a wider interval of time
##' or frequency than in the experimental case. The Master Curve is estimated
##' for each selected reference temperaure using TTS procedures. Three TTS
##' methodologies are implemented in this package: the two wider used methods,
##' Arrhenius based methods and WLF, and the newer methodology based on
##' derivatives procedure. The Master Curve is smoothed by B-splines basis.
##' 
##' El New method for estimating shift factors in time-temperatura
##' superposition models (Naya et al., 2013) opens the possibility to perform
##' the TTS function. The horizontal and vertical shifts are calculated.
##' Namely, the different methods are differenciated due to the expression for
##' estimating the horizontal shifts, aT. The "derivated" method is based on
##' the application of horizontal shifts to the moduli derivatives (depending
##' on the frequency) and thus obtaining the Master Curve at the selected
##' temperature:\cr
##' 
##' (dE')/dx(x+aT) -> (dE')/dx(x)\cr
##' 
##' WLF method is defined by the parametric expression:\cr
##' 
##' Log10(aT)=-C1*(T-To)/(C2+(T-To))\cr
##' 
##' Where C1 and C2 are constants to be estimated, T is the temperature and To
##' the reference temperature.\cr
##' 
##' Arrhenius method is defined by the parametric expression:\cr
##' 
##' Log10(aT)=Ea*((1/T)-(1/To))*log10(2.718282)/R
##' 
##' Where Ea is the activation energy, R = 8.314 J/mol is the universal gas
##' constant, T is the absolute temperature (Kelvin degrees), and To the
##' reference temperature (Celsius degrees).\cr
##' 
##' The vertical shifts, bT, are calculated taking into account the vertical
##' distance between the moduli curves.
##' 
##' @param x Matrix or data frame composed of three columns: a numeric column
##' vector with the experimental frequencies (in logarithmic scale, base-ten),
##' the modulus (E' or G') base-ten logarithm vector and, finally the
##' corresponding temperatures vector.
##' @param reference.temperature Value of the selected reference temperatura at
##' which the Master Curve of modulus will be obtained, the default value of
##' temperature is 150.
##' @param n Number of partitions in the frequency domain in order to fit the
##' B-spline basis. The default value is 100.
##' @param nB Number of bootstrap replicates to estimate confidence intervals
##' of master curve fitting. The default is 100.
##' @param method A string vector composed of one of the following options:
##' "derived" (by default), "WLF" and "Arrhenius".
##' @return The function returns a list composed of the following outputs:
##' 
##' \item{data}{Input experimental data.} \item{aT}{Numerical vector of
##' horizontal shifts between the modulus curves.} \item{bT}{Numerical vector
##' of vertical shifts between the modulus curves.} \item{TTS.data}{Master
##' Curve Data frame defined by three columns: log10frequency, log10module and
##' temperature.} \item{ref.temp}{Reference temperature value.}
##' \item{TTS.gam}{Data frame of the Generalized Additive Model with B-splines
##' (GAM) estimate of the Master Curve. It contains two columns: frequency and
##' Prediction.} \item{I.lower}{Lower limit of bootstrap confidence interval
##' corresponding to the estimated B-splines Master Curve.}
##' \item{I.upper}{Upper limit of bootstrap confidence interval corresponding
##' to the estimated B-splines Master Curve.} \item{residuals}{Residuals
##' corresponding to the GAM with B-splines Master Curve fitting.}
##' @author Antonio Meneses \email{antoniomenesesfreire@@hotmail.com}, Salvador
##' Naya \email{salva@@udc.es} and Javier Tarrio-Saavedra
##' \email{jtarrio@@udc.es}
##' @references Naya, S., Meneses A., Tarrio-Saavedra, J., Artiaga R.,
##' Lopez-Beceiro, J. and Gracia-Fernandez C. (2013) New method for estimating
##' shift factors in time-temperatura superposition models. Journal of Thermal
##' Analysis and Calorimetry. ISSN 1388-6150. DOI 10.1007/s10973-013-3193-1.\cr
##' 
##' Williams, M. L. (1964) Structural analysis of Viscoelastic materials. AIAA
##' Journal, 785-808.\cr
##' 
##' Artiaga R., Garcia A. Fundamentals of DMA. In: 'Thermal analysis.
##' Fundamentals and applications to material characterization' (ed.: Artiaga
##' R.) Publicaciones de la Universidade da Coruna, A Coruna, Spain, 183-206
##' (2005).\cr
##' @keywords TTS
##' @examples
##' 
## library(TTS)
##' ## Polycarbonate dataset
##' data(PC)
##' x=PC
##' ## TTS function applied to polycarbonate.
##' derive=TTS(x,reference.temperature=150, method=c("derived","WLF","Arrhenius"))
##' names(derive)
##' ##[1] "data"      "aT"        "bT"        "TTS.data"  "ref.temp"  "TTS.gam"
##' ##[7] "I.lower"   "I.upper"   "residuals"
##' ## Horizontal shifts vector of modulus versus frequency curves.
##' derive$aT
##' ## Reference temperature
##' derive$ref.temp
##' 
##' @export TTS
TTS <-
function(x,reference.temperature=150,n=100,nB=100,method=c("derived","WLF","Arrhenius"))
  {
    BS1 <- c();BS2 <- c();BS3 <- c()
    for(i in 1:length(unique(x[,3])))
        {
          BS1 <- c(BS1,spline(x[,1][x[,3]==unique(x[,3])[i]],x[,2][x[,3]==unique(x[,3])[i]],n)$x)
          BS2 <- c(BS2,spline(x[,1][x[,3]==unique(x[,3])[i]],x[,2][x[,3]==unique(x[,3])[i]],n)$y)
          BS3 <- c(BS3,rep(unique(x[,3])[i],n))
        }
    nCOL <- n*length(unique(x[,3]))
    BS <- matrix(c(BS1,BS2,BS3),nrow=nCOL,ncol=3)
    D2 <- c()
    for(i in 1:length(unique(x[,3])))
        {
          D2 <- c(D2,D1ss(spline(x[,1][x[,3]==unique(x[,3])[i]],x[,2][x[,3]==unique(x[,3])[i]],n)$x,
                         spline(x[,1][x[,3]==unique(x[,3])[i]],x[,2][x[,3]==unique(x[,3])[i]],n)$y))
        }
    D <- matrix(c(BS1,D2,BS3),nrow=nCOL,ncol=3)
    j <- c()
    for(k in 1:(length(unique(D[,3]))-1))
    {
      MM <- c()
      for(i in 1:n)
        {
        LL <- mean(D[,2][D[,3]==unique(D[,3])[k]][1:(n-i+1)]-D[,2][D[,3]==unique(D[,3])[k+1]][i:n],trim=0.1)
        MM <- c(MM,LL)
        }
      for(i in 2:n)  if(MM[i-1]*MM[i]<0) j <- c(j,i)
    }
    LH <- c()
    for(i in 1:(length(unique(D[,3]))-1))
    {
      lh <- ( D[,1][D[,3]==unique(D[,3])[i]][j[i]] ) - ( D[,1][D[,3]==unique(D[,3])[i]][1] )
      LH <- c(LH,lh)
    }
    AH <- c()
    for(i in 1:(length(LH)-1))
        {
          AH <- c(AH,c(unique(D[,3])[(i+1)],cumsum(LH[i:1])[i:1],0,cumsum(LH[(i+1):length(LH)])))
        }
    LHH <- c(c(unique(D[,3])[1],0,cumsum(LH[1:length(LH)])),AH,c(unique(D[,3])[(length(LH)+1)],cumsum(LH[length(LH):1])[length(LH):1],0))
    MATRIX.LH <- matrix(LHH,nrow=(length(LH)+2),ncol=(length(LH)+1))
    for(i in 1:(length(LH)+1))
        {
          if(reference.temperature==MATRIX.LH[1,i]) Mref <- MATRIX.LH[2:(length(LH)+2),i]
         }
    if(x[1,2]<=x[length(x[,2][x[,3]==unique(x[,3])[1]]),2])
       { m <- 1
    }else m <- 2
    MrefIS <- c((-1)^(m+1)*Mref[1:which(Mref==0)],(-1)^m*Mref[which(Mref==0):length(Mref)][-1])
    MrefISH <- c()
    for(i in 1:length(MrefIS)) MrefISH <- c(MrefISH,rep(MrefIS[i],n))
    BSHH <- BS; BSHH[,1] <- BS[,1]+MrefISH
    BSHV <- function(BSH,m,reference.temperature)
     {
          if(m==1)
          {
            LV <- c()
            for(i in 1:length(LH))
             {
            L1 <- max(which(BSH[,1][BSH[,3]==unique(BSH[,3])[i+1]]<=BSH[,1][BSH[,3]==unique(BSH[,3])[i]][n^(m-1)]))
            L2 <- min(which(BSH[,1][BSH[,3]==unique(BSH[,3])[i+1]][n^(m*(2-m))]<=BSH[,1][BSH[,3]==unique(BSH[,3])[i]]))
            L1S <- spline(BSH[,1][BSH[,3]==unique(BSH[,3])[i]][1:L2],BSH[,2][BSH[,3]==unique(BSH[,3])[i]][1:L2],n)
            L1I <- spline(BSH[,1][BSH[,3]==unique(BSH[,3])[i+1]][L1:n],BSH[,2][BSH[,3]==unique(BSH[,3])[i+1]][L1:n],n)
            LV <- c(LV,mean(L1S$y-L1I$y))
             }
          }else
          {
            LV <- c()
            for(i in 1:(length(unique(BSH[,3]))-1))
             {
            L1 <- max(which(BSH[,1][BSH[,3]==unique(BSH[,3])[i+1]]<=BSH[,1][BSH[,3]==unique(BSH[,3])[i]][n^(m-1)]))
            L2 <- min(which(BSH[,1][BSH[,3]==unique(BSH[,3])[i+1]][n^(m*(2-m))]<=BSH[,1][BSH[,3]==unique(BSH[,3])[i]]))
            L1S <- spline(BSH[,1][BSH[,3]==unique(BSH[,3])[i]][1:L2],BSH[,2][BSH[,3]==unique(BSH[,3])[i]][1:L2],n)
            L1I <- spline(BSH[,1][BSH[,3]==unique(BSH[,3])[i+1]][L1:n],BSH[,2][BSH[,3]==unique(BSH[,3])[i+1]][L1:n],n)
            LV <- c(LV,mean(L1S$y-L1I$y))
            }
          }
          AV <- c()
          for(i in 1:(length(LV)-1))
              {
                AV <- c(AV,c(unique(BSH[,3])[(i+1)],cumsum(LV[i:1])[i:1],0,cumsum(LV[(i+1):length(LV)])))
              }
          LVV <- c(c(unique(BSH[,3])[1],0,cumsum(LV[1:length(LV)])),AV,c(unique(BSH[,3])[(length(LV)+1)],cumsum(LV[length(LV):1])[length(LV):1],0))
          MATRIX.LV <- matrix(LVV,nrow=(length(LV)+2),ncol=(length(LV)+1))
          for(i in 1:(length(LV)+1))
              {
                if(reference.temperature==MATRIX.LV[1,i]) Vref <- MATRIX.LV[2:(length(LV)+2),i]
               }
          VrefIS <- c(-Vref[1:which(Vref==0)],Vref[which(Vref==0):length(Vref)][-1])
          VrefISV <- c()
          for(i in 1:length(VrefIS)) VrefISV <- c(VrefISV,rep(VrefIS[i],n))
          BSV <- BSH; BSV[,2] <- BSH[,2]+VrefISV
       list(BSV1=BSV,VrefIS1=VrefIS)
     }
    VrefIS <- BSHV(BSHH,m,reference.temperature)$VrefIS1
    VrefISX <- c()
    for(i in 1:length(VrefIS)) VrefISX <- c(VrefISX,rep(VrefIS[i],length(x[,2][x[,3]==unique(x[,3])[i]])))
    MrefISX <- c()
    for(i in 1:length(MrefIS)) MrefISX <- c(MrefISX,rep(MrefIS[i],length(x[,1][x[,3]==unique(x[,3])[i]])))
    XX <- x; XX[,1] <- x[,1]+MrefISX; XX[,2] <- x[,2]+ VrefISX
    ii <- order(XX[,1])
    X.gam <- gam(XX[,2][ii]~s(XX[,1][ii]))
    X.pred <- predict(X.gam,type="response")
    X.resid <- resid(X.gam)
    ICB <- function(X,B=nB)
     {
            ii <- order(X[,1])
            x.PC <- X[,1][ii]; y.PC <- X[,2][ii]
            quantiles.PC <- function(x.PC,probs=c(0.025,0.975))
                               {
                                quantile(x.PC,probs)
                               }
            ajuste_gam.PC <- function(y.PC)
                                  {
                                   predict.gam(gam(y.PC~s(x.PC)),type='response')
                                  }
            muhat.PC <- predict.gam(gam(y.PC~s(x.PC)),type="response")
            e.PC <- y.PC-muhat.PC; e.PC <- e.PC-mean(e.PC); sigma.PC <- sd(e.PC)
            nn <- length(x.PC)
            yboot.PC <- matrix(nrow=nn,ncol=B)
            for (icol in 1:B)
                {
                 yboot.PC[,icol] <- muhat.PC+rnorm(nn)*sigma.PC
                }
            muhatboot.PC <- apply(yboot.PC,2,ajuste_gam.PC)
            ic.PC <- apply(muhatboot.PC,1,quantiles.PC)
            lib <- ic.PC[1,];lsb <- ic.PC[2,]
            list(X.PC=x.PC,LIB=lib,LSB=lsb)
    }
    WLF.x <- unique(x[,3])[unique(x[,3])!=reference.temperature]
    WLF.Lh <- MrefIS[MrefIS!=0]
          WLF2 <- 1/WLF.Lh; WLF1 <- 1/(WLF.x-reference.temperature)
          ajuste <- lm(WLF2 ~ WLF1)
          WLF.a <- -1/ajuste$coef[1]
          WLF.b <- -WLF.a*ajuste$coef[2]
    log.WLF.at <- -WLF.a*(unique(x[,3])-reference.temperature)/(WLF.b+unique(x[,3])-reference.temperature)
    WLF.MrefISH <- c()
    for(i in 1:length(MrefIS)) WLF.MrefISH <- c(WLF.MrefISH,rep(log.WLF.at[i],n))
    WLF.BSH <- BS; WLF.BSH[,1] <- BS[,1]+WLF.MrefISH
    WLF.VrefIS <- BSHV(WLF.BSH,m,reference.temperature)$VrefIS1
    WLF.VrefISX <- c()
    for(i in 1:length(WLF.VrefIS)) WLF.VrefISX <- c(WLF.VrefISX,rep(WLF.VrefIS[i],length(x[,2][x[,3]==unique(x[,3])[i]])))
    WLF.MrefISX <- c()
    for(i in 1:length(log.WLF.at)) WLF.MrefISX <- c(WLF.MrefISX,rep(log.WLF.at[i],length(x[,1][x[,3]==unique(x[,3])[i]])))
    WLF.X <- x; WLF.X[,1] <- x[,1]+WLF.MrefISX; WLF.X[,2] <- x[,2]+ WLF.VrefISX
    iii <- order(WLF.X[,1])
    WLF.X.gam <- gam(WLF.X[,2][iii]~s(WLF.X[,1][iii]))
    WLF.X.pred <- predict(WLF.X.gam,type="response")
    WLF.X.resid <- resid(WLF.X.gam)
    AR.x <- unique(x[,3])+ 273.15
    AR <- MrefIS
          AR.z <- 1/AR.x
          ajusteAR <- lm(AR~AR.z)
          AR.A <- ajusteAR$coef[1]
          AR.BB <- ajusteAR$coef[2]
    R <- 8.314
    Ea <- AR.BB*R/log10(exp(1))
    log.AR.at <- Ea*(1/AR.x-1/(reference.temperature+273.15))/R * log10(exp(1))
    AR.MrefISH <- c()
    for(i in 1:length(MrefIS)) AR.MrefISH <- c(AR.MrefISH,rep(log.AR.at[i],n))
    AR.BSH <- BS; AR.BSH[,1] <- BS[,1]+AR.MrefISH
    AR.VrefIS <- BSHV(AR.BSH,m,reference.temperature)$VrefIS1
    AR.VrefISX <- c()
    for(i in 1:length(AR.VrefIS)) AR.VrefISX <- c(AR.VrefISX,rep(AR.VrefIS[i],length(x[,2][x[,3]==unique(x[,3])[i]])))
    AR.MrefISX <- c()
    for(i in 1:length(log.AR.at)) AR.MrefISX <- c(AR.MrefISX,rep(log.AR.at[i],length(x[,1][x[,3]==unique(x[,3])[i]])))
    AR.X <- x; AR.X[,1] <- x[,1]+AR.MrefISX; AR.X[,2] <- x[,2]+ AR.VrefISX
    iv <- order(AR.X[,1])
    AR.X.gam <- gam(AR.X[,2][iv]~s(AR.X[,1][iv]))
    AR.X.pred <- predict(AR.X.gam,type="response")
    AR.X.resid <- resid(AR.X.gam)
    ICB.X <- ICB(XX)
    deriv.icb.x <- ICB.X$X.PC
    d.lower1 <- ICB.X$LIB
    d.lower2 <- ICB.X$LSB
    ICB.WLF <- ICB(WLF.X)
    WLF.icb.x <- ICB.WLF$X.PC
    W.lower1 <- ICB.WLF$LIB
    W.lower2 <- ICB.WLF$LSB
    ICB.AR <- ICB(AR.X)
    AR.icb.x <- ICB.AR$X.PC
    A.lower1 <- ICB.AR$LIB
    A.lower2 <- ICB.AR$LSB
    der.gam <- data.frame(frequency=deriv.icb.x,prediction=X.pred)
    W.gam <- data.frame(frequency=WLF.icb.x,prediction=WLF.X.pred)
    A.gam <- data.frame(frequency=AR.icb.x,prediction=AR.X.pred)
      if(method[1]=="derived") DERIV <- list(data=x,aT=MrefIS,bT=VrefIS,TTS.data=XX,
                                     ref.temp=reference.temperature, TTS.gam=der.gam,
                                     I.lower=d.lower1,I.upper=d.lower2,residuals=X.resid)
      else if(method[1]=="WLF") W_L_F <- list(data=x,aT=log.WLF.at,bT=WLF.VrefIS,TTS.data=WLF.X,
                                ref.temp=reference.temperature, TTS.gam=W.gam,
                                     I.lower=W.lower1,I.upper=W.lower2,residuals=WLF.X.resid)
      else if(method[1]=="Arrhenius") ARRHEN <- list(data=x,aT=log.AR.at,bT=AR.VrefIS,
                                      TTS.data=AR.X,ref.temp=reference.temperature,
                                      TTS.gam=A.gam,I.lower=A.lower1,I.upper=A.lower2,
                                      residuals=AR.X.resid)
}
