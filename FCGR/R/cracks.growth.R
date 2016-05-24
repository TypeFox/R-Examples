cracks.growth <-
function(x, aF, T_c, method=c("SEP-lme_bkde","SEP-lme_kde","PB-nlme"),
                        nMC=5000, nBKDE=5000, nKDE=5000)
{
  a0=c()
  for(i in 1:length(unique(x[,3])))
   {
      a0=c(a0,rep(subset(x[,2],x[,3]==unique(x[,3])[i])[1],
                  length(x[,3][x[,3]==unique(x[,3])[i]])))
      Length=length(rep(subset(x[,2],x[,3]==unique(x[,3])[i])[1],
                  length(x[,3][x[,3]==unique(x[,3])[i]])))
    }
  x[,4]= a0
  colnames(x)=c("cycles","cracks","sample","a0")

 
        
               Paris.F <- function(cycles, C, m, a0, F.cte=1, S = 1)
                     (a0^(1 - m / 2) + (1 - m / 2) * C * F.cte^m * S^m * pi^(m / 2) * cycles)^(2 / (2 - m))
        
                     
              DerivInit <- function(data, cycles, y){
                    xx <- data$cycles
                   y <- data$cracks
                    deriv.ii=sfsmisc::D1D2(xx[-1],y[-1])
                    derivada <-deriv.ii$D1
                    grieta.media <- y[-1]
                    data.deriv <- data.frame(xx[-1], derivada, grieta.media)
                    colnames(data.deriv)[1] <- "cycles"
                    row.names(data.deriv) <- NULL
                    data.deriv
                }
          
              ParisInit <- function(mCall, data, LHS){
                     data.derivative <- DerivInit(data = data, cycles = mCall[["cycles"]], y = LHS)
                     fitted <- lm(log(derivada) ~ log(sqrt(pi * grieta.media)), data = data.derivative)
                     result <- c(exp(coefficients(fitted)[1]), coefficients(fitted)[2])
                     names(result) <- mCall[c("C", "m")]
                     result
                }

       SSParis <- selfStart(Paris.F, initial = ParisInit, parameters = c("C", "m"))

       environment(SSParis)<- .GlobalEnv

       assign("SSParis", SSParis, envir = environment(SSParis))   
            
       x <<- x

  model <- nlsList(cracks ~ SSParis(cycles, C, m, a0,
        F.cte =1, S=1)| sample, na.action = na.omit,data=x)                   
  model.nlme=nlme(model, random = C + m ~ 1, na.action = na.omit)
  COEF=coef(model.nlme); CC=mean(COEF[,1]); mm=mean(COEF[,2]); COV=var(COEF)
  sigma.r=model.nlme$sigma
  resid1=resid(model.nlme)
  ORDER.COEF=as.numeric(row.names(coef(model.nlme)))
  ii=order(ORDER.COEF)
  ORDER1.COEF=data.frame(C=COEF[,1][ii],m=COEF[,2][ii])
  grieta.P.nlme=list()
  for(i in 1:length(unique(x[,3]))){
    C=ORDER1.COEF[i,1];m=ORDER1.COEF[i,2]
    F.cte=1 ; S=1
    t=subset(x[,1],x[,3]==unique(x[,3])[i])
    a00=subset(x[,4],x[,3]==unique(x[,3])[i])[1]
    grieta.P.nlme[[i]]=(a00^(1 - m / 2) + (1 - m / 2) * C * F.cte^m * S^m * pi^(m / 2) * t)^(2 / (2 - m))
  }
  growt.est1=c()
  for(i in 1:length(unique(x[,3]))) growt.est1=c(growt.est1,grieta.P.nlme[[i]])
  Grieta.P.nlme=data.frame(time=x[,1],growt.est=growt.est1,sample=x[,3])
  T=function(a_0,aF,C,m){
      ( 1- (aF/a_0)^(-(m-2)/2))/((a_0^(m/2 -1))*( m / 2 -1) * C * pi^(m / 2))
      }

      T.M=c()
  for(i in 1:length(unique(x[,3]))){
      a00=subset(x[,4],x[,3]==unique(x[,3])[i])[1]
      T1=T(a_0=a00,aF,C=ORDER1.COEF[i,1],m=ORDER1.COEF[i,2])
      T.M=c(T.M,T1)
  }
  L2=seq(min(x[,1]),max(T.M),abs(max(T.M)-min(x[,1]))/1000)
  pred.P.nlme=list()
  for(i in 1:length(unique(x[,3]))){
      C=ORDER1.COEF[i,1];m=ORDER1.COEF[i,2];F.cte=1 ; S=1
      tt=L2
      a00=subset(x[,4],x[,3]==unique(x[,3])[i])[1]
      g.P.nlme=(a00^(1 - m / 2) + (1 - m / 2) * C * F.cte^m * S^m * pi^(m / 2) * tt)^(2 / (2 - m))
      pred.P.nlme[[i]]=data.frame(tt,g.P.nlme)
  }
  growt.pred1=c()
  time1=c()
  sample1=c()
  for(i in 1:length(unique(x[,3]))){
      growt.pred1=c(growt.pred1,pred.P.nlme[[i]]$g.P.nlme)
      time1=c(time1,L2)
      sample1=c(sample1,rep(i,length(L2)))
  }
  Grieta.P.nlme.pred=data.frame(time=time1,growt.pred=growt.pred1,
                                    sample=sample1)
  N=numeric(length(unique(x[,3])))
  for(i in 1:length(unique(x[,3]))){
      N[i]=min(which(round(pred.P.nlme[[i]][,2],2)==aF))
  }
  T.FALLO=c()
  for(i in 1:length(unique(x[,3]))) T.FALLO=c(T.FALLO,pred.P.nlme[[i]][,1][N[i]])
  set.seed(22)
  T=function(a_0,aF,C,m){
      ( 1- (aF/a_0)^(-(m-2)/2))/((a_0^(m/2 -1))*( m / 2 -1) * C * pi^(m / 2))
  }
  t_fallo=numeric(nMC)
  for(i in 1:length(t_fallo))
  {
        param= MASS::mvrnorm(mu = c(CC,mm), Sigma = COV)
        if(param[1]>0 & param[2]>0)
        {
            a_0=min(x[,2])
            t_fallo[i]=T(a_0,aF,C=param[1],m=param[2] )
        }
  }
  BOXP=boxplot.stats(t_fallo,coef=3)$stats
  t.fallo.2=t_fallo[0<t_fallo & t_fallo<=max(BOXP)]
  F=quantile(t.fallo.2, probs = seq(0, 1, 0.01))
  F1=as.vector(F)
  F1.F=seq(0,1,0.01)
  Fe.1=c()
  for(i in 1:length(unique(x[,3]))) Fe.1[i]=(i-0.5)/length(unique(x[,3]))
  Fe.1=round(Fe.1,3)
  iii=order(T.FALLO);T.FALLO.1=T.FALLO[iii]
  Dist.T.PB_nlme=data.frame(time=T.FALLO.1,Fe=Fe.1)
  Dist.PB_nlme=data.frame(time=F1,F=F1.F)
        
  x.gamm=data.frame(cycles=x$cycles,cracks=x$cracks,sample=as.factor(x$sample),a0=x$a0)
  ajuste.gamm =mgcv::gamm((cracks-a0)~s(cycles,by=sample)+sample,random=list(sample=~1),data=x.gamm)
  resid2=resid(ajuste.gamm$lme)
  pred.lme=predict(ajuste.gamm$lme,type="response")
  be.pred=data.frame(cycles=x$cycles, pred=as.vector(pred.lme),
                     cracks=x$cracks, sample=x$sample)
  spline.pr=list()
  for(i in 1:length(unique(be.pred$sample))){
      spline.pr[[i]] <- pspline::smooth.Pspline(be.pred$cycles[be.pred$sample==unique(be.pred$sample)[i]],
                                     be.pred$pred[be.pred$sample==unique(be.pred$sample)[i]],
                                     norder = 2,method=3)
  }
  pred.pr=list()
  for(i in 1:length(unique(be.pred$sample))){
      pred.pr[[i]] <- as.vector(predict(spline.pr[[i]],
                      be.pred$cycles[be.pred$sample==unique(be.pred$sample)[i]],
                      nderiv = 1))
  }
  coef.lme=list()
  for(i in 1:length(unique(be.pred$sample))){
      lin.mod <- lm(log(pred.pr[[i]][-1]) ~ log(sqrt(pi * be.pred$cracks[be.pred$sample==unique(be.pred$sample)[i]][-1])))
      coef.lme[[i]] <- coefficients(lin.mod)
  }
  C.lme=c()
  m.lme=c()
  for(i in 1:length(unique(be.pred$sample))){
      C.lme[i]=exp(as.vector(coef.lme[[i]][1]))
      m.lme[i]=as.vector(coef.lme[[i]][2])
  }
  
  COEF.lme=data.frame(C=C.lme,m=m.lme)
  grieta.P.lme=list()
  for(i in 1:length(unique(be.pred$sample))){
      Ce=COEF.lme[i,1];me=COEF.lme[i,2];F.cte=1 ; S=1
      t=be.pred$cycles[be.pred$sample==unique(be.pred$sample)[i]]
      aoo=subset(x[,4],x[,3]==unique(x[,3])[i])[1]
      grieta.P.lme[[i]]=(aoo^(1 - me / 2) + (1 - me / 2) * Ce * F.cte^me * S^me * pi^(me / 2) * t)^(2 / (2 - me))
  }
  growt.est=c()
  for(i in 1:length(unique(be.pred$sample))) growt.est=c(growt.est,grieta.P.lme[[i]])
  Grieta.P.lme=data.frame(time=x[,1],growt.est,sample=x[,3])
  
  T=function(a_0,aF,C,m){
      ( 1- (aF/a_0)^(-(m-2)/2))/((a_0^(m/2 -1))*( m / 2 -1) * C * pi^(m / 2))
  }
  TT.M=c()
  for(i in 1:length(unique(x[,3]))){
      a000=subset(x[,4],x[,3]==unique(x[,3])[i])[1]
      T11=T(a_0=a000,aF,C=COEF.lme[i,1],m=COEF.lme[i,2])
      TT.M=c(TT.M,T11)
  }
  LL2=seq(min(x[,1]),max(TT.M),abs(max(TT.M)-min(x[,1]))/1000)
  pred.P.lme=list()
  for(i in 1:length(unique(x[,3]))){
      C=COEF.lme[i,1];m=COEF.lme[i,2];F.cte=1 ; S=1;ttt=LL2;a0=subset(x[,4],x[,3]==unique(x[,3])[i])[1]
      g.P.lme=(a0^(1 - m / 2) + (1 - m / 2) * C * F.cte^m * S^m * pi^(m / 2) * ttt)^(2 / (2 - m))
      pred.P.lme[[i]]=data.frame(ttt,g.P.lme)
  }
  growt.pred=c()
  time2=c()
  sample2=c()
  for(i in 1:length(unique(be.pred$sample))){
      growt.pred=c(growt.pred,pred.P.lme[[i]]$g.P.lme)
      time2=c(time2,LL2)
      sample2=c(sample2,rep(i,length(LL2)))
  }
  Grieta.P.lme.pred=data.frame(time=time2,growt.pred,sample=sample2)
  NN=numeric(length(unique(x[,3])))
  for(i in 1:length(unique(x[,3]))){
      NN[i]=min(which(round(pred.P.lme[[i]][,2],2)==aF))
  }
  TT.FALLO=c()
  
  for(i in 1:length(unique(x[,3]))) TT.FALLO=c(TT.FALLO,pred.P.lme[[i]][,1][NN[i]])
  
  hop <- dpik(TT.FALLO);hop
  dens.lme <- bkde(TT.FALLO, gridsize = 10000, range.x = range(TT.FALLO),bandwidth=hop)
  bin <- (range(TT.FALLO)[2] - range(TT.FALLO)[1]) / nBKDE
            F.xwj <- cumsum(abs(dens.lme$y)) * bin
            xwj <- dens.lme$x[!duplicated(F.xwj)]
            F.xwj <- unique(F.xwj) / max(F.xwj)
            
  dens.lme1=data.frame(xwj,F.xwj)
  Fe.11=c()
  for(i in 1:length(unique(x[,3]))){Fe.11[i]=(i-0.5)/length(unique(x[,3]))}
  Fe.11=round(Fe.11,3)
  ii=order(TT.FALLO);TT.FALLO.1=TT.FALLO[ii]
  dens.lme11=spline(dens.lme1$F.xwj, dens.lme1$xwj, n =101, method = "fmm",
             xmin = 0, xmax = 100, xout=seq(0, 1, length.out = 101), ties = mean)
  x.lme.1=dens.lme11$y
  Fx.lme.1=dens.lme11$x
  Dist.T.SEP_lme=data.frame(time=TT.FALLO.1,Fe=Fe.11)
  Dist.SEP_lme=data.frame(time=x.lme.1,F=Fx.lme.1)
  range.F=range(TT.FALLO)
  FFF=kde(type_kernel = "e", vec_data=TT.FALLO,y=seq(range.F[1],range.F[2], length.out =nKDE))
  xpb=FFF$grid
  F.xpb=FFF$Estimated_values
  dens.lme2=data.frame(xpb,F.xpb)
  FFF$bw
  Fe.12=c()
  for(i in 1:length(unique(x[,3]))){Fe.12[i]=(i-0.5)/length(unique(x[,3]))}
  Fe.12=round(Fe.12,3)
  jj=order(TT.FALLO);TT.FALLO.2=TT.FALLO[jj]
  dens.lme3=spline(dens.lme2$F.xpb, dens.lme2$xpb, n =101, method = "fmm",
             xmin = 0, xmax = 100, xout=seq(0, 1, length.out = 101), ties = mean)
  x.lme.2=dens.lme3$y
  Fx.lme.2=dens.lme3$x
  Dist.T.SEP_lme_kde=data.frame(time=TT.FALLO.2,Fe=Fe.12)
  Dist.SEP_lme_kde=data.frame(time=x.lme.2,F=Fx.lme.2)
  if(method[1]=="PB-nlme")
  {
      PB.NLME=list(data=x[,-4], a.F=aF, Tc=T_c, param=ORDER1.COEF, crack.est=Grieta.P.nlme,
                 sigma=sigma.r, residuals=resid1, crack.pred=Grieta.P.nlme.pred,
                 F.emp=Dist.T.PB_nlme, F.est=Dist.PB_nlme, nMC=nMC)
  }
  else if(method[1]=="SEP-lme_bkde")
  {
  
     SEP.LME=list(data=x[,-4], a.F=aF, Tc=T_c, param=COEF.lme, crack.est=Grieta.P.lme,
                sigma=ajuste.gamm$lme$sigma,  residuals=as.numeric(resid2),
                crack.pred=Grieta.P.lme.pred, F.emp=Dist.T.SEP_lme, bw=hop,
                F.est=Dist.SEP_lme, nBKDE=nBKDE)
  }
  else if(method[1]=="SEP-lme_kde")
  {
   SEP.LME=list(data=x[,-4], a.F=aF, Tc=T_c, param=COEF.lme, crack.est=Grieta.P.lme,
                sigma=ajuste.gamm$lme$sigma,  residuals=as.numeric(resid2),
                crack.pred=Grieta.P.lme.pred, F.emp=Dist.T.SEP_lme_kde,
                bw=FFF$bw, F.est=Dist.SEP_lme_kde, nKDE=nKDE)
 }
}
