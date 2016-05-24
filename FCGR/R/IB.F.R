IB.F <-
function(z, nB, alpha=0.05, method=c("SEP-lme_bkde","SEP-lme_kde","PB-nlme"))
{
  if(method[1]=="SEP-lme_bkde"){
      set.seed(1802)
      Mat.B=function(z, nB){
        x=z$data; sigma.rr=z$sigma; aF =z$a.F; param.gc=z$param; nBKDE=z$nBKDE
        bw=z$bw
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
        FF.BB=matrix(0,101,nB,byrow=F)
        for(ncol.B in 1:nB){
          f.B=function(ncol.B){
              bb=sample(1:length(unique(x[,3])), size=length(unique(x[,3])), replace = TRUE)
              param=as.matrix(param.gc); param.b=param[bb,]
              paramm1=param.b[,1]; paramm2=param.b[,2]
              grieta.BB=c()
              for(i in 1:length(unique(x[,3])) ){
                  C=paramm1[i];m=paramm2[i]
                  c=x[,1][x[,3]==unique(x[,3])[i]]
                  a_0=x[,4][x[,3]==unique(x[,3])[i]]
                  grieta.BB=c(grieta.BB,(a_0^(1 - m / 2) + (1 - m / 2) * C * pi^(m / 2) * c)^(2 / (2 - m))  + rnorm(length(c),0,sigma.rr))
              }
              DIF=data.frame(grieta.BB,x[,2],x[,3])
              g.Temp=c()
              for(i in 1:length(unique(x[,3])) ){
                  TEMP=DIF[,1][DIF[,3]==unique(DIF[,3])[i]]+
                  (x[,2][x[,3]==unique(x[,3])[i]][1]-
                  DIF[,1][DIF[,3]==unique(DIF[,3])[i]][1])
                  g.Temp=c(g.Temp,TEMP)
              }
              grieta.BB=g.Temp
              x.gamm=data.frame(cycles=x$cycles,cracks=grieta.BB,sample=as.factor(x$sample),a0=x$a0)
              ajuste.gamm <- mgcv::gamm((cracks-a0)~s(cycles,by=sample)+sample,random=list(sample=~1),data=x.gamm)
              pred.lme.B=try(predict(ajuste.gamm$lme,type="response"))
              y=data.frame(ciclos=x[,1], pred.B=as.vector(pred.lme.B),      #  y=base.escobar.pred.B
                           grieta=x[,2],muestra=x[,3])
              spline.pr.B=list()
              for(i in 1:length(unique(x[,3])) ){
                  spline.pr.B[[i]] <- pspline::smooth.Pspline(y$ciclos[y$muestra==unique(y$muestra)[i]],
                                                 y$pred.B[y$muestra==unique(y$muestra)[i]],
                                                 norder = 2, method=3)
              }
              pred.pr.B=list()
              for(i in 1:length(unique(x[,3])) ){
                  pred.pr.B[[i]] <- as.vector(predict(spline.pr.B[[i]],
                                    y$ciclos[y$muestra==unique(y$muestra)[i]],
                                    nderiv = 1))
              }
              coef.lme.B=list()
              for(i in 1:length(unique(x[,3])) ){
                  lin.mod.B <- lm(log(pred.pr.B[[i]][-1]) ~ log(sqrt(pi * x.gamm[,2][y$muestra==unique(y$muestra)[i]][-1])))
                  coef.lme.B[[i]] <- coefficients(lin.mod.B)
              }
              C.lme.B=c()
              m.lme.B=c()
              for(i in 1:length(unique(x[,3])) ){
                  C.lme.B[i]=exp(as.vector(coef.lme.B[[i]][1]))
                  m.lme.B[i]=as.vector(coef.lme.B[[i]][2])
              }
              COEF.lme.B=data.frame(C=C.lme.B,m=m.lme.B)
              COEF.lme.B
              T=function(a_0,aF,C,m){
                  ( 1- (aF/a_0)^(-(m-2)/2))/((a_0^(m/2 -1))*( m / 2 -1) * C * pi^(m / 2))
              }
              TT.M=c()
              for(i in 1:length(unique(x[,3]))){
                  a000=subset(x[,4],x[,3]==unique(x[,3])[i])[1]
                  T11=T(a_0=a000,aF,C=COEF.lme.B[i,1],m=COEF.lme.B[i,2])
                  TT.M=c(TT.M,T11)
              }
              LL2=seq(min(x[,1]),max(TT.M),abs(max(TT.M)-min(x[,1]))/1000)
              pred.P.lme.B=list()
              for(i in 1:length(unique(x[,3]))){
                  C=COEF.lme.B[i,1];m=COEF.lme.B[i,2];F.cte=1 ; S=1;ttt=LL2;a0=subset(x[,4],x[,3]==unique(x[,3])[i])[1]
                  g.P.lme=(a0^(1 - m / 2) + (1 - m / 2) * C * F.cte^m * S^m * pi^(m / 2) * ttt)^(2 / (2 - m))
                  pred.P.lme.B[[i]]=data.frame(ttt,g.P.lme)
              }
              NN.B=numeric(length(unique(x[,3])))
              for(i in 1:length(unique(x[,3]))){
                  NN.B[i]=min(which(round(pred.P.lme.B[[i]][,2],2)==aF))
              }
              TT.FALLO.B=c()
              for(i in 1:length(unique(x[,3])) ) TT.FALLO.B=c(TT.FALLO.B,pred.P.lme.B[[i]][,1][NN.B[i]])
              TT.FALLO.BB <- TT.FALLO.B[!is.na(TT.FALLO.B)]
              hop=bw
              dens.lme.B=bkde(TT.FALLO.BB, gridsize = nBKDE, range.x = range(TT.FALLO.BB),bandwidth=hop)
                    bin <- (range(TT.FALLO.BB)[2] - range(TT.FALLO.BB)[1]) / nBKDE
                    F.x.B <- cumsum(abs(dens.lme.B$y)) * bin
                    x.B <- dens.lme.B$x[!duplicated(F.x.B)]
                    F.x.B <- unique(F.x.B) / max(F.x.B)
              dens.lme1=data.frame(x.B,F.x.B)
              dens.lme11=spline(dens.lme1$F.x.B, dens.lme1$x.B, n =101, method = "fmm",
                         xmin = 0, xmax = 100, xout=seq(0, 1, length.out = 101), ties = mean)
              FF.B1=dens.lme11$y
              FF.B=round(FF.B1,6)
              FF.B
          }
          FF.BB[,ncol.B]=try(f.B(ncol.B))
        }
        return(FF.BB)
      }
      ic.b=Mat.B(z,nB)
      ICB=t(ic.b)
      quantiles<-function(x,probs=c(alpha/2,1-alpha/2))
      {quantile(x,probs)}
      ic=apply(ICB,2,quantiles)
      I.Bootstrap=data.frame(low=ic[1,],up=ic[2,])
      return(list(Mat.F.B=ic.b, I.Bootstrap=I.Bootstrap))
    }
  else if(method[1]=="SEP-lme_kde"){
    set.seed(180215)
    Mat.B=function(z, nB){
        x=z$data; sigma.rr=z$sigma; aF =z$a.F; param.gc=z$param; nKDE=z$nKDE
        baw=z$bw
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
        FF.BB=matrix(0,101,nB,byrow=F)
        for(ncol.B in 1:nB){
          f.B=function(ncol.B){
            bb=sample(1:length(unique(x[,3])), size=length(unique(x[,3])), replace = TRUE)
            param=as.matrix(param.gc); param.b=param[bb,]
            paramm1=param.b[,1]; paramm2=param.b[,2]
            grieta.BB=c()
            for(i in 1:length(unique(x[,3])) ){
                C=paramm1[i];m=paramm2[i]
                c=x[,1][x[,3]==unique(x[,3])[i]]
                a_0=x[,4][x[,3]==unique(x[,3])[i]]
                grieta.BB=c(grieta.BB,(a_0^(1 - m / 2) + (1 - m / 2) * C * pi^(m / 2) * c)^(2 / (2 - m))  + rnorm(length(c),0,sigma.rr))
            }
            DIF=data.frame(grieta.BB,x[,2],x[,3])
            g.Temp=c()
            for(i in 1:length(unique(x[,3])) ){
                TEMP=DIF[,1][DIF[,3]==unique(DIF[,3])[i]]+
                (x[,2][x[,3]==unique(x[,3])[i]][1]-
                DIF[,1][DIF[,3]==unique(DIF[,3])[i]][1])
                g.Temp=c(g.Temp,TEMP)
            }
            grieta.BB=g.Temp
            x.gamm=data.frame(cycles=x$cycles,cracks=grieta.BB,sample=as.factor(x$sample),a0=x$a0)
            ajuste.gamm <- mgcv::gamm((cracks-a0)~s(cycles,by=sample)+sample,random=list(sample=~1),data=x.gamm)
            pred.lme.B=try(predict(ajuste.gamm$lme,type="response"))
            y=data.frame(ciclos=x[,1], pred.B=as.vector(pred.lme.B),      #  y=base.escobar.pred.B
                     grieta=x[,2],muestra=x[,3])
            spline.pr.B=list()
            for(i in 1:length(unique(x[,3])) ){
                spline.pr.B[[i]] <- pspline::smooth.Pspline(y$ciclos[y$muestra==unique(y$muestra)[i]],
                                               y$pred.B[y$muestra==unique(y$muestra)[i]],
                                               norder = 2, method=3)
            }
            pred.pr.B=list()
            for(i in 1:length(unique(x[,3])) ){
                pred.pr.B[[i]] <- as.vector(predict(spline.pr.B[[i]],
                            y$ciclos[y$muestra==unique(y$muestra)[i]],
                            nderiv = 1))
            }
            coef.lme.B=list()
            for(i in 1:length(unique(x[,3])) ){
                lin.mod.B <- lm(log(pred.pr.B[[i]][-1]) ~ log(sqrt(pi * x.gamm[,2][y$muestra==unique(y$muestra)[i]][-1])))
                coef.lme.B[[i]] <- coefficients(lin.mod.B)
            }
            C.lme.B=c()
            m.lme.B=c()
            for(i in 1:length(unique(x[,3])) ){
                C.lme.B[i]=exp(as.vector(coef.lme.B[[i]][1]))
                m.lme.B[i]=as.vector(coef.lme.B[[i]][2])
            }
            COEF.lme.B=data.frame(C=C.lme.B,m=m.lme.B)
            COEF.lme.B
            T=function(a_0,aF,C,m){
              ( 1- (aF/a_0)^(-(m-2)/2))/((a_0^(m/2 -1))*( m / 2 -1) * C * pi^(m / 2))
            }
            TT.M=c()
            for(i in 1:length(unique(x[,3]))){
                a000=subset(x[,4],x[,3]==unique(x[,3])[i])[1]
                T11=T(a_0=a000,aF,C=COEF.lme.B[i,1],m=COEF.lme.B[i,2])
                TT.M=c(TT.M,T11)
            }
            LL2=seq(min(x[,1]),max(TT.M),abs(max(TT.M)-min(x[,1]))/1000)
            pred.P.lme.B=list()
            for(i in 1:length(unique(x[,3]))){
                C=COEF.lme.B[i,1];m=COEF.lme.B[i,2];F.cte=1 ; S=1;ttt=LL2;a0=subset(x[,4],x[,3]==unique(x[,3])[i])[1]
                g.P.lme=(a0^(1 - m / 2) + (1 - m / 2) * C * F.cte^m * S^m * pi^(m / 2) * ttt)^(2 / (2 - m))
                pred.P.lme.B[[i]]=data.frame(ttt,g.P.lme)
            }
            NN.B=numeric(length(unique(x[,3])))
            for(i in 1:length(unique(x[,3]))){
                NN.B[i]=min(which(round(pred.P.lme.B[[i]][,2],2)==aF))
            }
            TT.FALLO.B=c()
            for(i in 1:length(unique(x[,3])) ) TT.FALLO.B=c(TT.FALLO.B,pred.P.lme.B[[i]][,1][NN.B[i]])
            TT.FALLO.BB <- TT.FALLO.B[!is.na(TT.FALLO.B)]
            range.F=range(TT.FALLO.BB)
            hop=baw
            FFF=kde(type_kernel = "e", vec_data=TT.FALLO.BB,y=seq(range.F[1],range.F[2],
                    length.out =nKDE),bw=hop)
            xpb1=FFF$grid; ii=order(xpb1);xpb=FFF$grid[ii]
            F.xpb=FFF$Estimated_values[ii]
            dens.lme2=data.frame(xpb,F.xpb)
            dens.lme3=spline(dens.lme2$F.xpb, dens.lme2$xpb, n =101, method = "fmm",
                      xmin = 0, xmax = 100, xout=seq(0, 1, length.out = 101), ties = mean)
            FF.B1=dens.lme3$y
            FF.B2=round(FF.B1,6)
            iii=order(FF.B2)
            FF.B=FF.B2[iii]
        }
        FF.BB[,ncol.B]=try(f.B(ncol.B))
      }
      return(FF.BB)
    }
    ic.b=Mat.B(z,nB)
    ICB=t(ic.b)
    quantiles<-function(x,probs=c(alpha/2,1-alpha/2))
    {quantile(x,probs)}
    ic=apply(ICB,2,quantiles)
    I.Bootstrap=data.frame(low=ic[1,],up=ic[2,])
    return(list(Mat.F.B=ic.b, I.Bootstrap=I.Bootstrap))
  }
  else if(method[1]=="PB-nlme"){
      set.seed(18)
      Mat.B=try(function(z, nB){
        x=z$data; sigma.r=z$sigma; aF =z$a.F; param.gc=z$param; nMC=z$nMC
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
        Mat.F.B=matrix(0,101,nB,byrow=F)
        k=1
        while(min(Mat.F.B[1,])==0){
            j=1
            while(j==1){
              bb=sample(1:length(unique(x[,3])), size=length(unique(x[,3])), replace = TRUE)
              param=as.matrix(param.gc); param.b=param[bb,]
              param1=param.b[,1]; param2=param.b[,2]
              grieta.B=c()
              for(i in 1:length(unique(x[,3])) ){
                  Cc=param1[i];mM=param2[i]; cyc=x[,1][x[,3]==unique(x[,3])[i]]
                  a_0=x[,4][x[,3]==unique(x[,3])[i]]
                  grieta.B=c(grieta.B,(a_0^(1 - mM / 2) + (1 - mM / 2) * Cc * pi^(mM / 2) * cyc)^(2 / (2 - mM))+ rnorm(length(cyc),0,sigma.r))
              }
              DIFF=data.frame(grieta.B,x[,3])
              d=0
              for(i in 1:length(unique(x[,3])) ){
              if(min(diff(DIFF[,1][DIFF[,2]==unique(DIFF[,2])[i]]))>=0) d=d+1
              }
              if(d==length(unique(x[,3])))j=2
            }
        DIF=data.frame(grieta.B,x[,2],x[,3])
        g.Temp=c()
        for(i in 1:length(unique(x[,3])) ){
            TEMP=DIF[,1][DIF[,3]==unique(DIF[,3])[i]]+
            (x[,2][x[,3]==unique(x[,3])[i]][1]-
            DIF[,1][DIF[,3]==unique(DIF[,3])[i]][1])
            g.Temp=c(g.Temp,TEMP)
        }
        grieta.B=g.Temp
        x.gD.BB=data.frame(cycles=x$cycles,cracks=grieta.B,sample=x$sample,a0=x$a0)
        Paris.F <- function(cycles, C, m, a0, F.cte=1, S = 1)
                 (a0^(1 - m / 2) + (1 - m / 2) * C * F.cte^m * S^m * pi^(m / 2) * cycles)^(2 / (2 - m))
        DerivInit <- function(data, cycles, y){
                xx <- data$cycles
                y <- data$cracks
                derivada <- diff(y)/diff(xx)
                grieta.media <- (y[-1]+y[-length(y)])/2
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
        
        environment(SSParis) <- .GlobalEnv
        assign("SSParis", SSParis, envir = environment(SSParis))       
        x.gD.BB <<- x.gD.BB
        
        
        model<-NULL
        try(model<-try(nlsList(cracks ~ SSParis(cycles, C, m, a0,
                               F.cte =1, S=1)| sample, na.action = na.omit,
                               data=x.gD.BB),silent=TRUE),silent=TRUE)
        if(!is.null(model)==TRUE){
               MB1=try(model,silent=TRUE)
               NN=try(min(coef(MB1)[,2]), silent=TRUE)
               MM=try(min(coef(MB1)[,1]),silent=TRUE)
            if( NN>0 & MM>0 & (!is.na(NN)==TRUE) & (!is.na(MM)==TRUE)){
                model.nlme.B<-try(nlme(MB1,random = C + m ~ 1,na.action = na.omit),silent=TRUE)
                if(class(model.nlme.B)[1]!="try-error"){
                    C.B=mean(coef(model.nlme.B)[,1])
                    m.B=mean(coef(model.nlme.B)[,2])
                    COV.B=var(coef(model.nlme.B))
                    T.B=function(a_0,aF,C,m){
                                ( 1- (aF/a_0)^(-(m-2)/2))/((a_0^(m/2 -1))*( m / 2 -1) * C * pi^(m / 2))
                    }
                    t_fallo.B=numeric(nMC)
                    for(i in 1:length(t_fallo.B)){
                        param= MASS::mvrnorm(mu = c(C.B,m.B), Sigma = COV.B)
                        t_fallo.B[i]=T.B(a_0=min(x[,4]),aF,C=param[1],m=param[2] )
                    }
                    BOXP.B=boxplot.stats(t_fallo.B,coef=3)$stats
                    t.fallo.2.B=t_fallo.B[0<t_fallo.B & t_fallo.B<=max(BOXP.B)]
                    Mat.F.B[,k]=as.vector(quantile(t.fallo.2.B, probs = seq(0, 1, 0.01)))
                    k=k+1
                }
            }
          }
        }
          return(Mat.F.B)
      }, silent=TRUE)
      ic.b=try(Mat.B(z,nB),silent=TRUE)
      ICB.nlme=t(ic.b)
      quantiles<-function(x,probs=c(alpha/2,1-alpha/2))
      {quantile(x,probs)}
      ic.nlme=apply(ICB.nlme,2,quantiles)
      I.Bootstrap=data.frame(low=ic.nlme[1,],up=ic.nlme[2,])
      return(list(Mat.F.B=ic.b, I.Bootstrap=I.Bootstrap))
    }
}
