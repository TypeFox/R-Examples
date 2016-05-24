Survgini <-
function(asymptotic=TRUE, permutation=TRUE, M=500, linearRank=TRUE, lastEvent=1, datasetOrig){
    library(survival)
    VarGinicensor<-function(data, Tmax){
        data1<-data.frame(data)
        info<- survfit(Surv(data1[,1], data1[,2])~1 , type='kaplan-meier', data=data1)
        S<-matrix(info$surv)
        T<-matrix(info$time)
        indices <- sum(1*(T<=Tmax))
        if (indices==0){
            var<-0
            return(var=var)
        }
        else{
            Vt<-matrix(ncol=1,nrow=indices)
            Vt[1]<-1*T[1]
            if(indices>=2){
                for (i in 2:indices){
                    Vt[i]<-Vt[i-1]+((S[i-1])^2)*(T[i]-T[i-1])
                }
            }
            lastpiecesVt<-((S[indices])^2)*(Tmax-T[indices])
            VtMax<-Vt[indices]+lastpiecesVt
            Wt<-matrix(ncol=1,nrow=indices)
            Wt[1]<-1*T[1]
            if(indices>=2){
                for (i in 2:indices){
                    Wt[i]<-Wt[i-1]+S[i-1]*(T[i]-T[i-1])
                }
            }
            lastpiecesWt<-(S[indices])*(Tmax-T[indices])
            WtMax<-Wt[indices]+lastpiecesWt
            mu2<-matrix(ncol=1,nrow=indices)
            mu2<-VtMax-Vt
            mu<-matrix(ncol=1,nrow=indices)
            mu<-WtMax-Wt
            n<-length(data[,1])
            event<-matrix( ncol=1, nrow=indices)
            atrisk<-matrix(ncol=1, nrow=indices)
            for (i in 1:indices){
                event[i]<-info$n.event[i]
                atrisk[i]<-info$n.risk[i]
            }
            dsigma<-matrix(0, ncol=1,nrow=indices)
            for (i in 1:indices){
                if (atrisk[i]>0)
                dsigma[i]<-(n*event[i])/(atrisk[i]^2)
            }
            varistant<-matrix(ncol=1, nrow=indices)
            for (i in 1:indices){
                varistant[i]<-(4*exp(2*log(mu2[i])-2*log(WtMax))+exp(2*log(mu[i])+
                2*log(VtMax)-4*log(WtMax))-4*exp(log(mu[i])+log(mu2[i])+log(VtMax)-
                3*log(WtMax)))*(dsigma[i])
            }
            var<-matrix(ncol=1,nrow=indices)
            var[1]<-varistant[1]
            if(indices>=2){
                for (i in 2:indices){
                   var[i]<-var[i-1]+varistant[i]
                }
            }
            return(var[indices]/n)
        }
    }
    Gcensor2<-function(data,Tmax=max(data[,1])){
        data1 <- data.frame(data)
        info <- survfit(Surv(data[,1], data[,2])~1 , type='kaplan-meier', data=data1)
        K<-length(info$time)
        S<-matrix(info$surv, ncol=1, nrow=K)
        num<-matrix(ncol=1,nrow=K)
        T<-matrix(info$time,ncol=1, nrow=K)
        num[1]<-1*T[1]
        if (K>=2){
            for (i in 2:K)
            num[i]<-num[i-1]+((S[i-1])^2)*(T[i]-T[i-1])
        }
        den <- matrix(ncol=1,nrow=K)
        den[1] <- 1*T[1]
        if (K>=2){
            for (i in 2:K)
            den[i]<-den[i-1]+S[i-1]*(T[i]-T[i-1])
        }
        G <- 1-(num/den)
        indices <- sum(1*(T<Tmax))
        lastpiecesnum <- ((S[indices])^2)*(Tmax-T[indices])
        lastpiecesden <- (S[indices])*(Tmax-T[indices])
        GTmax <- 1-(num[indices]+lastpiecesnum)/(den[indices]+lastpiecesden)
        return(list( G=G,info=info,GTmax=GTmax))
    }
    X <- datasetOrig[,1]
    delta <- datasetOrig[,2]
    Tx <- datasetOrig[,3]
    dataset1Orig <- datasetOrig[Tx==1,]
    dataset2Orig <- datasetOrig[Tx==2,]
    N<-length(datasetOrig[,1])
    N1<-length(dataset1Orig[,1])
    N2<-length(dataset2Orig[,1])
    if (lastEvent==1) A1<-dataset1Orig[,1]*1*(dataset1Orig[,2]==1) else A1<-dataset1Orig[,1]
    if (lastEvent==1) A2<-dataset2Orig[,1]*1*(dataset2Orig[,2]==1) else A2<-dataset2Orig[,1]
    Tmax1<-max(A1)+0.001
    Tmax2<-max(A2)+0.001
    Tmaxnum<-Tmax1*1*(N1>=N2)+Tmax2*1*(N2>N1)
    Ginisurv1 <- Gcensor2(dataset1Orig[,-3], Tmax=Tmaxnum)
    Ginisurv2 <- Gcensor2(dataset2Orig[,-3], Tmax=Tmaxnum)
    Ginis1Orig<- Ginisurv1$GTmax
    Ginis2Orig<- Ginisurv2$GTmax
    teststatGiniOrig<-(Ginis1Orig - Ginis2Orig)^2
    if (asymptotic==TRUE){
        VarasintGini1 <- VarGinicensor(dataset1Orig[,-3], Tmax1)
        VarasintGini2 <- VarGinicensor(dataset2Orig[,-3], Tmax2)
        teststatGiniAs<-(Ginis1Orig - Ginis2Orig)/sqrt(VarasintGini1 + VarasintGini2)
        pGiniAs<-1-pchisq((teststatGiniAs^2), df=1)
    }
    teststatGiniPerm<-rep(NA, M)
    if (permutation==TRUE){
        for(msims in 1:M){
            label<-sample(Tx)
            dataset1<-matrix(ncol=2, nrow=N1)
            dataset2<-matrix(ncol=2, nrow=N2)
            dataset1[,1]<-datasetOrig[label==1,1]
            dataset1[,2]<-datasetOrig[label==1,2]
            dataset2[,1]<-datasetOrig[label==2,1]
            dataset2[,2]<-datasetOrig[label==2,2]
            if (lastEvent==1) A1Perm<-dataset1[,1]*1*(dataset1[,2]==1) else A1Perm<-dataset1[,1]
            if (lastEvent==1) A2Perm<-dataset2[,1]*1*(dataset2[,2]==1) else A2Perm<-dataset2[,1]
            Tmax1Perm<-max(A1Perm)+0.001
            Tmax2Perm<-max(A2Perm)+0.001
            TmaxnumPerm<-Tmax1Perm*1*(N1>=N2)+Tmax2Perm*1*(N2>N1)
            Ginisurv1Perm <- Gcensor2(dataset1, Tmax=TmaxnumPerm)
            Ginisurv2Perm <- Gcensor2(dataset2, Tmax=TmaxnumPerm)
            Ginis1Perm<- Ginisurv1Perm$GTmax
            Ginis2Perm<- Ginisurv2Perm$GTmax
            teststatGiniPerm[msims]<-(Ginis1Perm-Ginis2Perm)^2
        }
        pGiniPerm<- sum( 1*(teststatGiniPerm > teststatGiniOrig) )/M
    }
    if (linearRank==TRUE){
        datatemp <- list(X=X,delta=delta, Tx=Tx)
        teststatGT <- survdiff(Surv(X,delta)~Tx,datatemp,rho=-1)$chisq
        teststatLR <- survdiff(Surv(X,delta)~Tx,datatemp,rho=0)$chisq
        teststatW <- survdiff(Surv(X,delta)~Tx,datatemp,rho=1)$chisq
        pLR<-1-pchisq(teststatLR, df=1)
        pGT<-1-pchisq(teststatGT, df=1)
        pW<-1-pchisq(teststatW, df=1)
    }
    if (permutation==TRUE & asymptotic==TRUE & linearRank==TRUE) 
	print(c(pGiniAs=pGiniAs, pGiniPerm=pGiniPerm, pGT=pGT,pLR=pLR,pW=pW), digits=5)
    else{
        if (permutation==FALSE & asymptotic==TRUE & linearRank==TRUE) 
            print(c(pGiniAs=pGiniAs, pGT=pGT,pLR=pLR,pW=pW), digits=5)
        if (permutation==TRUE & asymptotic==FALSE & linearRank==TRUE) 
            print(c(pGiniPerm=pGiniPerm, pGT=pGT,pLR=pLR,pW=pW), digits=5)
        if (permutation==TRUE & asymptotic==TRUE & linearRank==FALSE) 
            print(c(pGiniAs=pGiniAs, pGiniPerm=pGiniPerm), digits=5)
        if (permutation==FALSE & asymptotic==FALSE & linearRank==TRUE) 
            print(c(pGT=pGT,pLR=pLR,pW=pW), digits=5)
        if (permutation==FALSE & asymptotic==TRUE & linearRank==FALSE) 
            print(c(pGiniAs=pGiniAs), digits=5)
        if (permutation==TRUE & asymptotic==FALSE & linearRank==FALSE) 
            print(c(pGiniPerm=pGiniPerm ), digits=5)
        if (permutation==FALSE & asymptotic==FALSE & linearRank==FALSE) 
            print("No test chosen", digits=5)
    }
  }

