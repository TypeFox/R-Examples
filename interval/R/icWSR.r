wsrMC<-function(A,p,x,group,alternative,nwsr,np,digits=12){
    if (is.numeric(group)){ ng<-0
    } else { 
        ug<-unique(group)
        ng<-length(ug)
    }
    ## sample within subject in proportion to conditional probability density for that subject given subjects interval
    wsrResample<-function(Arow){
        sample(x,1,replace=TRUE,prob=Arow*p/sum(Arow*p) )
    }
    n<-dim(A)[1]
    ## use a different test Statatistic (testStat function) for different type of group objects 
    if (ng==2){
        ## two-sample tests
        n1<-length(group[group==ug[1]])
        n2<-n-n1
        z<- rep(NA,n)
        z[group==ug[1]]<- 1/n1
        z[group==ug[2]]<- -1/n2    
        testStat<-function(x){ sum(x*z) }
    } else if (ng==0){
        ## correlation tests
        testStat<-function(x){ sum(x*group) }
    } else {
        ## k-sample tests
        testStat<-function(x){
            out<-0
            x<- x- mean(x)
            for (i in 1:ng){
                I<- group==ug[i]
                out<-out+ length(x[I])*(mean(x[I]))^2
            }
            out
        }
    }
    ## Tij is a matrix of test statistics, 
    ## Tij[i,j] is the ith within subject resample and the jth permutation
    Tij<-matrix(NA,nwsr,np+1)
    for (i in 1:nwsr){  
        Cj<-apply(A,1,wsrResample)
        Tij[i,1]<- testStat(Cj)
        for (j in 2:(np+1)){
            Tij[i,j]<-testStat(Cj[sample(1:n,n,replace=FALSE)]) 
        }
    }
    ## round to avoid problems with ties
    Tij<-signif(Tij,digits=digits)
    if (ng<=2){
        ## one-sided tests make sense for ng<=2
        cnt.gte<-function(x){ length(x[x>=x[1]]) -1 }
        total.gte<-sum(apply(Tij,1,cnt.gte))
        p.gte<- (total.gte+1)/(nwsr*np+1)
        cnt.lte<-function(x){ length(x[x<=x[1]]) -1 }
        total.lte<-sum(apply(Tij,1,cnt.lte))
        p.lte<- (total.lte+1)/(nwsr*np+1)
        cnt.abs<-function(x){ length(x[abs(x)>=abs(x[1])]) -1 }
        total.abs<-sum(apply(Tij,1,cnt.abs))
        p.twosidedAbs<- (total.abs+1)/(nwsr*np+1)
        p.twosided<-min(1,2*min(p.lte,p.gte))
        p.values<-c(p.twosided=p.twosided,p.lte=p.lte,p.gte=p.gte,p.twosidedAbs=p.twosidedAbs)
        p.value<-switch(alternative,
            less=p.lte,
            greater=p.gte,
            two.sided=p.twosided,
            two.sidedAbs=p.twosidedAbs)
    } else {
        ## for k-sample tests one-sided tests do not make sense
        cnt.gte<-function(x){ length(x[x>=x[1]]) -1 }
        total.gte<-sum(apply(Tij,1,cnt.gte))
        p.twosided<- (total.gte+1)/(nwsr*np+1)
        p.values<-c(p.twosided=p.twosided,p.twosidedAbs=p.twosided)
        p.value<-p.twosided
    }
    out<-list(p.value=p.value,p.values=p.values,statistic=NULL,parameter=NULL)
    out
}

wsrHLYpclt<-function(A,p,x,group,alternative,nwsr,type,RHO){
    n<-dim(A)[1]
    ng<-length(unique(group))
    ug<-unique(group)
    if (type=="HLY"){
        N<-length(p)
        wsrResample<-function(Arow){
            sample(1:N,1,replace=TRUE,prob=Arow*p/sum(Arow*p) )
        }
    } else {
        wsrResample<-function(Arow){
            sample(x,1,replace=TRUE,prob=Arow*p/sum(Arow*p) )
        }
    }
    ## Ui is the efficient score 
    ## Vi is the associated variance
    ## getUiVi.HLY uses usual martingale methods from survdiff
    getUiVi.HLY<-function(x,group,ng,RHO){
        lrtmp<-survdiff(Surv(x,rep(1,n))~group,rho=RHO)     
        Ui<-lrtmp$obs[1:(ng-1)] - lrtmp$exp[1:(ng-1)]
        Vi<-lrtmp$var[1:(ng-1),1:(ng-1)]
        list(Ui=Ui,Vi=Vi)
    }
    ## getUiVi.pclt2 and getUiVi.pclt3 use permutational central limit theorem
    ##     pclt2 for two-sample or less
    ##     pclt3 for three-sample or more
    getUiVi.pclt2<-function(x,z,zbar,n){
        xbar<-mean(x)
        Ui<-sum(x*z) - n*xbar*zbar
        Vi<-(1/(n-1)) * sum( (x-xbar)^2 ) * sum( (z-zbar)^2 )
        list(Ui=Ui,Vi=Vi)
    }
    getUiVi.pclt3<-function(x,zMat,zMatbar,n,zSSE){
        xbar<-mean(x)
        Ui<-apply(zMat,2,function(zi){ sum(x*zi)}) - n*xbar*zMatbar
        Vi<-(1/(n-1)) * sum( (x-xbar)^2 ) * zSSE
        list(Ui=Ui,Vi=Vi)
    }


    if (ng>2){
        Ui<-matrix(NA,nwsr,ng-1)
        Vi<-array(NA,c(nwsr,ng-1,ng-1))
        if (type=="pclt"){    
            zMat<-model.matrix(~-1+factor(group),contr="contr.treatment")[,1:(ng-1)]
            zMatbar<-apply(zMat,2,mean)
            zSSE<-matrix(0,ng-1,ng-1)
            for (i in 1:n){
                zSSE<-zSSE + matrix(zMat[i,]-zMatbar,ng-1,1) %*% 
                    matrix(zMat[i,]-zMatbar,1,ng-1)
            }
        }
        for (i in 1:nwsr){  
            Cj<-apply(A,1,wsrResample)
            if (all(Cj==Cj[1])){
                Ui[i,]<-0
                Vi[i,,]<-0
            } else{
                if (type=="HLY"){
                   tmp<-getUiVi.HLY(Cj,group,ng,RHO)
                   Ui[i,]<-tmp$Ui
                   Vi[i,,]<-tmp$Vi
                } else {
                   tmp<-getUiVi.pclt3(Cj,zMat,zMatbar,n,zSSE)
                   Ui[i,]<-tmp$Ui
                   Vi[i,,]<-tmp$Vi
                }
            }
        }
        Ubar<-apply(Ui,2,mean)
        Vhat<-matrix(0,ng-1,ng-1)
        for (i in 1:nwsr){
            Vhat<-Vhat+ (1/nwsr)*Vi[i,,] - 
               (1/(nwsr-1))* matrix(Ui[i,]-Ubar,ng-1,1) %*% matrix(Ui[i,]-Ubar,1,ng-1)
        }
        chisq.value <- matrix(Ubar,1,ng-1) %*% solve(Vhat) %*% matrix(Ubar,ng-1,1)
        df<-ng-1
        p.twosided <- 1 - pchisq(chisq.value, df)
        p.values<-c(p.twosided=p.twosided,p.twosidedAbs=p.twosided)
        p.value<-p.twosided
        statistic<-chisq.value
        names(statistic)<-"Chi Square"
        parameter<-df
        names(parameter)<-"df" 
    } else {
        Ui<-Vi<-rep(NA,nwsr)
        if (type=="pclt"){    
            n1<-length(group[group==ug[1]])
            n2<-n-n1
            z<- rep(NA,n)
            z[group==ug[1]]<- 1/n1
            z[group==ug[2]]<- -1/n2 
            zbar<-mean(z)
        }
        for (i in 1:nwsr){  
            Cj<-apply(A,1,wsrResample)
            if (all(Cj==Cj[1])){
                Ui[i]<-0
                Vi[i]<-0
            } else{
                if (type=="HLY"){
                   tmp<-getUiVi.HLY(Cj,group,ng,RHO)
                   Ui[i]<-tmp$Ui
                   Vi[i]<-tmp$Vi
                } else {
                   tmp<-getUiVi.pclt2(Cj,z,zbar,n)
                   Ui[i]<-tmp$Ui
                   Vi[i]<-tmp$Vi
                }
            }
        }
        Ubar<-mean(Ui)
        Vhat<- sum(Vi)/nwsr - sum( (Ui-Ubar)^2 )/(nwsr-1)
        Z<-Ubar/sqrt(Vhat)
        p.lte<-pnorm(Z)
        p.gte<-1-pnorm(Z)
        p.twosidedAbs<- 1-pchisq(Z^2,1)
        # Note for normal theory p-values, p.twosided=p.twosidedAbs
        p.values<-c(p.twosided=p.twosidedAbs,p.lte=p.lte,p.gte=p.gte,p.twosidedAbs=p.twosidedAbs)
        if (alternative=="less" | alternative=="greater"){
            statistic<-Z
            names(statistic)<-"Z"
            parameter<-NULL
        } else {
            statistic<-Z^2
            names(statistic)<-"Chi Square"
            parameter<-1
            names(parameter)<-"df"
        }
        p.value<-switch(alternative,
            less=p.lte,
            greater=p.gte,
            two.sided=p.twosidedAbs,
            two.sidedAbs=p.twosidedAbs)
    }
    out<-list(p.value=p.value,p.values=p.values,statistic=statistic,parameter=parameter)
    out
}
## main function for Within Subject Resampling
icWSR<-function(fit,group,scores,alternative,type, control){
    nwsr<-control$nwsr
    np<-control$np
    digits<-control$digits
    seed<-control$seed
    setSEED<-control$setSEED

    if (setSEED) set.seed(seed)

    n<-dim(fit$A)[1]
    m<-length(fit$pf)
    if (is.factor(group) | length(unique(group))==2) group<-as.character(group)
    if (is.character(group)){
        ug<-unique(group)
        ng<-length(ug)
    } else if (is.numeric(group)){
        ng<-0
    } else stop("group must be either a factor, character or numeric vector")
    calc.scorej<-function(p,scores){
        Sj<-1-cumsum(p)
        Sj[Sj<0]<-0
        SL<-c(1,Sj[-m])
        SR<-Sj
        if (scores == "logrank1"){
            Lambdaj<-  - cumsum(p/SL) 
            LL<-c(0,Lambdaj[-m])
            LR<-Lambdaj
            sout<-(SL*LL - SR*LR)/(SL-SR)
        } else if (scores=="logrank2"){
            logSR<-log(SR)
            logSR[logSR==-Inf]<-0
            sout<-(SL*log(SL) - SR*logSR)/(SL-SR) 
        } else if (scores=="wmw"){
            sout<-SL+SR-1
        } else stop("scores must be 'logrank1', 'logrank2', or 'wmw' ")
        sout
    }
    scorej<-calc.scorej(fit$pf,scores)
    if (type=="wsr.mc"){
        wout<-wsrMC(fit$A,fit$pf,scorej,group,alternative,nwsr,np,digits)
    } else if (type=="wsr.pclt"){
        wout<-wsrHLYpclt(fit$A,fit$pf,scorej,group,alternative,nwsr,"pclt",RHO=NA)   
    } else if (type=="wsr.HLY"){
        if (ng==0) stop("group variable interpreted as nummeric, no trend test for method 'wsr.HYL'")
        if (scores=="logrank1"){ RHO<-0
        } else if (scores=="wmw"){ RHO<-1
        } else stop("only scores 'logrank1' or rho=0 and 'wmw' or rho=1 are supported for method='wsr.HLY'")
        wout<-wsrHLYpclt(fit$A,fit$pf,scorej,group,alternative,nwsr,"HLY",RHO)
    }
    wout
}
