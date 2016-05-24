schapt<-function(x, n.boots = 0, replace = FALSE, alternative = c("one.change", "epidemic"), conf.level = 0.95, 
        adj.Wn = FALSE, model.test = FALSE, n.model.boots = 0, tol=1.0e-7, maxit=50,trace=FALSE,... )
{
    n.boots<-as.integer(n.boots)
    n.model.boots<-as.integer(n.model.boots)
    if(n.boots>10 || n.model.boots>10) cat("It will take a while to do the bootstrap test(s).\n Please be patient.\n")
    DNAME <- deparse(substitute(x))    
    alternative <- match.arg(alternative)
    
    if(!is.matrix(x)) x<-as.matrix(x)
    if(is.matrix(x)){
        N<-nrow(x)
        r<-ncol(x)
    }
    if(is.vector(x)){
        N<-length(x); r<-1
        x<-matrix(x,N,1)
    }
    df <- r; nc<-0; 
    ifelse(alternative == "epidemic", nc<-2, nc<-1)
    res<-SemiparChangePoint(x, alternative = alternative, adj.Wn = adj.Wn, tol=tol, maxit=maxit,trace=trace)
    k.hat <-res$k.hat; alpha.hat <- res$alpha.hat;  beta.hat <- res$beta.hat
    if(nc == 1){
        Boots.Sn<-NULL
        Sn<-res$Sn
        pvalue<-p.OneChange(n=N,d=r,Sn)
        m.hat<-N
        CV<-Sn.alfa(1-conf.level,N,r)
#        cat("CV=",CV,"\n")
#        CI<-(1:(N+1))[2*(res$ll-res$ll[1]) > CV]
#        cat("CI=",CI,"\n")
#        if(length(CI)==0) cint<-c(NA,NA)
#        else cint<-c(min(CI)-1, max(CI)-1)
        #### this is not correct.
    }
    if(nc == 2){
        Vn<-res$Vn; Wn<-res$Wn
        pvalue.Vn<-p.Epidemic.Vn(Vn,r, tol=0.000001)
        pvalue.Wn<-p.Epidemic.Wn(Wn, tol=0.000001)
        m.hat<-res$m.hat
    }
    if(n.boots < 0 ) stop("incorrect number of bootstrap (permutation) samples")
    if(n.boots >0 && n.boots < 10 ) warning("too small number of bootstrap (permutation) samples")
    if(n.boots >= 1){
        if(nc == 1) Freq<-0
        else{ Freq.Vn<-0; Freq.Wn<-0}
        for(b in 1:n.boots){
            y<-x[sample(1:N, N, replace = replace ),]
            Temp<-SemiparChangePoint(y, alternative = alternative, adj.Wn = adj.Wn, tol = tol, maxit = maxit,trace = F)
            if(nc == 1){
#                Boots.Sn[b]<-Temp$Sn
                Freq<-Freq+(Temp$Sn>=Sn)
            }
            else{
                Freq.Vn<-Freq.Vn+(Temp$Vn>=Vn); 
                Freq.Wn<-Freq.Wn+(Temp$Wn>=Wn)}
        }
        if(nc==1){
#            p.boots<-mean(Boots.Sn>=Sn)
#            cat("qtl=",quantile(Boots.Sn, conf.level),"\n")
#            bCI<-(1:(N+1))[2*(res$ll-res$ll[1]) > quantile(Boots.Sn, conf.level) ]
#            if(length(bCI)==0) boot.cint<-c(NA,NA)
#            else boot.cint<-c(min(bCI)-1, max(bCI)-1)
            #### this is not correct.
            p.boots<-Freq/n.boots
        }
        else{ p.boots.Vn<-Freq.Vn/n.boots; p.boots.Wn<-Freq.Wn/n.boots }
    }
    if(model.test && n.model.boots > 0){
        model.res<-BootsModelTest(x,k.hat,m.hat,B=n.model.boots,alpha.hat, beta.hat, tol=tol, maxit=maxit,trace=F)
        p.value.model<-model.res$Pvalue
#        cat("p=",p.value.model,"\n")
        Delta = model.res$Delta
    }
    if(nc ==1){ para<-c(N,df); names(para)<-c("sample size", "df")
        chpt<-k.hat
        names(chpt)<-"k"
        est<-c(alpha.hat,  beta.hat)
        names(est)<-c("alpha","beta")
        if(!model.test || n.model.boots == 0 ){ 
            stats <-Sn; names(stats)<-c("Sn")
            if(n.boots >= 1)  p.value<-c(pvalue,p.boots)
            else p.value<-pvalue
        }
        else{
            stats <-c(Sn,Delta); names(stats)<-c("Sn", "Delta")
            if(n.boots >= 1 ) p.value<-c(pvalue,p.boots,p.value.model)
            else p.value<-c(pvalue,p.value.model)
        }
    }
    if(nc ==2){ para<-N; names(para)<-"sample size"
        chpt <-c(k.hat, m.hat)
        names(chpt)<-c("k","m")
        est <-c(alpha.hat, beta.hat)
        names(est)<-c("alpha","beta")
        if(!model.test || n.model.boots == 0 ){
            stats <-c(Vn, Wn); names(stats)<-c("Vn", "Wn")
            if(n.boots >= 1 ) p.value<-c(pvalue.Vn,pvalue.Wn,p.boots.Vn, p.boots.Wn)
            else p.value<-c(pvalue.Vn,pvalue.Wn)
        }
        else{
            stats <-c(Vn, Wn, Delta); names(stats)<-c("Vn", "Wn", "Delta")
            if(n.boots >= 1 ) p.value<-c(pvalue.Vn,pvalue.Wn,p.boots.Vn, p.boots.Wn,p.value.model)
            else p.value<-c(pvalue.Vn,pvalue.Wn,p.value.model)
        }
    }
    if(nc == 2){
        RVAL <- list(
                    data.name = DNAME,
                    parameter = para,
                    alternative = alternative,
                    statistic = stats,
                    change.points = chpt,
                    estimate = est,
                    p.value = p.value
                    )
    }
    if(nc == 1){# && length(CI)>0){
#        names(cint)<-c(paste(as.character(conf.level), "\% confidence interval: chpt.l",sep="\0"), "chpt.u")    
        RVAL <- list(
                 data.name = DNAME,
                 parameter = para,
                 alternative = alternative,
                 statistic = stats,
                 change.point = chpt,
                 estimate = est,
                 p.value = p.value)#, 
#                 conf.int = cint)
        if(n.boots>0){# && length(bCI)>0){
#            names(boot.cint)<-c(paste(as.character(conf.level), "\% bootstrap confidence interval: chpt.l",sep="\0"), "chpt.u")
            RVAL <- list(
                    data.name = DNAME,
                    parameter = para,
                    alternative = alternative,
                    statistic = stats,
                    change.point = chpt,
                    estimate = est,
                    p.value = p.value)#, 
#                    conf.int = c(cint,boot.cint))
        }
    }
    return(RVAL)
}
