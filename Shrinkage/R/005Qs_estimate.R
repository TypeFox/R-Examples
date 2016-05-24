
# Zahra Montazeri, February 20, 2009
# M. Padilla, Oct. 2013 (modifications)
#--------------------
# x is always treatment
# y is always control
k.xy.prnset2matrix<-function(x,y=NULL,paired=F){
    if(is(x,"nxprnSet")){x<-nexprs(x)}
    else if(is(x,"nxprnSetPair")){zaux<-nexprs(x);x<-zaux$x;y<-zaux$y}
    assert.is(x,"matrix")
    if(paired){nmvar.x<-nmvar.y<-c("X","I")} else {nmvar.x<-c("X","I");nmvar.y<-c("X","S")}
    x<-MakeNames(x,nmvar=nmvar.x)
    
    if(length(y)>0){
       if(is(y,"nxprnSet")){y<-nexprs(y)}
        assert.is(y,"matrix")
        y<-MakeNames(y,nmvar=nmvar.y)
        if(nrow(x)!=nrow(y)) {
            message ("number of features (rows no. in x and y) are not same: keeping the common ones")
            fint<-intersect(rownames(x),rownames(y))
            if(length(fint)==0){stop ("no common sample names betwween x and y")}
            }
        else{
            fint<-rownames(x)
        }
        if(!identical(rownames(y),fint)){y<-y[fint,,drop=F]}
        if(paired){
            fint<-intersect(colnames(x),colnames(y))
            if(length(fint)==0){stop ("no common sample names betwween x and y")}
            y<-y[,fint,drop=F];x<-x[,fint,drop=F]
        }

        if(!identical(rownames(x),rownames(y))) {stop ("different feature names for x and y")}
    }
    list(x=x,y=y,paired=paired)
}

nQs.est<-function(x,y=NULL,opt="Q1", mu0 = 0, c = 0.5, a = 0.4, b = 0.01, h = 1){opt<-opt[1]
    k.x.find.q1<-function(x,mu0=mu0,h=h,c=c){
        n<-length(x)-sum(is.na(x))
        xbar<-mean(x, na.rm=TRUE)
        s2<-var(x,na.rm=TRUE)
        if((xbar-mu0)>0) (sgn=1)
        if((xbar-mu0)==0) (sgn=0)
        if((xbar-mu0)<0) (sgn=-1)
        R<-s2/(n*(xbar-mu0)^2)
        V<-(sqrt(s2)/(sqrt(n)*(1+R)))*sgn
        q1<-mu0+((xbar-mu0)/(1+h*R))+c*V
        q1}
    k.xy.find.q1.unpaired<-function(x,y,mu0,h,c){
        xbar <- mean(x,na.rm = TRUE)
        n1<-length(x)- sum(is.na(x))
        ybar <- mean(y,na.rm = TRUE)
        n2<-length(y) - sum(is.na(y))
        s2x<-var(x, na.rm=TRUE)
        s2y<-var(y, na.rm=TRUE)
        s2<-((n1-1)*s2x+(n2-1)*s2y)/n1+n1-2
        d<-xbar-ybar
        n<-n1+n2
        if((d-mu0)>0) (sgn=1)
        if((d-mu0)==0) (sgn=0)
        if((d-mu0)<0) (sgn=-1)
        R<-s2/(n*(d-mu0)^2)
        V<-(sqrt(s2)/(sqrt(n)*(1+R)))*sgn
        q1<-mu0+((d-mu0)/(1+h*R))+c*V
        q1}
    k.x.find.q2<-function(x,mu0,a,b,c) {
        n<-length(x)-sum(is.na(x))
        xbar<-mean(x, na.rm=TRUE)
        s2<-var(x,na.rm=TRUE)
        if((xbar-mu0)>0) (sgn=1)
        if((xbar-mu0)==0) (sgn=0)
        if((xbar-mu0)<0) (sgn=-1)
        R<-s2/(n*(xbar-mu0)^2)
        V<-(sqrt(s2)/(sqrt(n)*(1+R)))*sgn
        q2<-xbar-a*(xbar-mu0)*exp(-b/R)+c*V
        q2
        }
    k.xy.find.q2.unpaired<-function(x,y,mu0,a,b,c){
        xbar <- mean(x,na.rm = TRUE)
        n1<-length(x)- sum(is.na(x))
        ybar <- mean(y,na.rm = TRUE)
        n2<-length(y) - sum(is.na(y))
        s2x<-var(x, na.rm=TRUE)
        s2y<-var(y, na.rm=TRUE)
        s2<-((n1-1)*s2x+(n2-1)*s2y)/n1+n1-2
        d<-xbar-ybar
        n<-n1+n2
        if((d-mu0)>0) (sgn=1)
        if((d-mu0)==0) (sgn=0)
        if((d-mu0)<0) (sgn=-1)
        R<-s2/(n*(d-mu0)^2)
        V<-(sqrt(s2)/(sqrt(n)*(1+R)))*sgn
        q2<-d-a*(d-mu0)*exp(-b/R)+c*V
        q2}
    #---------------
    k.x<-function(x,fun,...){
        kfct<-function(i){
                zaux<-try(fun(x=x[i,,drop=T],...))
                if(is_err(zaux)){return(as.numeric(NA))}
                zaux}
            Q<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=kfct)
        names(Q)<-rownames(x)
        Q
    }
    k.xy<-function(x,y,fun,...){
         kfct<-function(i){
                zaux<-try(fun(x=x[i,,drop=T],y=y[i,,drop=T],...))
                if(is_err(zaux)){return(as.numeric(NA))}
                zaux}
            Q<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=kfct)
        names(Q)<-rownames(x)
        Q
    }
    #---------------
    z<-k.xy.prnset2matrix(x=x,y=y,paired=F)
    x<-z$x;y<-z$y
    
    if(length(y)==0){
        if(opt%in%c("Q1","q1","1")){
            Q<-k.x(x=x,fun=k.x.find.q1,mu0=mu0,h=h,c=c)
        }
        else if(opt%in%c("Q2","q2","2")){
            Q<-k.x(x=x,fun=k.x.find.q2,mu0=mu0,a=a,b=b,c=c)
        }
        else{stop("bad opt")}
    }
    else{
       
        if(opt%in%c("Q1","q1","1")){
            Q<-k.xy(x=x,y=y,fun=k.xy.find.q1.unpaired,mu0=mu0,h=h,c=c)
        }
        else if(opt%in%c("Q2","q2","2")){
            Q<-k.xy(x=x,y=y,fun=k.xy.find.q2.unpaired,mu0=mu0,a=a,b=b,c=c)
        }
        else{stop("bad opt")}
    }
    Q
    
}
#
nQ1.est<-function(x,y=NULL,mu0 = 0, c = 0.5, h = 1){
    nQs.est(x=x,y=y,opt="Q1",mu0=mu0,h=h,c=c)
}
nQ2.est<-function(x,y=NULL,mu0=0,c = 0.5, a = 0.4, b = 0.01){
    nQs.est(x=x,y=y,opt="Q2",mu0=mu0,a=a,b=b,c=c)
}
##-----------

