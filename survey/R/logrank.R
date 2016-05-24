svylogrank<-function(formula, design,rho=0,gamma=0,method=c("small","large","score"),...){
	UseMethod("svylogrank",design)
}

print.svylogrank<-function(x,...){
	m<-t(x)
	rownames(m)=""
	printCoefmat(m,has.Pvalue=TRUE,P.values=TRUE)
	invisible(NULL)
	}
	
.logrank<-function(formula, design,rho=0,gamma=0){
	nullformula<-update(formula,.~1)
	S<-svykm(nullformula,design,se=FALSE)
	epsilon<-min(diff(sort(unique(S$time))))/10
	w<-approxfun(S$time+epsilon,S$surv^rho*(1-S$surv)^gamma,method="constant",rule=2)
	environment(formula)<-environment()
	coxmodel<-coxph(formula,data=model.frame(design), weights=weights(design,"sampling"),iter.max=0)
	x<-model.matrix(coxmodel)

	detail<-coxph.detail(coxmodel,riskmat=TRUE)	
	Y<-t(detail$riskmat)
	dLambda<-detail$hazard
	E<-as.matrix(detail$means)
	N<-coxmodel$y[,"status"]

	times<-coxmodel$y[,"time"]
	U<-matrix(nrow=nrow(x),ncol=ncol(x))
	index<-match(times[N==1],detail$time)
	ZmEdN<- matrix(0,nrow=nrow(x),ncol=ncol(x))
	ZmEdN[N==1,]<-x[N==1,,drop=FALSE]-E[index,]
	for(p in 1:ncol(x)){
            ZmE <- -outer(E[,p], x[,p], "-")  ##times are rows, people are columns
            U[,p]<- ZmEdN[,p]*w(times)- colSums(w(detail$time)*ZmE*dLambda*Y)
	}
	means <- svytotal(U,design)
	zstat<-coef(means)/SE(means)
	chisqstat<-coef(means)%*%solve(vcov(means),coef(means))

	rval<-list(cbind(score=coef(means),se=SE(means),z=coef(means)/SE(means),p= 2*pnorm(-abs(coef(means)/SE(means)))),
		 c(chisq=chisqstat,p=pchisq(chisqstat,df=ncol(x),lower.tail=FALSE)))
        class(rval)<-"svylogrank"
	rval
	}
	

.biglogrank<-function(formula, design,rho=0,gamma=0){
	nullformula<-update(formula,.~1)
	S<-svykm(nullformula,design,se=FALSE)
	epsilon<-min(diff(sort(unique(S$time))))/10
	w<-approxfun(S$time+epsilon,S$surv^rho*(1-S$surv)^gamma,method="constant",rule=2)
	environment(formula)<-environment()
	coxmodel<-coxph(formula,data=model.frame(design), weights=weights(design,"sampling"),iter.max=0)
	x<-model.matrix(coxmodel)

	detail<-coxph.detail(coxmodel)	

	dLambda<-detail$hazard
	E<-as.matrix(detail$means)
	N<-coxmodel$y[,"status"]

	times<-coxmodel$y[,"time"]
	U<-matrix(nrow=nrow(x),ncol=ncol(x))
	index<-match(times[N==1],detail$time)
	ZmEdN<- matrix(0,nrow=nrow(x),ncol=ncol(x))
	ZmEdN[N==1,]<-x[N==1,,drop=FALSE]-E[index,]
	for(p in 1:ncol(x)){
		U[,p]<- ZmEdN[,p]*w(times)
		for (j in seq_along(detail$time)){
		  thistime<-detail$time[j]
		  ZmE <- x[,p]-E[j,p]
		  U[,p] <- U[,p]  - w(thistime)*ZmE*dLambda[j]*(times>=thistime)
		}
	}
	means <- svytotal(U,design)
	zstat<-coef(means)/SE(means)
	chisqstat<-coef(means)%*%solve(vcov(means),coef(means))

	rval<-list(data.frame(score=coef(means),se=SE(means),z=coef(means)/SE(means),p= 2*pnorm(-abs(coef(means)/SE(means)))),
		 c(chisq=chisqstat,p=pchisq(chisqstat,df=ncol(x),lower.tail=FALSE)))
        class(rval)<-"svylogrank"
        rval
	
	}
		
svylogrank.survey.design2<-function(formula, design,rho=0,gamma=0,
                                    method=c("small","large","score"),
                                    ...){
	require(survival) || stop("requires the survival package")
        method<-match.arg(method)
        if (method=="small")
            return(.logrank(formula,design, rho,gamma,...))
        else if (method=="large")
            return(.biglogrank(formula,design,rho,gamma,...))
        if (rho!=0 || gamma!=0){
            return(expandlogrank(formula,design,rho,gamma,...))
        }
        
	tms<-delete.response(terms(formula,specials="strata"))
	findstrat<-untangle.specials(tms,"strata")
	if(length(findstrat$terms))
	   tms<-tms[-findstrat$terms]
	mf<-model.frame(tms,model.frame(design))
	if(length(mf)>1)
	   stop("Only one grouping variable allowed")
	if(!is.factor(mf[[1]]) && length(unique(mf[[1]]))>2)
	   stop("Grouping variable with more than 2 levels must be a factor")
		
	b<-coef(svycoxph(formula,design,iter=1))
	v<-vcov(svycoxph(formula,design,iter=0))
	x2<-sum(b*solve(v,b))
	rval<-c(z=b/sqrt(diag(v)), Chisq=x2, p=pchisq(x2,length(b),lower.tail=FALSE))
	class(rval)<-"svylogrank"
	rval
	}

svylogrank.twophase<-svylogrank.survey.design2
svylogrank.twophase2<-svylogrank.survey.design2

svylogrank.DBIsvydesign<-function (formula, design, ...) 
{
    design$variables <- dropFactor(getvars(formula, design$db$connection, 
        design$db$tablename, updates = design$updates, subset = design$subset), 
        weights(design))
    NextMethod("svylogrank", design)
}

svylogrank.ODBCsvydesign<-function (formula, design, ...) 
{
    design$variables <- dropFactor(getvars(formula, design$db$connection, 
        design$db$tablename, updates = design$updates), weights(design))
    NextMethod("svylogrank", design)
}	
	
svylogrank.svyrep.design<-function(formula, design,rho=0,gamma=0,method=c("small","large","score"), ...){
	require(survival) || stop("requires the survival package")
        method<-match.arg(method)
        if (method=="small")
            return(.logrank(formula,design, rho,gamma,...))
        else if (method=="large")
            return(.biglogrank(formula,design,rho,gamma,...))
        if (rho!=0 || gamma!=0){
            return(expandlogrank(formula,design,rho,gamma,...))
        }
	tms<-delete.response(terms(formula,specials="strata"))
	findstrat<-untangle.specials(tms,"strata")
	if(length(findstrat$terms))
	   tms<-tms[-findstrat$terms]
	mf<-model.frame(tms,model.frame(design))
	if(length(mf)>1)
	   stop("Only one grouping variable allowed")
	if(!is.factor(mf[[1]]) && length(unique(mf[[1]]))>2)
	   stop("Grouping variable with more than 2 levels must be a factor")
  
	rr<-withReplicates(design, function(w,df){
		  environment(formula)<-environment()
		  coef(coxph(formula,data=df,weights=w+1e-8,iter=1))
		})
	
   b<-unclass(rr)
   attr(b,"var")<-NULL
	v<-attr(rr,"var")
	x2<-sum(b*solve(v,b))
   rval<- c(z=b/sqrt(diag(as.matrix(v))), Chisq=x2, p=pchisq(x2,length(b),lower.tail=FALSE))
   class(rval)<-"svylogrank"
	rval
	}	



expandlogrank<-function(formula, design, rho=0, gamma=0){
	nullformula<-update(formula,.~1)
	S<-svykm(nullformula,design,se=FALSE)
 	epsilon<-min(diff(sort(unique(S$time))))/10
	w<-approxfun(S$time+epsilon,S$surv^rho*(1-S$surv)^gamma,method="constant",rule=2)
	environment(formula)<-environment()
	coxmodel<-coxph(formula,data=model.frame(design), weights=weights(design,"sampling"),iter.max=0)
	mf<-model.frame(design)
	detail<-coxph.detail(coxmodel)	
	
	if(attr(coxmodel$y,"type")=="right"){
		mf$.time<-coxmodel$y[,"time"]
		mf$.status<-coxmodel$y[,"status"]
		mfsplit <- survSplit(mf, cut=detail$time, end=".time", event=".status", start=".start", id=".id", episode=".episode")
		} else {
		mf$.start<-coxmodel$y[,"start"]
		mf$.time<-coxmodel$y[,"stop"]
		mf$.status<-coxmodel$y[,"status"]
		mfsplit <- survSplit(mf, cut=detail$time, end=".time", event=".status", start=".start", id=".id", episode=".episode")
		}	
	
	formula[[2]]<-quote(Surv(.start,.time,.status))

	mfsplit$.weights<-weights(design,"sampling")[match(mfsplit$.id, rownames(mf))]*w(mfsplit$.time)
	expdesign<-svydesign(ids=eval(design$call$id), strata=eval(design$call$strata), data=mfsplit, weights=~.weights)
	#svylogrank(formula,expdesign)
	tms<-delete.response(terms(formula,specials="strata"))
	findstrat<-untangle.specials(tms,"strata")
	if(length(findstrat$terms))
	   tms<-tms[-findstrat$terms]
	mf<-model.frame(tms,model.frame(expdesign))
	if(length(mf)>1)
	   stop("Only one grouping variable allowed")
	if(!is.factor(mf[[1]]) && length(unique(mf[[1]]))>2)
	   stop("Grouping variable with more than 2 levels must be a factor")
		
	b<-coef(svycoxph(formula,expdesign,iter=1))
	v<-vcov(svycoxph(formula,expdesign,iter=0))
	x2<-sum(b*solve(v,b))
	rval<-c(z=b/sqrt(diag(v)), Chisq=x2, p=pchisq(x2,length(b),lower.tail=FALSE))
	class(rval)<-"svylogrank"
	rval
	}
	
