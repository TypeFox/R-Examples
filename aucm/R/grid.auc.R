# eval auc functions at a grid of betas
# approx.param is s for rauc and h for sauc
# approx.param=1;lambda=0;loss=TRUE;t0=NULL; t1=NULL; approx.type=NULL # default
grid.auc=function(formula, dat, beta, approx.type=NULL, approx.param=1, lambda=0, loss=TRUE, t0=NULL, t1=NULL, ret.vcov=FALSE){
    
    if (lambda>0) loss=TRUE
    
    tmp=model.frame(formula, dat)
    x1=model.matrix(formula, tmp[tmp[,1]==1,])[,-1]
    x2=model.matrix(formula, tmp[tmp[,1]==0,])[,-1]
    n1=nrow(x1)
    n2=nrow(x2)
    dat.y=tmp[,1]
    X=as.matrix(tmp[,-1,drop=FALSE])
    
    pAUC=ifelse(!is.null(t0) | !is.null(t1), TRUE, FALSE)
    if (pAUC) {
        if (is.null(t0)) t0=0
        if (is.null(t1)) t1=1            
    }
    
    
    if ( ncol(beta) == 1 ) {
        beta1=cbind(1,beta)
    } else {
        beta1 = beta
    }
        
    aux=apply(beta1, 1, function(coef.){
    
        if (pAUC) {            
            eta.j.init=drop(x2 %*% coef.)            
            ind2=order(eta.j.init)[ceiling((1-t1)*n2):floor((1-t0)*n2)] 
            A=rep(FALSE,n2)
            A[ind2]=TRUE            
        } else {
            A=rep(TRUE, n2)         
        }            
        
        lin=X%*%coef.
        dat1=lin[dat.y==1]
        dat0=lin[dat.y==0] [A] 
        
        if (is.null(approx.type)) {
            inner=outer(dat1, dat0, "-")
            tmp=mean(inner>0)
            #tmp=fast.auc(lin, dat.y)
            if (loss) tmp=1-tmp
            out=tmp
            if (pAUC) out=out*length(dat1)*length(dat0)/(n1*n2)
            c(out, 1)
        } else if (approx.type %in% c("rauc","normal","logistic")) {
            inner=outer(dat1, dat0, "-")
            if (approx.type=="rauc") {
                tmp=ramp.f(inner,approx.param,loss=loss)
            } else if (approx.type=="normal"){
                tmp=phi.f(inner,approx.param,loss=loss)
            } else if (approx.type=="logistic"){
                tmp=logistic.f(inner,approx.param,loss=loss)
            }

            out=mean(tmp) 
            if (pAUC) out=out*nrow(tmp)*ncol(tmp)/(n1*n2)
            bounded=mean(tmp==1 | tmp==0)
            c(out, bounded)
        } else stop("approx.type not supported: "%+% approx.type)
        
        
    }) 
    
    vals=aux[1,]*n1*n2 + .5*lambda*apply(beta, 1, function(x) sum(x**2))
    if (loss) {
        val=min(vals)
        par=beta[which.min(vals),]
    } else {
        val=max(vals)
        par=beta[which.max(vals),]
    }
    
    out=list("coefficient"=par, "val"=val/n1/n2, "pars"=t(beta1), "vals"=vals/n1/n2, "bounded.ratio"=aux[2,], "formula"=formula)
        
    if (ret.vcov) out$vcov = get.vcov.sauc (sigm.f=ifelse(is.null(approx.type),"logistic",approx.type), beta.hat=par, h=approx.param, X1=x1, X2=x2)     
    
    class(out)=c("auc","gridsearch", class(out))
    return(out)
    
}


vcov.gridsearch=function(object, ...) {
    object$vcov
}

plot.gridsearch=function(x, type="l", xlab="par", ylab="val", ...){
    plot(x$pars[2,], x$vals, type=type, xlab=xlab, ylab=ylab, ...)
}
