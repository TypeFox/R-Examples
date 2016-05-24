"mreg" <-
function(formula, data, patid, start.theta=NULL, modify=unity, modify.p=0, mod.formula=~1, density.name="negbin",link="log", iterlim=100, gradtol=1e-6,steptol=1e-6, na.action=NULL,print.level=2, zero.start=FALSE){
	cl<-match.call()  
joint.formula <- form.add(substitute(formula), substitute(mod.formula))
        mf <- model.frame( joint.formula, data=data, na.action=na.action)
        Y <- model.response(mf)
        Z <- model.matrix(formula,data=mf)
        mod.Z <- model.matrix(mod.formula, data=mf)
        
        off <- model.offset(mf)
        if(is.null(off)){ off <- rep(0, length(Y))}
        pat <- data[,deparse(substitute(patid)),drop=TRUE]
        p <- dim(Z)[2]
        nuisance.p <- family.dim(density.name)
        density.name <- as.character(substitute(density.name))
        if( !is.null( start.theta)){
          if( length(start.theta)!= p+nuisance.p+modify.p){
            stop("incorrect length of starting parameters")
          }
        }
        
        link <- as.character(substitute(link))
        modify <- as.character(substitute(modify))
        if( is.null(start.theta)){
          start.theta <- rep(0,nuisance.p+p+modify.p)
        }
          answer <- nlm(minimand, start.theta, hessian=TRUE, Y=Y,Z=Z,off=off, pat=pat, print.level=print.level, nuisance.p=nuisance.p, iterlim=iterlim,gradtol=gradtol,steptol=steptol, modify=modify,density.name=density.name,link=link, mod.Z=mod.Z,zero.start=zero.start)
        
        answer$contrasts<-attr(Z,"contrasts")
        coef.index <- (nuisance.p+1):(nuisance.p+p+modify.p)
	dn<-c(colnames(Z),
              attributes(eval(call(modify, x=0, y=rep(0, modify.p),mod.Z=mod.Z[1,] )))$par.names)
	answer$coefficients<-answer$estimate[coef.index]
	names(answer$coefficients)<-dn
	answer$se<-sqrt( diag( solve(answer$hessian)))[coef.index]
	names(answer$se)<-dn
	answer$deviance<-2*answer$minimum
	answer$call<-cl
	answer$model<-mf
        answer$Z <- Z
        answer$offset<-off
        if( modify=="unity"){
          answer$linear.predictors<-Z%*%answer$coefficients[1:p]+off
          answer$fitted.values<-linkfn(answer$linear.predictors, link)
          answer$residuals<-Y-answer$fitted.values
        }
	answer$formula<-formula
	answer$y<-Y
	answer$density.name <- density.name
	answer$link <- link
        answer$modify <- modify
        if( nuisance.p>0){
        answer$nuisance<-list( estimate=answer$estimate[1:nuisance.p], 
                                  se=sqrt( diag( solve(answer$hessian)))[1:nuisance.p])
        }
	class(answer)<-c("mreg","lm")
	answer
	
}

