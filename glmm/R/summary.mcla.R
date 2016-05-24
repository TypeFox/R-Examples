summary.glmm <-
function(object,...){
    mod.mcml<-object
    stopifnot(inherits(mod.mcml, "glmm"))

	trust.converged<-mod.mcml$trust.converged
	fixedcall<-mod.mcml$fixedcall
	randcall<-mod.mcml$randcall
	call<-mod.mcml$call
	x<-mod.mcml$x
	y<-mod.mcml$y
	z<-mod.mcml$z
	
	#the coefficients matrix for fixed effects
	beta<-mod.mcml$beta
	nbeta<-length(beta)
	hessian<-mod.mcml$likelihood.hessian
	if(det(hessian)==0) {
            warning(paste("estimated Fisher information matrix not positive",
               "definite, making all standard errors infinite"))
            all.ses <- rep(Inf, nrow(hessian))
        }
	else{	all.ses<-sqrt(diag(solve(-hessian)))}
	
	beta.se<-all.ses[1:nbeta]
	zval<-beta/beta.se
	coefmat<-cbind(beta,beta.se,zval,2*pnorm(abs(zval),lower.tail=F))
	colnames(coefmat)<-c("Estimate","Std. Error", "z value", "Pr(>|z|)")
	rownames(coefmat)<-colnames(mod.mcml$x)
	coefficients<-coefmat[,1]
	
	nu<-mod.mcml$nu
	nu.se<-all.ses[-(1:nbeta)]
	nuzval<-nu/nu.se
	nucoefmat<-cbind(nu,nu.se,nuzval,pnorm(abs(nuzval),lower.tail=F))
	colnames(nucoefmat)<-c("Estimate","Std. Error", "z value", "Pr(>|z|)/2")
	rownames(nucoefmat)<-mod.mcml$varcomps.names
	link<-mod.mcml$family.glmm$link

	
	return(structure(list(x=x,y=y, z=z, coefmat=coefmat, fixedcall = fixedcall, randcall = randcall, coefficients = coefficients,
family.mcml = mod.mcml$family.mcml, call = call, nucoefmat = nucoefmat, link = link, trust.converged = trust.converged),class="summary.glmm"))
}

print.summary.glmm <-
    function (x, digits = max(3, getOption("digits") - 3),
        signif.stars = getOption("show.signif.stars"), ...)
{
    summ<-x	
    stopifnot(inherits(summ, "summary.glmm"))

	if(summ$trust.converged==FALSE)  cat("\nWARNING: the optimizer trust has not converged to the MCMLE. The following estimates are not maximum likelihood estimates, but they can be used in the argument par.init when rerunning glmm. \n\n",sep="")

    cat("\nCall:\n", paste(deparse(summ$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

 cat("\nLink is: ", paste(deparse(summ$link)),"\n\n",sep="")
   
cat("Fixed Effects:")
   cat("\n")

	printCoefmat(summ$coefmat, digits = digits,
        signif.stars = signif.stars, na.print = "NA", ...)
   cat("\n")

   cat("\n")
cat("Variance Components for Random Effects (P-values are one-tailed):")
   cat("\n")

	printCoefmat(summ$nucoefmat, digits = digits,
        signif.stars = signif.stars, na.print = "NA", ...)
   cat("\n")

}



coef.glmm <-
function(object,...){
	mod<-object
   	stopifnot(inherits(mod, "glmm"))
	coefficients<-mod$beta
	names(coefficients)<-colnames(mod$x)
	coefficients
}


vcov.glmm <-
function(object,...){
	mod<-object
   	stopifnot(inherits(mod, "glmm"))
	vcov <- -solve(mod$likelihood.hessian)

	#get names for vcov matrix
	rownames(vcov)<-colnames(vcov)<-rep(c("blah"),nrow(vcov))
	nbeta<-length(mod$beta)
	rownames(vcov)[1:nbeta] <- colnames(vcov)[1:nbeta] <- colnames(mod$x)
	rownames(vcov)[-(1:nbeta)] <- colnames(vcov)[-(1:nbeta)] <- mod$varcomps.names

	vcov
}

varcomps<-function(object,...){

	mod<-object
   	stopifnot(inherits(mod, "glmm"))
	coefficients<-mod$nu
	coefficients
}

confint.glmm<-function(object,parm,level=.95,...){
   	stopifnot(inherits(object, "glmm"))

	cf<-c(coef(object),varcomps(object))
	pnames<-names(cf)
	if(missing(parm)) parm<-pnames
	else if (is.numeric(parm)) parm<-pnames[parm]
	stopifnot(parm %in% pnames)

	a<-(1-level)/2
	a<-c(a,1-a)
	fac<-qnorm(a)
	pct<-a*100
	
	betaandnu<-c(object$beta,	nu<-object$nu)
	hessian<-object$likelihood.hessian
	if(det(hessian)==0) {
            warning(paste("estimated Fisher information matrix not positive",
               "definite, making all standard errors infinite"))
            all.ses <- rep(Inf, nrow(hessian))
        }
	else{	all.ses<-sqrt(diag(solve(-hessian)))}
	
	names(all.ses)<-pnames
	ci<-matrix(data=NA,nrow=length(parm),ncol=2)
	colnames(ci)<-a
	rownames(ci)<-parm
	ci[,1]<-betaandnu[parm]+a[1]*all.ses[parm]
	ci[,2]<-betaandnu[parm]+a[2]*all.ses[parm]

	ci
}

