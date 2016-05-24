ignorableWarnings=function(w)
{
	harmless=c(
		"step size truncated due to divergence",
		"step size truncated: out of bounds", 
		"step size truncated due to increasing deviance", 
		"inner loop 3: cannot correct step size",
		"inner loop 3; cannot correct step size",
		"glm.fit: algorithm did not converge",
		"glm.fit: algorithm stopped at boundary value",
		"glm.fit3: algorithm did not converge. Try increasing the maximum iterations",
		"glm.fit3: algorithm stopped at boundary value", 
		"glm.fit3: fitted rates numerically 0 occurred"
	)
	if (!(w$message %in% harmless)) return(w)
	
	invokeRestart("muffleWarning")
}

glmSolvers=list(
	glm=list(solve=expression(withCallingHandlers(glm(formula=formula, family=family, data=data, control=control, ...),simpleWarning=ignorableWarnings	)),
			 success=expression(ans$converged), 
			 result=expression(ans)),
	glm.fit3=list(solve=expression(withCallingHandlers(glm(formula=formula, family=family, data=data, control=control, method=glm.fit3, ...),simpleWarning=ignorableWarnings)),
			 success=expression(ans$converged), 
			 result=expression(ans)), 
	nlminb=list(solve=expression({
		this.ctrl=control; 
		this.ctrl$maxit=1L; this.initGlm=withCallingHandlers(glm(formula=formula, family=family, data=data, control=this.ctrl, ...), simpleWarning=ignorableWarnings);
		this.wt=this.initGlm$weights; 
		this.idxn=which(this.wt>0); this.idxp=which(!is.na(this.initGlm$coefficients)); 
		this.off=if(is.null(this.initGlm$offset)) 0 else this.initGlm$offset; 
		this.x=model.matrix(this.initGlm);  		
		dev.func=function(bet){
			sum(family$dev.resids(this.initGlm$y[this.idxn], family$linkinv(this.x[this.idxn,this.idxp,drop=FALSE]%*%bet[this.idxp]+ this.off[this.idxn]), this.wt[this.idxn]))
		}
		nlminb(this.initGlm$coefficients, dev.func, control=list(iter.max=control$maxit, eval.max=round(control$maxit/3L*4L), abs.tol=control$epsilon^2, rel.tol=control$epsilon, x.tol=control$epsilon, trace=FALSE))
		}), 
		success=expression(ans$convergence==0), 
		result=expression({
			this.rslt=ans$par; this.rslt[is.na(this.rslt)]=0; 
			this.ans=withCallingHandlers(tryCatch(glm(formula=formula, family=family, data=data, control=this.ctrl, start=this.rslt, method=glm.fit3, ...), simpleError=function(message, call = NULL)message), simpleWarning=ignorableWarnings)
			if(!inherits(this.ans, 'glm')) next
			this.ans
		})), 
	BFGS=list(solve=expression({
		this.ctrl=control; 
		this.ctrl$maxit=1L; this.initGlm=withCallingHandlers(glm(formula=formula, family=family, data=data, control=this.ctrl, ...), simpleWarning=ignorableWarnings);
		this.wt=this.initGlm$weights; 
		this.idxn=which(this.wt>0); this.idxp=which(!is.na(this.initGlm$coefficients)); 
		this.off=if(is.null(this.initGlm$offset)) 0 else this.initGlm$offset; 
		this.x=model.matrix(this.initGlm);  		
		dev.func=function(bet){
			sum(family$dev.resids(this.initGlm$y[this.idxn], family$linkinv(this.x[this.idxn,this.idxp,drop=FALSE]%*%bet[this.idxp]+ this.off[this.idxn]), this.wt[this.idxn]))
		}
		this.start=this.initGlm$coefficients; if(length(this.idxp)<length(this.start)) this.start[-this.idxp]=0
		optim(this.start, dev.func, method='BFGS', control=list(maxit=control$maxit, reltol=control$epsilon, trace=FALSE))
		}), 
		success=expression(ans$convergence==0), 
		result=expression({
			this.rslt=ans$par; this.rslt[is.na(this.rslt)]=0 ;
			this.ans=withCallingHandlers(tryCatch(glm(formula=formula, family=family, data=data, control=this.ctrl, start=this.rslt, method=glm.fit3, ...), simpleError=function(message, call = NULL)message), simpleWarning=ignorableWarnings)
			if(!inherits(this.ans, 'glm'))	next
			this.ans
		}))
)

glmsolve=function(formula, family, data, control, solvers=glmSolvers, ...)
{
	nsolvers=length(solvers)
	if(nsolvers==0L || !is.list(solvers))stop('solvers must be supplied as a list')
	solver.names=names(solvers)
	if(is.null(solver.names)) solver.names=paste('anonymous',seq_len(nsolvers),sep='.')
	tmp=solver.names==''
	if(any(tmp)) solver.names[tmp]=paste('anonymous',seq_len(sum(tmp)),sep='.')
	
	thisEnv=sys.frame(sys.nframe())
	lastNonError=NULL
	for(i in seq_len(nsolvers)){
		solver.name=solver.names[i]
		ans=try( eval(solvers[[i]]$solve, envir=thisEnv), silent=TRUE)
		if(inherits(ans,'try-error')) next 
		success=eval(solvers[[i]]$success, envir=thisEnv)
		ans=eval(solvers[[i]]$result, envir=thisEnv)
		attr(ans, 'glmsolve.method')=solver.name
		if(isTRUE(success)){
			attr(ans, 'glmsolve.success')=TRUE
			return(ans)
		}else {
			lastNonError=ans
			attr(lastNonError, 'glmsolve.success')=FALSE
		}
	}
	
	family=tryCatch(reinitialize.fbrNBfamily(family), simpleWarning=function(w)w$message)
	if(inherits(family, 'fbrNBfamily')) 
		return(Recall(formula, family, data, control, solvers, ...))
	
	do.call(if(is.null(lastNonError)) 'stop' else 'warning', list('None of the glm solvers succeeded.'))
	lastNonError
}