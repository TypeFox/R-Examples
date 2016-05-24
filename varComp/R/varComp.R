nlminb.control=function(
	eval.max=200L, 
	iter.max=150L,
	trace = 0L, 
	abs.tol = 0, 
	rel.tol = 1e-10, 
	x.tol = 1.5e-8, 
	xf.tol = 2.2e-14, 
	step.min = 1, 
	step.max = 1, 
	sing.tol =rel.tol,
	...)
{structure(list(
	eval.max=eval.max, 
	iter.max=iter.max,
	trace = trace, 
	abs.tol = abs.tol, 
	rel.tol = rel.tol, 
	x.tol = x.tol, 
	xf.tol = xf.tol, 
	step.min = step.min, 
	step.max = step.max, 
	sing.tol =sing.tol, 
	...) 
	, class = 'nlminb.control')
}

varComp.control=function(verbose = FALSE, start=NULL, REML = TRUE, 
information=informationTypes, boundary.eps=5e-4, 
  nlminb=nlminb.control(iter.max=200L, eval.max=500L), 
  plot.it=FALSE, keepXYK=TRUE)
{	optMethods=c('nlminb', 'optim', 'NRGD')
	optMethod=optMethods[match('nlminb', optMethods)]
	structure(
	  list(optMethod=optMethod,
		 verbose=verbose, 
		 starts=start, 
		 REML = REML, 
		 information = match.arg(information, informationTypes), 
		 boundary.eps= boundary.eps, 
		 nlminb = nlminb, 
		 plot.it=plot.it, 
		 keepXYK=keepXYK) #nStepHalving
	  , class = 'varComp.control'
	)
}

varComp=function(fixed, data, random, varcov, weights, subset, family = stats::gaussian('identity'), na.action, offset, control = varComp.control(...), 
     doFit = TRUE, normalizeTrace = TRUE, 
     contrasts = NULL, model = TRUE, X = TRUE, Y = TRUE, K = TRUE, ...)
{
	if(missing(fixed) || !is.formula(fixed) || length(as.list(fixed))!=3L)	stop('fixed needs to be a two-sided formula')
	if(!missing(random)) {
		if(!is.formula(random) || length(as.list(random))!=2L)  stop('random needs to be a right-sided formula, when non-missing')
		if(!identical(environment(fixed) , environment(random)) )warning("environment(random) does not match environment(fixed).")
		fixedRandom = update.formula(fixed, as.formula(paste0('.~.+',random[2],collapse=''))) # random[2] is a call to the RHS
		random = update.formula(random, as.formula('~.')) ## this will expand the formula (removing * shorthand)
	}else fixedRandom = fixed 
	
	if(!missing(varcov) && !is.matrix(varcov) && !is.list(varcov)) stop("varcov needs to be a matrix or a list, when non-missing")
	if(!missing(varcov) && is.matrix(varcov)) varcov=list(varcov)
    if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family))  family <- family()
	if(!all.equal(stats::gaussian('identity'), family)) { 
		warning("Currently only Gaussian family with identity link is supported")
		family=stats::gaussian(link='identity')
	}
	if (missing(data)) data <- environment(formula)
	ret.x = X; ret.y=Y; ret.k=K; ret.mod=model
	
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)

    m <- match(c("fixed", "data", "random", "weights" ,"subset", "na.action", "offset"
			), names(mf), 0L)
    mf <- mf[c(1L, m)]  ## mf[1] is varComp()
	
	if(m[5L]>0 && is.logical(subset)) {
		subset=which(subset)
		mf$subset=subset
	}
	if(m[4L]>0 && any(weights<=0)) {	## handling non-positive weights to subset argument
		if(m[5L]>0) subset = unique(setdiff(subset, which(weights<=0))) else subset = which(weights>0)
		mf$subset=subset
	}
	if(m[3L]>0L) mf['random']=NULL
	mf['fixed']=NULL
	mf$formula=fixedRandom		
	
	mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mfAll <- eval(mf, parent.frame())	## removing obs with missing data or not in subset
	
	if(!is.null( na.vec <- attr(mfAll, 'na.action'))){	## missing data removed. 
		if(m[5L]> 0 ){ ## there exists subset argument
			subset = sort(setdiff(subset, na.vec))
		}else   subset = sort(setdiff(seq_len(nrow(mfAll)+length(na.vec)), na.vec))
		mf$subset = subset
		mfAll = eval(mf, parent.frame())
	}
	mtAll = attr(mfAll, 'terms')
	
	# mf$data = as.name('mfAll')
	mf$formula=fixed
	eframe=sys.frame(sys.nframe())
	mf.fixed=eval(mf, parent.frame())
	
    mt.fixed <- attr(mf.fixed, "terms")
	fixedTerms = attr(mt.fixed, 'term.labels')
	fixedTerms = sapply(fixedTerms, sortTerm)
	fixedVars  = sapply(attr(mt.fixed, 'variables'), as.character)[-1L]
    Y <- model.response(mf.fixed, "numeric")
    w <- as.vector(model.weights(mf.fixed))
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf.fixed))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
		Y = Y - offset  ## CHECKME when family()$link is not identity
    }
    if (is.empty.model(mt.fixed)) {
        X <- matrix(0, length(Y), 0L)
    } else {
        X <- model.matrix(mt.fixed, mf.fixed, contrasts)
    }
	
	if(missing(random)){	
		if(missing(varcov)) varcov=list()
		K=varcov
		nK=length(K)
	}else{
	## preparing a vector of all random terms
		mf$formula=random
		# if(m[2L]>0L) mf$data = cl$data
		mf.random=eval(mf, parent.frame())
		mt.random=attr(mf.random, 'terms')
		rterms = attr(mt.random, 'term.labels')
		# rterms = unlist(strsplit(as.character(random[2]), ' *\\+ *'))
		rterms = setdiff(rterms, -1:1) ## removing all terms invoving intercept(s)
		rterms = sapply(rterms, sortTerm, priority=fixedVars)
		rterms = setdiff(rterms, fixedTerms)
		rterms = local({  ## expand all fixed-by-random interactions to ensure that later when fixed-by-random interactions are encountered, the marginal random effect does not need to be added. 
			ans=rterms
			for(iz in seq_along(rterms)){
				this.form=update.formula(random, as.formula(paste('~', rterms[iz], '+0', collapse='')))
				mf$formula = this.form
				mf.this <- eval(mf, parent.frame())
				# isFact <- vapply(mf.this, function(x) is.factor(x) || is.logical(x), NA)
				isFixed = names(mf.this)%in%fixedVars
				# idx = (isFact & isFixed)
				idx = isFixed
				if(any(idx) && !all(idx)) ans[iz] = gsub(":", "*", ans[iz], fixed=TRUE)
			}
			ans
		})
		mf$formula = update.formula(random, as.formula(paste0("~", paste0(rterms, collapse='+'))))
		mf.random=eval(mf, parent.frame())
		mt.random=attr(mf.random, 'terms')
		rterms = attr(mt.random, 'term.labels')
		rterms = setdiff(rterms, -1:1) ## removing all terms involving intercept(s)
		rterms = sapply(rterms, sortTerm, priority=fixedVars)
		rterms = setdiff(rterms, fixedTerms)

		#### processing each random term
		nRterms = length(rterms)
		Z=vector('list', nRterms)
		j=1L
		for(iz in seq_len(nRterms)){
			this.form=update.formula(random, as.formula(paste('~', rterms[iz], '+0', collapse='')))
			mf$formula = this.form
			mf.this <- eval(mf, parent.frame())
			mt.this <- attr(mf.this, "terms")
			tmp.model.matrix = suppressWarnings(model.matrix(mt.this, mf.this, contrasts))
			isFact <- vapply(mf.this, function(x) is.factor(x) || is.logical(x), NA)
			isFixed = names(mf.this)%in%fixedVars
			idx = which(isFact & isFixed)
			if(length(idx)>0L){ ## random interactions involving fixed variables: previously, it has been ensured that random marginal effects have already been included
				fterm=paste0(sort(names(mf.this)[idx]), collapse=':')
				if(fterm%in%fixedTerms){ ## truly fixed-by-random interaction
					this.form=update.formula(this.form, paste0("~",fterm))
					mf$formula=this.form
					mf.tmp=eval(mf, parent.frame())
					mt.tmp=attr(mf.tmp, 'terms')
					oldContr=getOption('contrasts')
					options(contrasts=c('contr.sum','contr.sum'))
					fixedX = suppressWarnings(model.matrix(mt.tmp, mf.tmp, contrasts)[,-1L,drop=FALSE])
					options(contrasts=oldContr)
					
					tmpRterm = substr(rterms[iz], nchar(fterm)+2L, nchar(rterms[iz]))
					if(tmpRterm=='') {warning("DEBUG ME: tmpRterm should not be empty"); browser()}
					this.form=update.formula(this.form, paste0("~", tmpRterm, '+0'))
					mf$formula=this.form
					mf.tmp=eval(mf, parent.frame())
					mt.tmp=attr(mf.tmp, 'terms')
					rdmZ=suppressWarnings(model.matrix(mt.tmp, mf.tmp, contrasts))
					
					for(tmpi in seq_len(ncol(fixedX))){
						Z[[j]] = fixedX[, tmpi] * rdmZ
						names(Z)[j] = paste(colnames(fixedX)[tmpi], tmpRterm, sep=':')
						j=j+1L
					}
				}else{	## the fixed interaction is actually random
					Z[[j]] = tmp.model.matrix
					names(Z)[j] = rterms[iz]
					j=j+1L
				}				
			}else{	## no fixed terms involved
				Z[[j]] = tmp.model.matrix
				names(Z)[j] = rterms[iz]
				j = j+1L
			}
		}
		nRterms = length(Z)
		if(missing(varcov)){
		  nK=nRterms
		  K=vector('list', nK)
		  for(j in seq_len(nK))  K[[j]]=tcrossprod(Z[[j]])
		}else{
			K=varcov
			nK=length(K)
		  if(nRterms!=nK) stop('The number of matrices in "varcov" needs to equal the number of random effect terms in "random", when both are provided.')
		  ## match the order of elements in Z with that in K the same way as matching arguments in function calls
			tmpf=function(){mget(ls())}
			tmpz=seq_along(Z); namesZ=names(tmpz)=names(Z)
			formals(tmpf)=as.pairlist(tmpz)
			tmpk=as.list(seq_along(K)); names(tmpk)=names(K)
			zk.idx=do.call('tmpf', tmpk)
			K.bak=K
		  for(j in 1:nK)  K[[j]]=tcrossprod(Z[[j]]%*%K.bak[[ zk.idx[[ namesZ[j] ]] ]], Z[[j]])
		}
 	    names(K)=names(Z)
 	}
	
	if(isTRUE(normalizeTrace)) K=lapply(K, normalizeTrace)
	
	if(!missing(weights)){
		keep.weights= which(w>0)
		w.5 = sqrt(w[keep.weights])
		rw.5 = w.5[rep(seq_along(keep.weights), each=length(keep.weights))]
		Y=Y[keep.weights]*w.5
		X=w.5*X[keep.weights, , drop=FALSE]
		for(ik in seq_len(nK)) K[[ik]] = w.5*K[[ik]][keep.weights, keep.weights, drop=FALSE]*rw.5
	}
	
	if(any(ret.x, ret.y, ret.k)) control$keepXYK=TRUE
	ansCall = call('varComp.fit', Y=Y, X=X, K=K, control=control)
	
	
	ans = list(doFit = ansCall)
	ans$fixef = matrix(NA_real_, ncol(X), 0L, dimnames=list(colnames(X),NULL))
	
    class(ans) <- "varComp"
    ans['na.action'] <- list(attr(mfAll, "na.action"))  ## could be null
    ans['offset'] <- list(offset) ## could be null
    ans['contrasts'] <- list(attr(X, "contrasts"))  ## if X is empty, this could be null
    ans$xzlevels <- .getXlevels(mtAll, mfAll)
    ans$call <- cl
    ans$terms <- mtAll
	ans$nobs = length(Y)
	ans$control=control
    ## CHEKME: 
	if (ret.mod)  ans$model <- mfAll
    if (!ret.x)     ans$X <- NULL else ans$X=X
    if (!ret.y)     ans$Y = NULL else ans$Y = Y
	if (!ret.k) 	 ans$K = NULL else ans$K = K
	ans$random.labels = if(is.null(names(K))) {if(nK>0L) paste("varComp", seq_along(K), sep='.') else character(0L)} else names(K)
	if(missing(weights)) ans$weights = NULL else ans$weights = w
	
	if(isTRUE(doFit)) doFit.varComp(ans) else ans
}

varComp.fit = function(Y, X=matrix(0,length(Y),0L), K, control=varComp.control())
{
	{
	#                  a.	varComp: The main fitting function using REML criterion
	#                      i.	Y: response
	#                      ii.	X: fixed effect design matrix. Treated as intercept if missing. A zero matrix can be supplied for REML. 
	#                      iii.	Z: random effect design matrix, or a list of such matrices. Treated as identity if missing. 
	#                      iv.	K: A matrix or a list of matrices, each determining the correlations among observations. This has to be at least p.s.d. in this implementation. 
	#                      v.	information: A character specifying the information matrix to be used during optimization. For doing score tests, this has to be EI; for fitting models, others are OK and the default AOI is usually fast enough.
	#                      vi.	optMethod: A character specifying the optimization method. The default should work well in most cases. 'NRGD' is only a naive implementation and has not been tested extensively. 
	#                      vii.	boundary.eps: A small numeric, below which estimates of tau will be checked for hitting boundary (0). This is a simple attempt at differentiating small but non-zero variance components vs truly zero variance components. Another way of doing this is to divide the K matrices by some large number. 
	#                      viii.	conv: A small numeric, the convergence criterion, only used when optMethod='NRGD'.
	#                      ix.	nStepHalving: A positive integer, giving the max number of step halving during NRGD. 
	#                      x.	nIter: A positive integer, giving the max number of iterations allowed for NRGD.
	#                      xi.	starts: A vector of nonnegative doubles of the same length as the number of variance components. This is the starting value for tau, i.e., the ratio of each variance components to the error variance. 
	#                      xii.	plot.it: A logical scalar, indicating if PREML surface will be plotted (for one variance components only). Please send Long Qu an Email with your data set if you find multiple local maxima, for possible improvements of this function. 
	#                      xiii.	verbose: A logical scalar, indicating if extra information is printed. 
	#                             keepXYK: Logical, indication if original X, Y and K matrices are stored in the results. 
	}
  #if(!is.matrix(X)) X=as.matrix(X)  ## does it matter?
  #if(diff(range(X))>0) X=model.matrix(~X)  ##  dealing with constant X including zero X
  #Y=as.numeric(as.vector(Y)) ## plan to move to varComp formula interface for efficiency
	{
		 optMethod=control$optMethod
		 verbose= control$verbose 
		 starts= control$starts 
		 REML = control$REML 
			if(!isTRUE(REML)) stop("Currently only REML method is implemented")
		 information = control$information 
		 boundary.eps= control$boundary.eps 
		 nlminb.control = control$nlminb 
		 plot.it= control$plot.it 
		 keepXYK= control$keepXYK
		 
		 environment(V)=
		 environment(updateLI)=
		 environment(updateLIkLI)=
		 environment(updateLIy)=
		 environment(PREML)=
		 environment(preprocPREML)=
		 environment(obj)=
		 environment(obj2)=
		 environment(updateNumsPart)=
		 environment(updateNums)=
		 environment(updateDenom)=
		 environment(updateTr1)=
		 environment(updateTr2)=
		 environment(updateNums2)=
		 environment(score)=
		 environment(updateNegGrad)=
		 environment(gradFunc)=
		 environment(OI)=
		 environment(hess)=
		 environment(EI)=
		 environment(AOI)=
		 environment(AEI)=
		 environment(WAI)=
		 environment(updateNegHess)=sys.frame(sys.nframe())
	}

  
  qrx=qr(X)
  # * this block checks if an intercept is in the model or not; 
  # * don't needed this anymore
	 # Q1=qr.Q(qrx)
	 # if(qrx$rank>0 && max(abs(Q1%*%crossprod(Q1, rep(1,length(Y))) - rep(1,length(Y)))) > sqrt(.Machine$double.eps)){  ## no intercept
	   # X=cbind(`(Intercept)`=1, X)
	   # qrx=qr(X)
	 # }
  
  Q2=if(qrx$rank>0) qr.Q(qrx, complete=TRUE)[, -seq_len(qrx$rank),drop=FALSE] else diag(1, length(Y)) #  qr.Q(qrx, complete=TRUE) == IdMat
  y=crossprod(Q2, Y)
  n=length(y)
  
  if(missing(K) || length(K)==0) {  ## linear model fit
		null.sig2=drop(crossprod(y))  / n   
		null.preml=.5*(
		  -n*log(crossprod(y)) -n -n*log(2*pi)-n*log(n)
		)

	  ans=list(
	    ## varComp.fit specific block
		parms=numeric(0L),
		gradients=numeric(0L), 
		hessian=matrix(numeric(0L), 0L, 0L),
		sigma2=drop(null.sig2), 
		varComps=numeric(0),
		n.iter=0L, PREML=drop(null.preml),
		X.Q2=Q2, 
		residual.contrast=y, 
		working.cor=vector('list', 0L),
		
		## varComp common block
		# na.action=NULL,
		# offset = NULL,
		# contrasts=NULL,
		# xzlevels = NULL,
		# terms = NULL, 
		call=match.call(), 
		nobs = length(Y), 
		control=control, 
		random.labels = character(0L), 
		doFit= TRUE,
		
		# frame=if(isTRUE(keepXYK)) NULL else parent.frame(), 
		X=if(isTRUE(keepXYK)) X else NULL,
		# qrx = if(isTRUE(keepXYK)) qrx else NULL, 
		Y=if(isTRUE(keepXYK)) Y else NULL, 
		K=if(isTRUE(keepXYK)) vector('list', 0L) else NULL
	  )
	  class(ans)='varComp'	
	  return(ans)
  }
  
  nK=length(K)
  if(nK>1) plot.it=FALSE
  
  k=vector('list',nK)
  for(j in 1:nK) k[[j]]=crossprod(Q2, K[[j]]%*%Q2)
  
  LIkLI=k; LIy=numeric(n); numsPart=matrix(NA_real_, n, nK); 
  tr2=negHess=nums2=matrix(NA_real_, nK, nK)
  tr1=nums=negGrad=numeric(nK); denom=numeric(1L)
  LI=matrix(NA_real_, n, n)
  diag.1.n=diag(1,n)
  if(nK==1) {
	eigK=eigen(k[[1]],TRUE)
	eigK$tvector=t(eigK$vector)
  }
  
  infoFunc=get(information)
  preprocScore=expression({
	updateLIkLI()
	updateNumsPart();   
	updateNums(); updateDenom(); updateTr1()
  })
  preprocInfo=switch(information,
	OI=expression({ updateNums2(); updateTr2() }),
	EI=expression({ updateTr2();  }),
	AOI=expression({updateNums2();    }),
	AEI=expression({updateNums2();   }),
	WAI=expression({updateNums2();   })
  )
  
  
  last.tau=tau=if(is.null(starts)) {
	minque(y, k, rep(0,nK), lower.bound=.Machine$double.eps^.5, restricted=TRUE)
  }else rep(starts, length=nK)
  ltau=log(max(.Machine$double.eps, tau))
  
  if(optMethod=='nlminb'){
	objNeg=function(tau)-obj(tau)
	gradNeg=function(tau)-gradFunc(tau)
	hessNeg=function(tau)-hess(tau)
	nlminb.fit=nlminb(tau, objNeg, gradNeg, hessNeg, lower=rep(0,nK), control=nlminb.control)
	tau=nlminb.fit$par
	n.nr=nlminb.fit$iterations
  }else if(optMethod=='optim'){
	optim.fit=optim(tau, obj, gradFunc, method='L-BFGS-B', lower=rep(0, nK), control=list(fnscale=-1))
	tau=optim.fit$par
	n.nr=optim.fit$counts
  }else if(optMethod=='NRGD'){
	nStepHalving = 20L
	preprocPREML(tau)
	last.func = PREML()
	n.nr=0L
	repeat{
	  if(verbose) cat('PREML=',last.func, 'tau=',tau,'\n')
	  eval(preprocScore)
	  eval(preprocInfo)
	  updateNegGrad()
	  updateNegHess()
	  
	  negGrad=tau*negGrad
	  negHess=negHess*outer(tau,tau)
	  diag(negHess)=diag(negHess)+negGrad
	  
	  adj.ltau=solve(negHess, negGrad)
	  
	  doNewton=TRUE
	  repeat{
		step.size=1
		n.step=0L
		ltau.new=ltau
		repeat{
		  ltau.new=ltau - step.size * adj.ltau
		  preprocPREML(exp(ltau.new))
		  this.func=PREML()
		  if( ( (this.func>last.func) & doNewton) | ((!doNewton) & (this.func>=last.func)) ) {
			doNewton=TRUE
			break
		  }
		  step.size=step.size/2
		  n.step=n.step+1L
		  if((n.step > nStepHalving) & doNewton ) {  ## try gradient descent
			#warning('Step halving failed')
			adj.ltau=negGrad
			doNewton=!doNewton
			break
		  }
		}
		if(isTRUE(doNewton)) break
	  }
	  
	  last.tau=tau
	  ltau=ltau.new
	  tau=exp(ltau)
	  
	  if(max(abs(last.tau-tau)) / max(abs(last.tau)+abs(tau)) < control$nlminb$x.tol){
		break
	  }
	  n.nr=n.nr+1L
	  last.func=this.func
	  if(n.nr >= control$nlminb$iter.max){
		warning(paste(control$nlminb$iter.max, 'iterations reached'))
		break
	  }
	}    
  }else stop('Method not implemented')
  

  bd.idx=which(tau<boundary.eps)
  if(FALSE){
	  if(length(bd.idx)>0L){#browser()
		preprocPREML(tau)
		eval(preprocScore)
		all.score=score()
		bd.idx=bd.idx[ all.score[bd.idx] < -boundary.eps ]
		
		if(length(bd.idx)==length(tau)) {
		  tau=rep(0,length(tau))
		}else if(length(bd.idx)>0L) {
		  bd.idx0= bd.idx [ all.score[bd.idx]==0 ]
		  if(length(bd.idx0) != length(bd.idx)){
			new.fit=Recall(y, rep(0,length(y)), , k[-bd.idx], information, method, boundary.eps, conv, nStepHalving, control$nlminb$iter.max, starts=tau[-bd.idx], plot.it, verbose)
			if(new.fit$PREML >= PREML()){
			  tau[-bd.idx]=new.fit$parms
			  tau[bd.idx]=0
			}
		  }
		}
	  }
  }

  if((nearZero = length(bd.idx))>0L){  # check boundary
	bd.idx.all= if(nearZero == 1) list(matrix(bd.idx)) else lapply(seq_len(nearZero), combn, x=bd.idx)
	cur.obj=obj(tau)
	tau.bak = tau
	for(i.0 in bd.idx.all){
		if(nrow(i.0)<nK){
			for(j.0 in seq_len(ncol(i.0))){
				this.0 = i.0[,j.0]
				setZero=function(tau){tau0=tau.bak; tau0[this.0]=0; tau0[-this.0]=tau; tau0}  # this.0
				objNeg0=function(tau) objNeg(setZero(tau))
				gradNeg0=function(tau) gradNeg(setZero(tau))[-this.0]
				hessNeg0=function(tau) hessNeg(setZero(tau))[-this.0, -this.0, drop=FALSE]
				nlminb.fit=nlminb(tau[-this.0], objNeg0, gradNeg0, hessNeg0, lower=rep(0,nK-length(this.0)), control=nlminb.control)
				if( -nlminb.fit$objective > cur.obj) {tau = setZero(nlminb.fit$par); cur.obj=-nlminb.fit$objective}
			}
		}else{
			if( (tmp = obj(rep(0, nK))) > cur.obj){ tau = rep(0, nK); cur.obj = tmp}
		}
	}
  }

  preprocPREML(tau)  
  eval(preprocScore)
  eval(preprocInfo)
  updateNegGrad()
  updateNegHess()
  sigma2=crossprod(LIy)/n
	
  if(isTRUE(plot.it)){
	taus=exp(seq(log(1e-9), log(2*tau), length=500))
	taus=sort(unique(c(taus, seq(0, 2*tau, length=500))))
	objs=sapply(taus, obj)
	# x11()
	plot(taus, objs, xlab='tau', ylab='objective', type='o',main='Profiled residual log likelihood')
	abline(v=tau)
  }
  
  nm=names(K)
  if(is.null(nm)) nm = if(nK>0L) paste('varComp', seq_along(K), sep='.') else character(0L)
  ans=list(
	## varComp.fit specific block
	parms=structure(tau, names=nm),
	gradients=structure(-negGrad, names=nm), 
	hessian=structure(hess(tau), dimnames=list(nm,nm)),
	sigma2=drop(sigma2), 
	varComps=structure(drop(tau*sigma2), names=nm),
	n.iter=n.nr, PREML=PREML(),
	X.Q2=Q2, 
	residual.contrast=y, 
	working.cor=structure(k, names=nm),
	
	## varComp common block
	# na.action=NULL,
	# offset = NULL,
	# contrasts=NULL,
	# xzlevels = NULL,
	# terms = NULL, 
	call=match.call(), 
	nobs = length(Y), 
	control=control, 
	random.labels = nm, 
	doFit= TRUE,
	# frame=if(isTRUE(keepXYK)) NULL else parent.frame(), 
	
	X=if(isTRUE(keepXYK)) X else NULL,
	# qrx = if(isTRUE(keepXYK)) qrx else NULL, 
	Y=if(isTRUE(keepXYK)) Y else NULL, 
	K=if(isTRUE(keepXYK)) K else NULL
  )
  class(ans)='varComp'
  ans
  
}

doFit.varComp=function(object)
{
	if(isTRUE(object$doFit)) return(object)
	ans=eval(object$doFit)

	if(ncol(ans$X)>0L){
		this.V=vcov(ans, what='Y')
		this.Vbet=ginv(crossprod(ans$X, solve(this.V,ans$X)))
		this.bet=drop(this.Vbet%*%crossprod(ans$X, solve(this.V, ans$Y)))
		names(this.bet)=colnames(ans$X)
		ans$fixef = this.bet
	}else ans$fixef = numeric(0L)
	ans$na.action = object$na.action
	ans$offset = ans$offset
	ans$contrasts = object$contrasts
	ans$xzlevels = object$xzlevels
	ans$call = object$call
	ans$terms = object$terms
	ans$model = object$model
	ans$X = object$X
	ans$Y = object$Y
	ans$K = object$K
	ans$random.labels = object$random.labels
	ans$weights = object$weights
	ans	
}

if(FALSE){  ## test
	library(SPA3G)
	library(compiler)
	library(quadprog)
	library(nlme)
	library(Matrix)
	library(MASS)
	library(CompQuadForm)
	library(RLRsim)
	library(mvtnorm)

	source("varComp.R")
	source("minque.R")
	source("kernels.R")
	source("coef.varComp.R")
	source("logLik.varComp.R")
	source("fixef.varComp.R")
	source("satterth.R")
	source("model.matrix.varComp.R")
	source("vcov.varComp.R")
	source("print.varComp.R")
	source("formulas.R")
	source("varComp-internal.R")

	mod=distance~age+Sex+age:Sex
	rdm=~Subject

	Z=model.matrix(~Subject, data=Orthodont)

	tmp=varComp(mod, Orthodont, rdm)
	lmef=lme(mod, Orthodont, ~1|Subject)
	VarCorr(lmef)

	# system.time(for(i in seq_len(1e2L)) tmp=varComp(mod, Orthodont, rdm))


	Z=model.matrix(~Subject, data=Orthodont)
	Z0=model.matrix(~0+Subject, data=Orthodont)

	Orthodont$Z=Z
	Orthodont$Z0=Z0

	undebug(varComp)
	tmp2=varComp(mod, Orthodont, random=~Z, normalizeTrace=FALSE)
	tmp2$varComp

	varComp(mod, Orthodont, rdm, normalizeTrace=TRUE)$varComp
	varComp(mod, Orthodont, random=~Z, normalizeTrace=TRUE)$varComp

	varComp(mod, Orthodont, random=~lin0(Z), normalizeTrace=TRUE)$varComp
	varComp(mod, Orthodont, random=~lin0(Z0), normalizeTrace=TRUE)$varComp
	varComp(mod, Orthodont, random=~ibs(Z0), )$varComp
	varComp(mod, Orthodont, random=~quad1(Z0), )$varComp
	varComp(mod, Orthodont, random=~am(Z0), )$varComp

	varComp(mod, Orthodont, random=~ibs(Z0)+lin0(Z0) )$varComp
	varComp(mod, Orthodont, random=~Z0+lin0(Z0) )$varComp
	varComp(mod, Orthodont, random=~Z+lin0(Z0) )$varComp
	varComp(mod, Orthodont, random=~Z+lin0(Z0) , normalizeTrace=TRUE)$varComp

	varComp(mod, Orthodont, random=~Z:lin0(Z0) , normalizeTrace=TRUE)$varComp

	tmp1=varComp(mod, Orthodont, random=~ibs(Z0)*lin0(Z0) , normalizeTrace=TRUE)
	tmp2=varComp(mod, Orthodont, random=~ibs(Z0)+lin0(Z0)+intxn2(ibs(Z0), lin0(Z0)) , normalizeTrace=TRUE)
	tmp1$varComps 
# [1] 51.5910016  0.6939282  0.6939280
	tmp2$varComps
# [1] 51.5921366  0.6939071  0.6939070
	range(tmp1$K[[1]]-tmp2$K[[1]])
# [1] 0 0
	range(tmp1$K[[2]]-tmp2$K[[2]])
# [1] 0 0
	range(tmp1$K[[3]]-tmp2$K[[3]])
# [1] -3.330669e-16  0.000000e+00

	varComp(mod, Orthodont, random=~ibs(Z0)+Sex:ibs(Z0) , normalizeTrace=TRUE)$varComp
	varComp(mod, Orthodont, random=~Sex:ibs(Z0) , normalizeTrace=TRUE)$varComp
	varComp(mod, Orthodont, random=~ibs(Z0):Sex , normalizeTrace=TRUE)$varComp
	(tmp3=varComp(mod, Orthodont, random=~Sex*ibs(Z0) , normalizeTrace=TRUE))$varComp
	(tmp4=varComp(mod, Orthodont, varcov=list(IBS(Z0), tcrossprod(model.matrix(~Sex, Orthodont)[,2L,drop=FALSE])*IBS(Z0)) , normalizeTrace=TRUE))$varComp
	range(model.matrix(tmp3, 'K')[[1]]-model.matrix(tmp4, 'K')[[1]])
	range(model.matrix(tmp3, 'K')[[2]]-model.matrix(tmp4, 'K')[[2]])



}