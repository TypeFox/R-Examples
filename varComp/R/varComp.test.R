varComp.test = function(object, ...) UseMethod("varComp.test", object)

varComp.test.formula = function(object, data, random1, varcov1, random2, varcov2, fit.control, test.control, ...)
{
	this.call=match.call()
	this.call[[1L]] = as.name('varComp')
	m = match(c('random1', 'varcov1', 'random2', 'varcov2','fit.control','object'), names(this.call), 0L)
	if(m[1L]>0L) names(this.call)[m[1L]] = 'random'
	if(m[2L]>0L) names(this.call)[m[2L]] = 'varcov'
	this.call$test.control=this.call$test=NULL
	if(m[5L]>0L) names(this.call)[m[5L]] = 'control'
	if(m[6L]>0L) names(this.call)[m[6L]] = 'fixed'
	this.call$random2 = this.call$varcov2 = NULL
	this.call$doFit = FALSE
	nul=eval.parent(this.call)
	
	this.call=match.call()
	this.call[[1L]] = as.name('varComp')
	m = match( c('random1', 'varcov1', 'random2', 'varcov2','fit.control','object'), names(this.call), 0L)
	if(m[3L]>0L) names(this.call)[m[3L]] = 'random'
	if(m[4L]>0L) names(this.call)[m[4L]] = 'varcov'
	this.call$test.control=this.call$test=NULL
	if(m[5L]>0L) names(this.call)[m[5L]] = 'control'
	if(m[6L]>0L) names(this.call)[m[6L]] = 'fixed'
	this.call$random1 = this.call$varcov1 = NULL
	this.call$doFit = FALSE
	alt=eval.parent(this.call)
	
	if(missing(test.control)) varComp.test(nul, alt, ...) else varComp.test(nul, alt, control=test.control, ...)
}  

varComp.test.varComp = function(object, object2, 
	additional.varcov, null, test='LinScore', 
  control=varCompTest.control(test),  ...)
{
  # information=match.arg(information, informationTypes)
  test=match.arg(test, varCompTests, several.ok=TRUE)
  test[test=='HP01']='SS95'; test=unique(test)
  ntest=length(test)
  if(ntest == 0L) stop(sprintf("Test needs to be given among [%s].", paste0(varCompTests, collapse=',')))

	if(!missing(object2)){  ## two object comparison
		if(!missing(additional.varcov) || !missing(null)) warning("'additional.varcov' and/or 'null' will be ignored when 'object2' is provided.")
		
		if(length(object2$random.labels) > length(object$random.labels)){
			nul=object
			alt=object2
		}else{
			nul=object2
			alt=object
		}
		if(!all(nul$random.labels %in% alt$random.labels)) stop("Two model ares not nested.")
		varComp.test.2modelDoTest(null.fit=nul, alt.fit=alt, test=test, control=control, ...)
	}else if(!missing(additional.varcov)){  ## null fit is available
		if(!missing(null)) warning("'null' is ignored when 'additional.varcov' is provided.")
		varComp.test.nulDoTest(null.fit=object, additional.varcov = additional.varcov, test=test, control=control, ...)
	}else if(!missing(null)){ ##  alt fit is available
		varComp.test.altDoTest(alt.fit=object, null=null, test=test,  control=control, ...)
	}else {  ## default testing all components of object
		varComp.test.altDoTest(alt.fit=object, null=integer(0L), test=test,  control=control, ...)
	}
}
varComp.test.2modelDoTest = function(null.fit, alt.fit, test='LinScore', control=varCompTest.control(test), ...)
{
  # information=match.arg(information, informationTypes)
  test=match.arg(test, varCompTests, several.ok=TRUE)
  test[test=='HP01']='SS95'; test=unique(test)
  ntest=length(test)
  if(ntest == 0L) stop(sprintf("Test needs to be given among [%s].", paste0(varCompTests, collapse=',')))
  
	if(!isTRUE(null.fit$doFit)) {
		if(any(test%in%c(varCompScoreTests, 'RLRT'))){ ## needs a full fitting of null model
			null.fit=doFit.varComp(null.fit)
		}else { ## needs a partial fitting only
			### FIXME: consider Wald and RWD88 separately? 
			bak.control = tmp.control = null.fit$control
			tmp.control$nlminb$iter.max=0L
			null.fit$control = tmp.control
			null.fit=doFit.varComp(null.fit)
			null.fit$control=bak.control
		}
	}
	if((!isTRUE(alt.fit$doFit)) && any(test %in%c('RLRT','Wald'))){	## needs a full alternative fit
		alt.fit = doFit.varComp(alt.fit) 
	}else{ ## needs a partial alternative fit
		tmp.control = bak.control = alt.fit$control
		tmp.control$nlminb$iter.max=0L
		alt.fit$control = tmp.control
		alt.fit = doFit.varComp(alt.fit)
		alt.fit$control = bak.control
	}
  if(!all(null.fit$random.labels %in% alt.fit$random.labels)) stop("Two models are not nested.")
  null=which(alt.fit$random.labels %in% null.fit$random.labels)
  nNull = length(null)
  nK = length(alt.fit$working.cor) 
  
  tau.idx=seq_len(nK)
  if(nNull>0L) tau.idx=tau.idx[-null]
  n.tau=length(tau.idx)	
  
	environment(varComp.test.Common)  = sys.frame(sys.nframe())
	varComp.test.Common()
}
varComp.test.nulDoTest = function(null.fit, additional.varcov, test='LinScore', control=varCompTest.control(test), alt.fit=NULL, ...)
{
  # information=match.arg(information, informationTypes)
  test=match.arg(test, varCompTests, several.ok=TRUE)
  test[test=='HP01']='SS95'; test=unique(test)
  ntest=length(test)
  if(ntest == 0L) stop(sprintf("Test needs to be given among [%s].", paste0(varCompTests, collapse=',')))

	if(!isTRUE(null.fit$doFit)) {
		if(any(test%in%c(varCompScoreTests, 'RLRT'))){ ## needs a full fitting of null model
			null.fit=doFit.varComp(null.fit)
		}else { ## needs a partial fitting only
			### FIXME: consider Wald and RWD88 separately? 
			bak.control = tmp.control = null.fit$control
			tmp.control$nlminb$iter.max=0L
			null.fit$control = tmp.control
			null.fit=doFit.varComp(null.fit)
			null.fit$control=bak.control
		}
	}
	
  null=seq_along(null.fit$working.cor)
  nNull = length(null)
  if(is.matrix(additional.varcov)) additional.varcov = list(additional.varcov)
  nK = length(additional.varcov) + nNull

  tau.idx=seq_len(nK)
  if(nNull>0L) tau.idx=tau.idx[-null]
  n.tau=length(tau.idx)
  
	if(!isTRUE(alt.fit$doFit)){
		if(is.null(alt.fit)){
			k = vector('list', nK)
			k[seq_len(nNull)] = null.fit$working.cor
			Q2 = null.fit$X.Q2
			for(i.k in (nNull+1L):nK) k[[i.k]] = crossprod(Q2, additional.varcov[[i.k-nNull]]%*%Q2)
		}
		if(any(test %in%c('RLRT','Wald'))){	## needs a full alternative fit
			if(!is.null(alt.fit)) {
				alt.fit = doFit.varComp(alt.fit) 
			}else {
				alt.fit = varComp.fit(Y=null.fit$residual.contrast, K = k, control = null.fit$control)
			}
		}else if(any(test%in%varCompScoreTests)){ ## needs a partial fit
			if(!is.null(alt.fit)){
				tmp.control=bak.control=alt.fit$control
				tmp.control$nlminb$iter.max=0L
				alt.fit$control=tmp.control
				alt.fit=doFit.varComp(alt.fit)
				alt.fit$control=bak.control
			}else{
				tmp.control=null.fit$control
				tmp.control$nlminb$iter.max=0L
				alt.fit=varComp.fit(Y=null.fit$residual.contrast, K=k, control=tmp.control)
			}
		}else stop("should not reach here")
	}
	environment(varComp.test.Common)  = sys.frame(sys.nframe())
	varComp.test.Common()
}

varComp.test.altDoTest = function(alt.fit, null=integer(0L), test='LinScore', control=varCompTest.control(test), null.fit=NULL, ...)
{
  # information=match.arg(information, informationTypes)
  test=match.arg(test, varCompTests, several.ok=TRUE)
  test[test=='HP01']='SS95'; test=unique(test)
  ntest=length(test)
  if(ntest == 0L) stop(sprintf("Test needs to be given among [%s].", paste0(varCompTests, collapse=',')))
  
  null=as.integer(round(null))
  nNull=length(null)
	
    if( !isTRUE(alt.fit$doFit) ) {
		if( any(test %in% c('RLRT','RWD88','Wald')) ) { ## needs a full fit of alternative
			alt.fit = doFit.varComp(alt.fit)
		}else if( any(test %in%varCompScoreTests)){ ## needs a quick fit of alternative
			tmp.control=bak.control=alt.fit$control
			tmp.control$nlminb$iter.max=0L
			tmp.call=alt.fit$doFit
			tmp.call[['control']] = tmp.control
			alt.fit=eval(tmp.call)
			alt.fit$control = bak.control
		}
	}

	nK = length(alt.fit$working.cor)
	if(any(null<1L | null>nK) | nNull>=nK) {
		stop('The number of null components should be less than the number of non-error variance components')		
	}
	tau.idx=seq_len(nK)
	if(nNull>0L) tau.idx=tau.idx[-null]
	n.tau=length(tau.idx)
  
	if( !isTRUE(null.fit$doFit)) {
		if(any(test%in%c(varCompScoreTests,'RLRT') ) ) { ## full null fit
			if(is.null(null.fit)) {
				null.fit = varComp.fit(Y=alt.fit$residual.contrast, K=alt.fit$working.cor[null], control=alt.fit$control)
			}else null.fit=doFit.varComp(null.fit)
		}else{
			null.fit=NULL
		}
	}

	environment(varComp.test.Common) = sys.frame(sys.nframe())
	varComp.test.Common()
}

varComp.test.Common = evalq(function()  ## not to be called directly! 
{	
	# This function requires external objects: 
	#	"test", "varCompScoreTests", "control", "alt.fit",  "null.fit", 
	#	"null", "updateNegGrad",  "updateNegHess", "updateTr2", "tau.idx"

	call.lists=list(LinScore = c('infoMat', 'all.scores','tr1','n','LIkLI','non.pd'),
				 SS95 = c('infoMat', 'all.scores', 'tr1','n','LIkLI','LIy'),
				 RLRT = c('alt.fit', 'null.fit'),
				 VM03 = c('infoMat', 'all.scores', 'tr1','n','LIkLI','LIy'))
	#environment(char2list) = sys.frame(sys.nframe())
	
	if(any(test%in%varCompScoreTests)){
		{
			 environment(V)=
			 environment(updateLI)=
			 # environment(updateLIk)=
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
			infoFunc=get(control$information)
			preprocScore=expression({
				updateLIkLI()
				updateNumsPart();   
				updateNums(); updateDenom(); updateTr1()
			})
			preprocInfo=switch(control$information,
				OI=expression({ updateNums2(); updateTr2() }),
				EI=expression({ updateTr2();  }),
				AOI=expression({updateNums2();    }),
				AEI=expression({updateNums2();   }),
				WAI=expression({updateNums2();   })
			)			 
		}

		y=alt.fit$residual.contrast
		n=length(y); diag.1.n = diag(1, n)
		k=alt.fit$working.cor
		nK=length(k)
		
		LIkLI=k; LIy=numeric(n); numsPart=matrix(NA_real_, n, nK); 
		tr2=negHess=nums2=matrix(NA_real_, nK, nK)
		tr1=nums=negGrad=numeric(nK); denom=numeric(1L)
		LI=matrix(NA_real_, n, n)
		if(nK==1) {
			eigK=eigen(k[[1]],TRUE)
			eigK$tvector=t(eigK$vector)
		}
		
		tau=numeric(nK); tau[null]=null.fit$parms
		preprocPREML(tau)
		eval(preprocScore)
		eval(preprocInfo)
		updateNegGrad()
		updateNegHess()
		non.pd=FALSE
		if(control$information!='EI' && (min(eigen(negHess,TRUE,TRUE)$value)< -1e-6 || any(diag(as.matrix(negHess))<0) ) ){
			non.pd=TRUE		
			updateTr2();
			negHess=EI()
		}
		infoMat = as.matrix(nearPD(negHess)$mat)    ## added line on 062212
		sigma2=null.fit$sigma2   #sigma2=crossprod(LIy)/n
		all.scores= - negGrad    #score()
	}else if("RLRT"%in%test){
	}else stop("should not reach this line")
  
  test.list=vector('list', length(test)); names(test.list)=test
  for(i.test in seq_along(test)){
    # this.control=control[[test[i.test]]] ## do NOT move this into the call below, because the name 'control' has conflict.
	# this.fun=paste("varComp", test[i.test], "test", sep='.')
	# thisToCall=call(this.fun, control=this.control
          # ,null=null, null.fit=null.fit, alt.fit=alt.fit, X=X, Y=Y, K=K, k=k, lin.form=lin.form, LIkLI=LIkLI, LIy=LIy, tr1=tr1
          # ,infoMat=negHess, negHess=negHess, all.scores=all.scores, n=n, tau.idx=tau.idx, non.pd=non.pd
          # ,obs.score=obs.score)
		  
    # test.list[[i.test]]=eval(thisToCall)
	test.list[[i.test]] = do.call(
		what = paste('varComp', test[i.test], 'test', sep='.'),
		args = c(mget(call.lists[[test[i.test]]], inherits=TRUE), list(tau.idx = tau.idx, control=control[[test[i.test]]]))
    )
  }
  class(test.list)='varComp.test'
  test.list
}, refugeEnvironment)


if(FALSE){
	library(SPA3G)
	library(compiler)
	library(quadprog)
	library(nlme)
	library(Matrix)
	library(MASS)
	library(CompQuadForm)
	library(RLRsim)
	library(mvtnorm)

	getwd("C:/Users/Longor/Dropbox/varComp/R")
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
	source("varComp.test.R")
	source("varComp.R")
	source('varCompTest.control.R')
	source("varComp.LinScore.test.R")
	source("varComp.LinScore.SSAS155.R")
	source("varComp.LinScore.Satterthwaite.R")
	source("varComp.LinScore.Normal.R")
	source("varComp.RLRT.test.R")
	source("varComp.SS95.test.R")
	source("varComp.VM03.test.R")
	source("get.seed.R")
	source("pchibarsq.R")

	mod=distance~age+Sex+age:Sex
	rdm=~Subject

	Z=model.matrix(~Subject, data=Orthodont)

	(tmp=varComp(mod, Orthodont, rdm))
	lmef=lme(mod, Orthodont, ~1|Subject)
	VarCorr(lmef)
	tmp0=varComp(mod, Orthodont, rdm, doFit=FALSE)
	tmp2 = doFit.varComp(tmp)
	identical(tmp,tmp2)		## should be true


	# system.time(for(i in seq_len(1e2L)) tmp=varComp(mod, Orthodont, rdm))


	Z=model.matrix(~Subject, data=Orthodont)
	Z0=model.matrix(~0+Subject, data=Orthodont)

	Orthodont$Z=Z
	Orthodont$Z0=Z0

	debug(varComp.test.2modelDoTest)
	varComp.test(mod, Orthodont, random1=~ibs(Z0), random2=~Sex:ibs(Z0))
	
	
	mod0=distance~Sex
	nulFit = varComp(mod0, Orthodont, random = ~ibs(Z0))
	altFit = varComp(mod0, Orthodont, random = ~Sex:ibs(Z0))
	debug(varComp.test); debug(varComp.test.nulDoTest)
	varComp.test(nulFit, additional.varcov=model.matrix(altFit, 'K')[[2]])

	debug(varComp.test); debug(varComp.test.altDoTest)
	varComp.test(altFit, null=1L)

	varComp.test(altFit, null=1L, control=varCompTest.control(test='LinScore',LinScore.method=c('AS155','SSAS155')))
	
	varComp.test(altFit, null=1L, control=varCompTest.control(test='LinScore',LinScore.method=c('AS155','Satterth')))

	varComp.test(altFit, null=1L, control=varCompTest.control(test='LinScore',LinScore.method=c('AS155','Normal')))

	varComp.test(altFit)
	
	varComp.test(nulFit, test='RLRT')
	varComp.test(altFit, null=1L, test='RLRT')
	
	varComp.test(nulFit, test='SS95')
	varComp.test(altFit, null=1L, test='SS95')
	varComp.test(altFit, test='SS95')
	
	varComp.test(nulFit, test='SS95', control=varCompTest.control(test='SS95', SS95.method=c('ChiBarSq', 'ChiBarSq')))
	varComp.test(altFit, null=1L, test='SS95', control=varCompTest.control(test='SS95', SS95.method=c('ChiBarSq', 'ChiBarSq')))
	varComp.test(altFit, test='SS95', control=varCompTest.control(test='SS95', SS95.method=c('ChiBarSq', 'ChiBarSq')))
	
	varComp.test(nulFit, test='VM03')
	varComp.test(altFit, null=1L, test='VM03')
	varComp.test(altFit, test='VM03')
	
	varComp.test(nulFit, test='VM03', control=varCompTest.control(test='SS95', SS95.method=c('ChiBarSq', 'ChiBarSq')))
	varComp.test(altFit, null=1L, test='VM03', control=varCompTest.control(test='SS95', SS95.method=c('ChiBarSq', 'ChiBarSq')))
	varComp.test(altFit, test='VM03', control=varCompTest.control(test='SS95', SS95.method=c('ChiBarSq', 'ChiBarSq')))
	
}
