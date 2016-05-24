#
#Code to fit some of the Spatial Ecconometrics models for a given rho
#


#Simultaneous Autorregresive Model
#
#Return an INLA model
#Formula is for the FIXED effects ONLY
#'...' are passed to the INLA call
#areaid is the column to be used as area index in the random spatial effect
#

sem.inla<-function(formula, d, W, rho, improve=TRUE, impacts=FALSE, fhyper=NULL, probit=FALSE,...)
{
	#require(INLA)

	IrhoW<-Matrix::Diagonal(nrow(W))-rho*W
	#IrhoW2<-t(IrhoW)%*%IrhoW
	IrhoW2<-Matrix::crossprod(IrhoW)

	#environment(formula)<-environment()
	#This is a fix to be able to use improve=TRUE later
        #environment(formula)<-environment()
        assign("IrhoW2", IrhoW2, environment(formula) )


	if(is.null(fhyper))
	{
	formula<-update(formula, 
	  . ~ . + f(idx, model="generic0", Cmatrix =IrhoW2) )
	}
	else
	{
	formula<-update(formula, 
	  . ~ . + f(idx, model="generic0", Cmatrix =IrhoW2, hyper=fhyper) )
	}


	res<-INLA::inla(formula, data=d, ...)

	if(improve)
		res<-INLA::inla.rerun(res)#inla.hyperpar(res, diff.logdens=20)

	#Compute log-determinat to correct the marginal-loglikelihood
	res$logdet<-as.numeric(Matrix::determinant(IrhoW2)$modulus)
	res$mlik<-res$mlik+res$logdet/2

	#Add impacts as the marginals of the covariate coefficients
	res$impacts<-FALSE
	if(impacts)
	{
		res$impacts<-TRUE
		#Add summary.
		res$impacts$summary.impacts<-res$summary.fixed[rownames(res$summary.fixed) != "(Intercept)",]

		#Add marginals
		res$impacts$marginals.impacts<-res$marginals.fixed[ names(res$marginals.fixed) !="(Intercept)"]

		if(probit)#Re-scale marginals
		{
			Dfm<-mean(dnorm(res$summary.linear.predictor[,1]))

			mm<-lapply(res$impacts$marginals.impacts,
			function(X){
	INLA::inla.tmarginal(function(x){Dfm*x}, X)
			})

	names(mm)<-names(res$impacts$marginals.impacts)
		res$impacts$marginals.impacts<-mm

		res$impacts$summary.impacts<-Dfm*res$impacts$summary.impacts

		}

		
	}

	return(res)
}



#Spatial lag model
#mm is model.matrix to speed up computations
#
#using family="binomial" with a probit link then we will be fitting a
#spatial probit model 
#
slm.inla<-function(formula, d, W, rho, mmatrix=NULL, improve=TRUE, 
   impacts=FALSE, fhyper=NULL, probit=FALSE,...)
{
	#require(INLA)

	IrhoW <- Matrix::Diagonal(nrow(W))-rho*W
	#IrhoW2<-t(IrhoW)%*%IrhoW
	IrhoW2<-Matrix::crossprod(IrhoW)

        #environment(formula)<-environment()
        #This is a fix to be able to use improve=TRUE later
        #environment(formula)<-environment()
        assign("IrhoW2", IrhoW2, environment(formula) )


	if(is.null(mmatrix))
		mmatrix<-model.matrix(formula, d)

	mm<-as.data.frame(as.matrix(solve(IrhoW)%*%mmatrix))
	names(mm)<-paste("x", 1:ncol(mm), sep="")
	xnam<-names(mm)

	#mm$z<-1:nrow(mm)

	d2<-cbind(d, mm)

	fmla <- paste(as.character(formula)[2], "~ -1+", paste(xnam, collapse= "+"))
	if(is.null(fhyper))
	fmla<-paste(fmla, "+f(idx, model=\"generic0\", Cmatrix=IrhoW2)", sep="")
	else
	fmla<-paste(fmla, "+f(idx, model=\"generic0\", Cmatrix=IrhoW2, hyper=fhyper)", sep="")
	fmla<-as.formula(fmla)


	res <- INLA::inla(fmla, data=d2, ...)


	if(improve)
		res <- INLA::inla.rerun(res)#inla.hyperpar(res, diff.logdens=20)

	#Compute log-determinat to correct the marginal-loglikelihood
	res$logdet<-as.numeric(Matrix::determinant(IrhoW2)$modulus)
	res$mlik<-res$mlik+res$logdet/2

	res$impacts<-FALSE
	if(impacts)
	{
		res$impacts<-TRUE
		#Compute weights for impacts
		if(!probit)
		{
		wtotal<-1/(1-rho)
		wdirect<-trIrhoWinv(W, rho)/nrow(W)
		}
		else
		{
		Df<-dnorm(res$summary.linear.predictor[,1])
		wtotal<-mean(Df)*1/(1-rho)
		wdirect<-trIrhoWinv(W, rho, Df = Matrix::Diagonal(x=Df))/nrow(W)
		}
		windirect<-wtotal-wdirect

		#
		#TOTAL IMPACTS
		#
		#Add summary.
		#Total impacts
		res$summary.total.impacts<-wtotal*res$summary.fixed[rownames(res$summary.fixed) != "(Intercept)",1:2]
		#Direct impacts
		res$summary.direct.impacts<-wdirect*res$summary.fixed[rownames(res$summary.fixed) != "(Intercept)",1:2]
		#INDirect impacts
		res$summary.indirect.impacts<-windirect*res$summary.fixed[rownames(res$summary.fixed) != "(Intercept)",1:2]


		#Add marginals
		res$marginals.total.impacts<-res$marginals.fixed[ names(res$marginals.fixed) !="(Intercept)"]
		res$marginals.direct.impacts<-res$marginals.fixed[ names(res$marginals.fixed) !="(Intercept)"]
		res$marginals.indirect.impacts<-res$marginals.fixed[ names(res$marginals.fixed) !="(Intercept)"]

		for(i in 1:length(res$marginals.total.impacts))
		{
	xx<-res$marginals.total.impacts[[i]]
	res$marginals.total.impacts[[i]]<-rescalemarg(xx, wtotal)
	res$marginals.direct.impacts[[i]]<-rescalemarg(xx, wdirect)
	res$marginals.indirect.impacts[[i]]<-rescalemarg(xx, windirect)
		}
	}

	return(res)
}


#Spatial Durbin Model
#intercept, whether the intercept is in the model or not

sdm.inla<-function(formula, d, W, rho, mmatrix=NULL, intercept=TRUE, 
   impacts=FALSE, improve=TRUE, fhyper=NULL, probit=FALSE, ...)
{

	#require(INLA)

	IrhoW <- Matrix::Diagonal(nrow(W))-rho*W
	#IrhoW2<-t(IrhoW)%*%IrhoW
	IrhoW2<-Matrix::crossprod(IrhoW)

        #environment(formula)<-environment()
        #This is a fix to be able to use improve=TRUE later
        #environment(formula)<-environment()
        assign("IrhoW2", IrhoW2, environment(formula) )


	if(is.null(mmatrix))
	{
		mmatrix<-model.matrix(formula, d)

		if(intercept)
			mmW<-W%*%mmatrix[,-1]#Remove intercept
		else
			mmW<-W%*%mmatrix

		mmatrix<-cbind(mmatrix, as.matrix(mmW))
	}
	mm<-as.data.frame(as.matrix(solve(IrhoW)%*% mmatrix))


	names(mm)<-paste("x", 1:ncol(mm), sep="")
	xnam<-names(mm)


	d2<-cbind(d, mm)

	fmla <- paste(as.character(formula)[2], "~ -1+", paste(xnam, collapse= "+"))

	if(is.null(fhyper))
	fmla<-paste(fmla, "+f(idx, model=\"generic0\", Cmatrix=IrhoW2)", sep="")
	else
	fmla<-paste(fmla, "+f(idx, model=\"generic0\", Cmatrix=IrhoW2, hyper=fhyper)", sep="")
	
	fmla<-as.formula(fmla)

	#if(exists(lincomb))
	#	lc.many<-lincomb
	#else
	#	lc.many<-list()
	if(impacts)#Compute impacts
	{

		#Compute weights for impacts: c(beta_k, gamma_k)
		if(!probit)
		{
                wtotal<-rep(1/(1-rho), 2)
                wdirect<-c(trIrhoWinv(W, rho), trIrhoWinv(W, rho, 1))/nrow(W)
		}
		else
		{
		stop("Cannot compute impacts because we need to run the model twice")
		#FIXME: Update the following code
		Df<-dnorm(res$summary.linear.predictor[,1])
                wtotal<-mean(Df)*rep(1/(1-rho), 2)
                wdirect<-c(trIrhoWinv(W, rho, Df = Matrix::Diagonal(x=Df)), trIrhoWinv(W, rho, 1, Df = Matrix::Diagonal(x=Df)))/nrow(W)
		}
                windirect<-wtotal-wdirect

		ncov<-(ncol(mm)-1)/2

		#Linear combinations for total impacts
		lc.total.impacts<-lapply(1:ncov, function(X){
			lc<-list(a=wtotal[1], b=wtotal[2])
			names(lc)<-xnam[1+c(X, X+ncov)]
			INLA::inla.make.lincomb(lc)
		})
		lc.total.impacts<-do.call(c, lc.total.impacts)
		names(lc.total.impacts)<-paste("totalimp.", xnam[1+1:ncov], sep="")


		#Linear combinations for direct impacts
		lc.direct.impacts<-lapply(1:ncov, function(X){
			lc<-list(a=wdirect[1], b=wdirect[2])
			names(lc)<-xnam[1+c(X, X+ncov)]
			INLA::inla.make.lincomb(lc)
		})
		lc.direct.impacts<-do.call(c, lc.direct.impacts)
		names(lc.direct.impacts)<-paste("directimp.", xnam[1+1:ncov], sep="")

		#Linear combinations for indirect impacts
		lc.indirect.impacts<-lapply(1:ncov, function(X){
			lc<-list(a=windirect[1], b=windirect[2])
			names(lc)<-xnam[1+c(X, X+ncov)]
			INLA::inla.make.lincomb(lc)
		})
		lc.indirect.impacts<-do.call(c, lc.indirect.impacts)
		names(lc.indirect.impacts)<-paste("indirectimp.", xnam[1+1:ncov], sep="")


		lc.impacts<-c(lc.total.impacts,lc.direct.impacts, 
			lc.indirect.impacts)
		#lc<-as.list(rep(1/(1-rho), length(xnam)-1))
		#names(lc)<-xnam[-1]
		#lc.impacts = inla.make.lincomb(lc=lc)

		res <- INLA::inla(fmla, data=d2, lincomb=lc.impacts, ...)
	}
	else
	{
		res <- INLA::inla(fmla, data=d2, ...)
	}



	if(improve)
		res <- INLA::inla.rerun(res)#inla.hyperpar(res, diff.logdens=20)

	#Compute log-determinat to correct the marginal-loglikelihood
	res$logdet<-as.numeric(Matrix::determinant(IrhoW2)$modulus)
	res$mlik<-res$mlik+res$logdet/2

	res$impacts<-FALSE
	if(impacts)
	{
		res$impacts<-TRUE

		#Summary statistics
		idx<-substr(rownames(res$summary.lincomb.derived), 1,3)
		res$summary.total.impacts<-res$summary.lincomb.derived[idx=="tot",]
		res$summary.direct.impacts<-res$summary.lincomb.derived[idx=="dir",]
		res$summary.indirect.impacts<-res$summary.lincomb.derived[idx=="ind",]


		#Marginals
		idx<- substr(names(res$marginals.lincomb.derived), 1, 3)
		res$marginals.total.impacts<-res$marginals.lincomb.derived[idx=="tot"]
		res$marginals.direct.impacts<-res$marginals.lincomb.derived[idx=="dir"]
		res$marginals.indirect.impacts<-res$marginals.lincomb.derived[idx=="ind"]

	}

	return(res)

}



#
#Log-prior for rho: logit(rho) ~ N(0, prec=.1)
#
logprrho<-function(rho)
{
	dnorm(log(rho/(1-rho)), 0, sqrt(1/.1), log=TRUE)-log(rho*(1-rho))
}


#Compute trace of (I-\rho*W)^{-1} using properties of trace and 
#Taylor expansion approximation of order 'order'
#
#W: Adjacency matrix
#rho: rho
#offset: Number of times (I-\rho*W)^{-1} is multiplied by W (for sdm model)
#order: order of Taylor expansion
#direct: User direct method, i.e., matrix multiplication, etc.
#Df:Diagonal matrix used to compute the impacts in the Probit model
#      only used if direct=TRUE.
trIrhoWinv<-function(W, rho, offset=0, order=20, direct=TRUE, Df = Matrix::Diagonal(nrow(W)))
{
	if(!direct)
	{
	lambdas<-eigen(W)$values

	tr<-0
	for(i in 0:order)
		tr<-tr+sum(lambdas^(i+offset))*(rho^i)
	}
	else
	{
		tr<-0
		WW <- Matrix::Diagonal(nrow(W))
		if(offset>0)
		{
			for(i in 1:(offset))
				WW<-WW%*%W
		}
		WW=Df%*%WW
		for(i in 0:order)
		{
			tr<-tr+sum(Matrix::diag(WW))*(rho^i)
			WW<-WW%*%W
		}
	}
	return(as.double(tr))
}


#Re-scale marginal to comkpute the ditribution of w*x
#Used when computing the impacts
#
#xx: nx2 matrix with x,y values of marginal
#w: weight
rescalemarg<-function(xx, w)
{
#	xx[,1]<-xx[,1]*w
#	xx[,2]<-xx[,2]/w
	#USe inla.marginal.transform
	xx <- INLA::inla.tmarginal(function(x){w*x}, xx)

	return(xx)
}



#
#Recompute impacts summaries using marginals
#
#obj: An object with impacts results
#
recompute.impacts<-function(obj, impacts=c("total", "direct", "indirect"))
{
	if(!obj$impacts)#Nothing to do...
		return(obj)
	
	for(ii in impacts)
	{
		sumname<-paste("summary.", ii,".impacts", sep="")
		margname<-paste("marginals.", ii,".impacts", sep="")

		stab<-obj[[sumname]]
		for(i in 1:nrow(stab))
		{
	xmean <- INLA::inla.emarginal(function(x){x}, obj[[margname]][[i]])
	xsd<-sqrt(INLA::inla.emarginal(function(x){(x-xmean)^2}, obj[[margname]][[i]]))

		stab[i,1:2]<-c(xmean, xsd)
		}
		obj[[sumname]]<-stab[,1:2]
	}

	return(obj)
}

