
gam.gradients <- function(
    mod,
    phenotype,
    covariates=NULL,
    standardized=FALSE,
    se.method='boot.para',
    n.boot=1000,
    parallel="no",
    ncpus=1,
    refit.smooth=FALSE
) {
	
	## sub-functions for gam.gradients()
	
	Wbar <- function(x, mod, phenotype,covariates) {
    		new.d <- as.data.frame(mod$model[,c(phenotype,covariates)])
    		names(new.d)<-c(phenotype,covariates)
    		new.d2<-new.d
    		for (i in 1:length(x)) {
    		    new.d[,as.character(phenotype[i])] <-
    		    new.d[,as.character(phenotype[i])]+x[i]
    		}
    
    		p <- predict.gam(
			   object=mod,
    	 	   newdata=new.d,
    	 	   newdata.guaranteed=TRUE,
    	 	   type="response"
    		)
    	return(mean(p))
	}

	gradients <- function(m, phenotype, covariates) {
    	nTraits = length(phenotype)
    	first.derivatives <- grad(func=Wbar, x=rep(0, nTraits), 
    	         mod=m, phenotype=phenotype, covariates=covariates)
    	second.derivatives <- hessian(func=Wbar, x=rep(0, nTraits), 
    	         mod=m, phenotype=phenotype, covariates=covariates)
    	denom <- Wbar(x=rep(0,nTraits), mod=m, phenotype=phenotype,
    	         covariates=covariates)
    	
    	beta <- first.derivatives  / denom
    	gamma <- second.derivatives / denom
    	
    	if(standardized){
    		sds<-apply(as.matrix(mod$model[,phenotype]),2,sd)
    		beta<-beta*sds
    		gamma<-gamma*outer(sds,sds)
    	}
    	return( list(
    	    beta  = beta,
    	    gamma = gamma
    	))
	}

	boot.gradients <- function(data, original, mod, 
	                             phenotype,covariates) {
  		d <- as.data.frame(data)[original,]

		mod.prime <- NULL
		if(refit.smooth){
		  mod.prime <- gam(
        		as.formula(mod$formula),
        		data=d,
        		family=mod$family,
        		start=mod$coefficients
    		  )
		}else{
		  mod.prime <- gam(
        		as.formula(mod$formula),
        		data=d,
        		family=mod$family,
        		start=mod$coefficients,
        		sp=mod$sp
    		  )
    		}

    	g <- gradients(mod.prime, phenotype,covariates)

    	return(c(g$beta,diag(g$gamma),
    	                         g$gamma[upper.tri(g$gamma)]))
	}


	simulate.gam.gsg<-function(mod){
		sim.lin.pred<-predict(mod)
		s<-NULL
		if(mod$family[[1]]=='poisson'&mod$family[[2]]=='log') {
			s<-rpois(length(sim.lin.pred),exp(sim.lin.pred))
		}
		if(mod$family[[1]]=='binomial') {
			s<-rbinom(length(sim.lin.pred),1,
			                     inv.logit(sim.lin.pred))
		}
		if(strsplit(mod$family[[1]],split="\\(")[[1]][1]==
		            'Negative Binomial') {
			s<-rnbinom(length(sim.lin.pred),size=mod$family$getTheta(),
			                     mu=exp(sim.lin.pred))
		}
		if(mod$family[[1]]=='gaussian') {
			s<-rnorm(length(sim.lin.pred),
			   sim.lin.pred,sqrt(mod$sig2))
		}		
		s
	}


	boot.gradients.para <- function(data, original, mod, 
	                                      phenotype, covariates) {
			# this function is set up for use with 'boot()',
			# but does not use argument 'original'.  This
			# allows each replicate to be executed in the
			# same way as for case bootstrapping, but with
			# the parametric bootstrap part of the analysis
			# done internally in this sub-function
		
    		d <- as.data.frame(data)
    		d[[attr(attr(mod$terms,'factors'),
    		            'dimnames')[[1]][1]]]<-simulate.gam.gsg(mod)

    		mod.prime <- NULL
    if(strsplit(mod$family[[1]],split="\\(")[[1]][1]==
		            'Negative Binomial'){
		if(refit.smooth){
			theta<-mod$family$getTheta()
		  mod.prime <- gam(
        		as.formula(mod$formula),
        		data=d,
        		family=negbin(theta),
        		start=mod$coefficients
    		  )
		}else{
			theta<-mod$family$getTheta()
		  mod.prime <- gam(
        		as.formula(mod$formula),
        		data=d,
        		family=negbin(theta),
        		start=mod$coefficients,
        		sp=mod$sp
    		  )
    	}
    }else{ #everything other than negbinom
		if(refit.smooth){
		  mod.prime <- gam(
        		as.formula(mod$formula),
        		data=d,
        		family=mod$family,
        		start=mod$coefficients
    		  )
		}else{
		  mod.prime <- gam(
        		as.formula(mod$formula),
        		data=d,
        		family=mod$family,
        		start=mod$coefficients,
        		sp=mod$sp
    		  )
    	}
    }
    	g <- gradients(mod.prime, phenotype, covariates)

    	return(c(g$beta,diag(g$gamma),
    	                     g$gamma[upper.tri(g$gamma)]))
	}


	permute.W.refit<-function(data, original, mod, 
	                              phenotype, covariates){

	    d <- as.data.frame(mod$model[,
	         c(attr(attr(mod$terms,"factors"),"dimnames")[[1]][1],
	         phenotype, covariates)])
	    names(d)<-c(attr(attr(mod$terms,"factors"),"dimnames")[[1]][1],
	         phenotype, covariates)
	    
	    d[ ,attr(attr(mod$terms,"factors"),
	         "dimnames")[[1]][1] ] <- d[ sample(1:(dim(d)[1]),
	         dim(d)[1],replace=FALSE) , attr(attr(mod$terms,"factors"),
	         "dimnames")[[1]][1] ]

		mod.prime <- gam(
        		as.formula(mod$formula),
        		data=d,
        		family=mod$family,
        		start=mod$coefficients
    		  )
    p<-predict(mod.prime,type="response")

    		g <- gradients(mod.prime, phenotype, covariates)
    
    		return(c(g$beta,diag(g$gamma),g$gamma[upper.tri(g$gamma)],var(p)))
	}

	posterior.sim.gradients<-function(m,n.boot, phenotype,covariates){
    	nTraits <- length(phenotype)
    	mu <- m$coefficients
    	Vp <- m$Vp
    	coef.boot <- rmvnorm(n.boot,mu,Vp)
    	grads.boot <- array(dim=c(n.boot,2*nTraits+(nTraits^2-nTraits)/2))
    	for (b in 1:n.boot) {
      		mod.prime <- m
      		m$coefficients <- coef.boot[b,]
      		g <- gradients(mod.prime, phenotype, covariates)
      		grads.boot[b,]<-c(g$beta,diag(g$gamma),g$gamma[upper.tri(g$gamma)])
    	}
    	return(grads.boot)
	}

	## end sub-functions for gam.gradients()	


	## some checks
		
    cls <- class(mod)
    rightType <- any(cls %in% c("gam"))
    if (!rightType) {
        stop("Argument 'mod' must be a generalized additive model.")
    }
    if (!is.character(phenotype)) {
        stop("Argument 'phenotype' must be a character vector of terms from the fitness model.")
    }
    if (!is.logical(standardized)) {
        stop("Argument 'standardized' must be logical.")
    }
    if (!is.character(se.method) ) {
        stop("Argument 'se.method' must be character and one of 'boot.para', 'boot.case, 'posterior', 'permute', or 'n'.")
    }
    if (!is.numeric(n.boot)) {
        stop("Argument 'n.boot' must be numeric.")
    }
    if (n.boot <= 1) {
        stop("Argument 'n.boot' must be greater than one.")
    }
    if (!is.character(parallel)) {
        stop("Argument 'parallel' must be a character vector, see ?boot.")
    }
    if (!is.numeric(ncpus)) {
        stop("Argument 'ncpus' must be numeric.")
    }

    termLabels <- attr(x=terms(mod), which='term.labels')
    nTerms <- length(termLabels)
    if (!all(phenotype %in% termLabels)) {
        stop("Some selected phenotypes not included in the model.")
    }


    nTraits<-length(phenotype)

    # Calculate estimates
    ests <- list()
    g <- gradients(mod, phenotype,covariates)
    ests$estimates<-c(g$beta,diag(g$gamma),g$gamma[upper.tri(g$gamma)])

    # Calculate standard errors
    if (se.method %in% c('posterior','boot.case','boot.para')) {
        boot.res<-NULL
        cat(paste("Calculating bootstrap standard errors..."))
        flush.console();
        
        
        ## MVN simulation from the Bayesian approximation to the 
        ## posterior dsitribution of the model parameters
        
        if (se.method=='posterior') {
            boot.res <- posterior.sim.gradients(mod, n.boot, 
                                      phenotype,covariates)
            cat(paste("done.",'\n')); flush.console();


	   ## case bootstrapping

        } else if (se.method == 'boot.case') {
            boot.rep.1 <- min(100,n.boot)
 
             a<-Sys.time()
            ttc <- system.time(
                boot.res.1 <- boot(
                    data = mod$model,
                    statistic = boot.gradients,
                    R = boot.rep.1,
                    parallel = parallel,
                    ncpus = ncpus,
                    mod = mod,
                    phenotype=phenotype,
                    covariates=covariates
                )
            )
            b<-Sys.time()
            if (n.boot <= 100) {
                cat(paste('\n',"     ... time used: ",
                                   ttc['sys.self'],"...",'\n')); 
                                  flush.console();
                boot.res <- boot.res.1$t
            } else {
                cat(paste('\n',"     ... estimated completion at ",
                             b+(b-a)/100*(n.boot-100),"..."
                                           )); flush.console();
                boot.res.2 <- boot(
                    data = mod$model,
                    statistic = boot.gradients,
                    R = n.boot-100,
                    parallel = parallel,
                    ncpus = ncpus,
                    mod = mod,
                    phenotype = phenotype,
                    covariates = covariates)
                boot.res <- rbind(boot.res.1$t,boot.res.2$t)
            }
            cat(paste("done.",'\n')); flush.console();

	   ## parametric bootstrap

        } else if (se.method == 'boot.para') {
            boot.rep.1 <- min(100,n.boot)
            a<-Sys.time()
            ttc <- system.time(
                boot.res.1 <- boot(
                    data = mod$model,
                    statistic = boot.gradients.para,
                    R = boot.rep.1,
                    parallel = parallel,
                    ncpus = ncpus,
                    mod = mod,
                    phenotype=phenotype,
                    covariates=covariates
                )
            )
            b<-Sys.time()
            if (n.boot <= 100) {
                cat(paste('\n',"     ... time used: ",
                                       ttc['sys.self'],"...",'\n')); 
                               flush.console();
                boot.res <- boot.res.1$t
            } else {
                cat(paste('\n',"     ... estimated completion at ",
                            b+(b-a)/100*(n.boot-100),"..."
                                                 )); flush.console();
                boot.res.2 <- boot(
                    data = mod$model,
                    statistic = boot.gradients.para,
                    R = n.boot-100,
                    parallel = parallel,
                    ncpus = ncpus,
                    mod = mod,
                    phenotype = phenotype,
                    covariates=covariates
                )
                boot.res <- rbind(boot.res.1$t,boot.res.2$t)
            }
            cat(paste("done.",'\n')); flush.console();
        }
        
        
        ## estimation of SEs from bootstrap (or simulation from)
        ## approximate posterior
        
        ests[['SE']]<-sapply(as.data.frame(boot.res),sd)
        ests[['P-value']]<-sapply(
            X = as.data.frame(boot.res),
            FUN = function(x) {
                2*min(sum((x>0)+0), sum((x<0)+0))/length(x)
            }
        )

	## permutation test (P-values, no SEs)

    } else if (se.method=='permute') {
        perm.res<-NULL
        cat(paste("Calculating permutation P-values..."))
        flush.console();
        
        if(is.null(covariates)==FALSE){
        	 cat(paste("Warning: fitness is permuted across",
        	                        "all values of covariate(s)."))
        	 flush.console()
        }
    	
            boot.rep.1 <- min(100,n.boot)
            perm.res<-NULL
            a<-Sys.time()
            ttc <- system.time(
                perm.res.1 <- boot(
                    data = mod$model,
                    statistic = permute.W.refit,
                    R = boot.rep.1,
                    parallel = parallel,
                    ncpus = ncpus,
                    mod = mod,
                    phenotype=phenotype,
                    covariates=covariates
                )
            )
            b<-Sys.time()
            if (n.boot <= 100) {
                cat(paste('\n',"     ... time used: ",
                                    ttc['sys.self'],"...",'\n')); 
                                flush.console();
                perm.res <- perm.res.1$t
            } else {
                cat(paste('\n',"     ... estimated completion at ",
                               b+(b-a)/100*(n.boot-100),"..."
                                            )); flush.console();
                perm.res.2 <- boot(
                    data = mod$model,
                    statistic = permute.W.refit,
                    R = n.boot-100,
                    parallel = parallel,
                    ncpus = ncpus,
                    mod = mod,
                    phenotype = phenotype,
                    covariates=covariates
                    )
                perm.res <- rbind(perm.res.1$t,perm.res.2$t)
          	  cat(paste("done.",'\n')); flush.console();
            }
        
        
        ests[['SE']]<-NA
            n.grads<-dim(perm.res)[2]-1
        ests[['P-value']]<-sapply(
            X = as.data.frame(perm.res[,1:n.grads]-matrix(ests[['estimates']],
                   (dim(perm.res)[1]),length(ests[['estimates']]),byrow=TRUE)),
            FUN = function(x) {
                2*min(sum((x>0)+0), sum((x<0)+0))/length(x)
            }
        )
    	
    
    } else if (se.method=='n') {
        ests[['SE']] <- as.numeric(rep(NA,length(ests$estimates)))
        ests[['P-value']] <- as.numeric(rep(NA,length(ests$estimates)))
    } else {
        warning("Method for calculating standard errors not recognized.")
        ests[['SE']] <- as.numeric(rep(NA,length(ests$estimates)))
        ests[['P-value']] <- as.numeric(rep(NA,length(ests$estimates)))
    }

    # Label results nicely:
    ests<-as.data.frame(ests)
    rownames(ests)[1:nTraits]<-paste("B-",phenotype,sep="")
    rownames(ests)[(nTraits+1):(2*nTraits)]<-paste("G-",phenotype,sep="")

    if (nTraits > 1) {
        cor.names<-function(x,y){paste("G",x,y,sep="-")}
        rownames(ests)[(2*nTraits+1):(2*nTraits+(nTraits^2-nTraits)/2)] <-
            outer(phenotype, phenotype,cor.names)[
                upper.tri(outer(phenotype, phenotype,cor.names))
            ]
    }
    
    res<-ests
    if(se.method %in% c('posterior','boot.case','boot.para')) {
    	                          res<-list(ests=ests,boot=boot.res)}
    if(se.method %in% c('permute')) res<-list(ests=ests,perm=perm.res)
    
    res
}
