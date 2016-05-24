
gppr.gradients <- function(
    mod,
    phenotype,
    covariates=NULL,
    standardized=FALSE,
    se.method='boot.para',
    n.boot=1000,
    parallel="no",
    ncpus=1
) {
	
  ## sub-functions for gam.gradients()
	
  Wbar <- function(x, mod, phenotype,covariates) {
    new.d <- as.data.frame(mod$data[,
            c(phenotype,covariates)])
    names(new.d)<-c(phenotype,covariates)
    new.d2<-new.d
    for (i in 1:length(x)) {
      new.d[,as.character(phenotype[i])] <-
         new.d[,as.character(phenotype[i])]+x[i]
    }
    p <- predict(
             mod,newdata=new.d,
             type="response")
    return(mean(p))
  }

  gradients <- function(m, phenotype, covariates) {
    nTraits = length(phenotype)

    sds<-apply(as.matrix(mod$data[,phenotype]),2,sd)
    h<-sds*0.02
    derivs<-finite.dif(Wbar,x=rep(0, nTraits), h=h,
      mod=m, phenotype=phenotype, covariates=covariates)

    beta<-derivs$grad.est/derivs$mu
    gamma<-derivs$hessian.est/derivs$mu

    	
    if(standardized){
    	sds<-apply(as.matrix(mod$gen$data[,phenotype]),2,sd)
    	beta<-beta*sds
    	gamma<-gamma*outer(sds,sds)
    }
    return(list(beta=beta,gamma=gamma))
  }

	# case bootstrap
	boot.gradients <- function(data, original, mod, 
	                             phenotype,covariates) {
    d <- as.data.frame(data)[original,]

    mod.prime <- gppr(
        y=mod$gen$y,
        xterms=mod$xterms,
        data=d,
        family=mod$family[[1]],
        nterms=mod$nterms,
        max.terms=mod$max.terms,
        tol=mod$tol
    )

    g <- gradients(mod.prime, phenotype,covariates)

    return(c(g$beta,diag(g$gamma),
        g$gamma[upper.tri(g$gamma)]))
	}

	# for parametric bootstrap
	simulate.gppr.gsg<-function(mod){
		sim.lin.pred<-predict(mod,type='link')
		s<-NULL
		if(mod$family[[1]]=='poisson'|mod$family[[1]]=='Poisson') {
			s<-rpois(length(sim.lin.pred),exp(sim.lin.pred))
		}
		if(mod$family[[1]]=='binomial'|mod$family[[1]]=='Binomial') {
			s<-rbinom(length(sim.lin.pred),1,
			                     inv.logit(sim.lin.pred))
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
    d[,as.character(mod$y)]<-simulate.gppr.gsg(mod)

    mod.prime <- gppr(
        y=mod$y,
        xterms=mod$xterms,
        data=d,
        family=mod$family[[1]],
        nterms=mod$nterms,
        gcvpen=mod$gcvpen,
        maxit=mod$maxit,
        max.terms=mod$max.terms
    )
    
    plot(mod.prime$ppr)

    g <- gradients(mod.prime, phenotype, covariates)

    return(c(g$beta,diag(g$gamma),
        g$gamma[upper.tri(g$gamma)]))
	}

	permute.W.refit<-function(data, original, mod, 
	                              phenotype, covariates){

    d <- mod$data[,c(mod$y,mod$ppr$xterms)]
    d[,mod$gen$y] <- d[ sample(1:(dim(d)[1]),dim(d)[1]),
	         	               mod$gen$y ]

    mod.prime <- gppr(
        y=mod$y,
        xterms=mod$xterms,
        data=d,
        family=mod$family[[1]],
        nterms=mod$nterms,
        gcvpen=mod$gcvpen,
        maxit=mod$maxit,
        max.terms=mod$max.terms
    )


    p<- predict(mod.prime,type="response")
    g <- gradients(mod.prime, phenotype, covariates)
    return(c(g$beta,diag(g$gamma),
            g$gamma[upper.tri(g$gamma)],var(p)))
	}


	## end sub-functions for gam.gradients()	


	## some checks
		
    cls <- class(mod)
    rightType <- any(cls %in% c("gppr"))
    if (!rightType) {
        stop("Argument 'mod' must be a generalized projection pursuit model.")
    }
    if (!is.character(phenotype)) {
        stop("Argument 'phenotype' must be a character vector of terms from the fitness model.")
    }
    if (!is.logical(standardized)) {
        stop("Argument 'standardized' must be logical.")
    }
    if (!is.character(se.method) ) {
        stop("Argument 'se.method' must be character and one of 'boot.para', 'boot.case, 'permute', or 'n'.")
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

    termLabels <- attr(x=terms(mod$ppr), which='term.labels')
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
    if (se.method %in% c('boot.case','boot.para')) {
        boot.res<-NULL
        cat(paste("Calculating bootstrap standard errors..."))
        flush.console();
        
        


	   ## case bootstrapping

        if (se.method == 'boot.case') {
            boot.rep.1 <- min(100,n.boot)
 
             a<-Sys.time()
            ttc <- system.time(
                boot.res.1 <- boot(
                    data = mod$data,
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
                    data = mod$data,
                    statistic = boot.gradients,
                    R = n.boot-100,
                    parallel = parallel,
                    ncpus = ncpus,
                    mod = mod,
                    phenotype = phenotype)
                boot.res <- rbind(boot.res.1$t,boot.res.2$t)
            }
            cat(paste("done.",'\n')); flush.console();

	   ## parametric bootstrap

        } else if (se.method == 'boot.para') {
            boot.rep.1 <- min(100,n.boot)
            a<-Sys.time()
            ttc <- system.time(
                boot.res.1 <- boot(
                    data = mod$data,
                    statistic = boot.gradients.para,
                    R = boot.rep.1,
                    parallel = parallel,
                    ncpus = ncpus,
                    mod = mod,
                    phenotype=phenotype,covariates=covariates
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
                    data = mod$data,
                    statistic = boot.gradients.para,
                    R = n.boot-100,
                    parallel = parallel,
                    ncpus = ncpus,
                    mod = mod,
                    phenotype = phenotype,covariates=covariates
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
                    data = mod$data,
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
                    data = mod$data,
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


