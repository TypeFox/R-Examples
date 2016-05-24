fitness.landscape <- function(
    mod,
    phenotype,
    covariates=NULL,
    points = NULL,
    plt.density = 25,
    PI.method = 'boot.para',
    PI.interval = c(0.25,0.75),
    n.boot = 50,
    refit.smooth=FALSE,
    parallel = "no",
    ncpus = 1
	) {
	
	Wbar <- function(x, mod, phenotype,covariates) {
    		new.d <- as.data.frame(mod$model[,c(phenotype,covariates)])
    		names(new.d)<-c(phenotype,covariates)
    		for (i in 1:length(x)) {
    		    new.d[,as.character(phenotype[i])] <-
    		    new.d[,as.character(phenotype[i])]+x[i]
    		}
    
    		p <- predict.gam(
    	 	   object= mod,
    	 	   newdata=new.d,
    	 	   newdata.guaranteed=TRUE,
    	 	   type="response"
    		)
    	return(mean(p))
	}

	fit.land<-function(mod, points, phenotype,covariates){
    	Wbar.predictions<-array(NA,dim(points)[1])
    	for(a in 1:(dim(points)[1])){
        	x<-array(0,dim=length(phenotype))
        	for(b in 1:length(phenotype)){
            	x[b]<-points[a,b]-mean(mod$model[,as.character(phenotype[b])])
        	}
        	Wbar.predictions[a]<-Wbar(x=x, mod=mod,
        				phenotype=phenotype,covariates=covariates)
    	}
    	return(Wbar.predictions)
	}


	case.boot.fit.land<-function(data, original, mod, points, phenotype,covariates){
    	d <- data[original,]
    	mod.prime <- gam(as.formula(mod$formula), data=d, family=mod$family)
    	l <- fit.land(mod=mod.prime, points=points, 
    				phenotype=phenotype,covariates=covariates)
    	return(l)
	}


	para.boot.fit.land<-function(data, original, mod, points, phenotype,covariates){

		simulate.gam.gsg<-function(mod){
			sim.lin.pred<-predict(mod,type="response")
			s<-NULL
			if(mod$family[[1]]=='poisson'&mod$family[[2]]=='log') {
#				s<-rpois(length(sim.lin.pred),exp(sim.lin.pred))
				s<-rpois(length(sim.lin.pred),sim.lin.pred)
			}
			if(mod$family[[1]]=='binomial') {
#				s<-rbinom(length(sim.lin.pred),1,inv.logit(sim.lin.pred))
				s<-rbinom(length(sim.lin.pred),1,sim.lin.pred)
			}		
			s
		}

		d <- as.data.frame(data)
    	d[[attr(attr(mod$terms,'factors'),'dimnames')[[1]][1]]]<-
    	                                     simulate.gam.gsg(mod)

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

    	l <- fit.land(mod=mod.prime, points=points, 
    			phenotype=phenotype,covariates=covariates)
    	return(l)
	}	
	
	## checks
	cls <- class(mod)
    rightType <- any(cls %in% c("gam"))
    if (!rightType) {
        stop("Argument 'mod' must be a generalized additive model.")
    }
    if (!is.character(phenotype)) {
        stop("Argument 'phenotype' must be a character vector of terms from the fitness model.")
    }
    if (!is.character(PI.method) ) {
        stop("Argument 'PI.method' must be character and one of 'boot.para', 'boot.case', or 'n'.")
    }
    if (!is.vector(x=PI.interval, mode='numeric')) {
        stop("Argument 'PI.interval' must be a numeric vector.")
    }
    if (any(PI.interval < 0 | PI.interval > 1)) {
        stop("Argument 'PI.interval' must have values on the closed interval [0,1].")
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



    if(is.null(points)){
        pts<-function(z){seq(mean(z)-sd(z),mean(z)+sd(z),length.out=plt.density)}
        if ( length(phenotype)==1 ) {
            points<-expand.grid(pts(mod$model[,as.character(phenotype)]))
        } else if (length(phenotype)==2) {
            points<-expand.grid(
                pts(mod$model[,as.character(phenotype[1])]),
                pts(mod$model[,as.character(phenotype[2])])
            )
        } else if (length(phenotype) > 2) {
            stop("We only do up to 2-D fitness landscapes.")
        }
    }

    res <- list(points=points,Wbar=fit.land(mod,points=points,phenotype,covariates))

    if (PI.method=='n') {
        res[['WbarPI']] <- NA
    } else if (PI.method=='boot') {
        cat(paste("Calculating bootstrap prediction intervals...")); flush.console();
        a<-Sys.time()
        ttc <- system.time(
            boot.res.1 <- boot(
                data = mod$model,
                statistic = case.boot.fit.land,
                R = min(100,n.boot),
                parallel = parallel,
                ncpus = ncpus,
                mod = mod,
                points = points,
                phenotype = phenotype,
                covariates = covariates
            )
        )
        b<-Sys.time()
        if (n.boot <= 100) {
            cat(paste('\n',"     ... time used: ",
                      ttc['sys.self'],"...")); flush.console();
            boot.res <- boot.res.1$t
        } else {
            cat(paste('\n',"     ... estimated completion at ",
                b+(b-a)/100*(n.boot-100),"...")); flush.console();
            boot.res.2 <- boot(
                data = mod$model,
                statistic = case.boot.fit.land,
                R = n.boot-100,
                parallel = parallel,
                ncpus = ncpus,
                mod = mod,
                points = points,
                phenotype = phenotype,
                covariates = covariates
            )
            boot.res <- rbind(boot.res.1$t,boot.res.2$t)
            cat(paste('\n',"     ... done.",'\n')); flush.console();
        }
        WbarPI <- sapply(as.data.frame(boot.res), quantile, probs=PI.interval)
        res[['WbarPI']] <- WbarPI
    }  else if (PI.method=='boot.para') {
        cat(paste("Calculating bootstrap prediction intervals..."
                        )); flush.console();
        a<-Sys.time()
        ttc <- system.time(
            boot.res.1 <- boot(
                data = mod$model,
                statistic = para.boot.fit.land,
                R = min(100,n.boot),
                parallel = parallel,
                ncpus = ncpus,
                mod = mod,
                points = points,
                phenotype = phenotype,
                covariates = covariates
            )
        )
        b<-Sys.time()
        if (n.boot <= 100) {
            cat(paste('\n',"     ... time used: ",
                    ttc['sys.self'],"...")); flush.console();
            boot.res <- boot.res.1$t
        } else {
            cat(paste('\n',"     ... estimated completion at ",
               b+(b-a)/100*(n.boot-100),"...")); flush.console();
            boot.res.2 <- boot(
                data = mod$model,
                statistic = para.boot.fit.land,
                R = n.boot-100,
                parallel = parallel,
                ncpus = ncpus,
                mod = mod,
                points = points,
                phenotype = phenotype,
                covariates = covariates
            )
            boot.res <- rbind(boot.res.1$t,boot.res.2$t)
            cat(paste('\n',"     ... done.",'\n')); flush.console();
        }
        WbarPI <- sapply(as.data.frame(boot.res), quantile, probs=PI.interval)
        res[['WbarPI']] <- WbarPI
    } else {
        msg <- paste("Argument 'PI.method' value, '", PI.method, 
              "', not recognized, PI not calculated.", sep='')
        warning(msg)
        res[['WbarPI']] <- NA
    }

    return(res)
}
