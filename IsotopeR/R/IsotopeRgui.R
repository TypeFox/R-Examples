#IsotopeR model. Bayesian Inference on stable istope analysis.
#code based on Semmens et al 2009
#see IsotopeR instructions.doc for details on running the model.
#see Hopkins and Ferguson 2011 for details on the model structure.

#for debugging
#source('IsotopeRModelsNoGroups.R')
#source('IsotopeRModelsGroups.R')


# source('IsotopeRgui.R'); IsoWrapper(Mixtures="/home/troutinthemilk/Desktop/Mixtures_yrs2.csv", Sources="/home/troutinthemilk/Desktop/Sources.csv", Concentrations="Optional File", Discrimination.Error="Optional File", output.name="SampleOutput.Rdata", mcmc.chains=3, mcmc.burn=1e3, mcmc.chainLength=1e3, mcmc.thin=1, plot.observations=F,  plot.mixing.estimates=T, plot.dietary.source.contributions=F, color.plots=T, run.parallel = TRUE)
#source('IsotopeRgui.R'); IsoWrapper(Mixtures="Data_Example/4-source/Mixtures.csv", Sources="Data_Example/4-source/Sources.csv", Concentrations="Data_Example/4-source/SourcesCD.csv", Discrimination.Error="Data_Example/4-source/DiscrimSD.csv", output.name="SampleOutput.Rdata", mcmc.chains=3, mcmc.burn=1e3, mcmc.chainLength=1e3, mcmc.thin=1, plot.observations=F,  plot.mixing.estimates=F, plot.dietary.source.contributions=F, color.plots=F, run.parallel = TRUE)
IsoWrapper <- function(Mixtures="Necessary File", Sources="Necessary File", Concentrations="Optional File", Discrimination.Error="Optional File", Measurement.Error="Optional File",  output.name="SampleOutput.Rdata", mcmc.chains=3, mcmc.burn=1e3, mcmc.chainLength=1e3, mcmc.thin=1, plot.observations=TRUE,  plot.mixing.estimates=TRUE, plot.dietary.source.contributions=TRUE, color.plots=TRUE, run.parallel = TRUE) {
   
    mcmc.chainLength    <- as.integer(mcmc.chainLength+mcmc.burn) #total number of iterations per chain (includes burnin)

    #name and location of the model file passed to JAGS
    model.loc   <- "IsotopeR.txt" 

    #parameters to be returned by JAGS
    jags.params <- c("mu.source", "sd.source", "rho.mat", "mu.conc", "sd.conc", "mu.mix", "p", "p.pop", "sd.me","sd.res") #modified 8/7/15
    file.flag <- ""
    noconc.flag = 0
    nome.flag = 0
    nodiscrim.flag = 0
    nodigest.flag=0

    #reads in the files               
    X           <- try(as.matrix(utils::read.table(Mixtures, sep='\t', header=TRUE)), silent=TRUE) #mixture data file
    if(class(X) == 'try-error') { stop("Mixture file not found") }
    if(dim(X)[2] == 1) {
        X           <- as.matrix(utils::read.table(Mixtures, sep=',', header=TRUE)) #mixture data file
    }
	if(dim(X)[2] == 1) {
        X           <- as.matrix(utils::read.table(Mixtures, sep=';', header=TRUE)) #mixture data file
    }
    N 		    <- dim(X)[1] #number of individuals in the sample
    num.iso     <- dim(X)[2]-2 #number of isotopes in the sample

    sources     <- try((utils::read.table(Sources,sep='\t',header=TRUE)), silent=TRUE) #source data
    if(class(sources) == 'try-error') {stop("Source file not found")}
    if(dim(sources)[2] == 1) { 
        sources     <- (utils::read.table(Sources,sep=',',header=TRUE)) #source data
    }
	if(dim(sources)[2] == 1) { 
        sources     <- (utils::read.table(Sources,sep=';',header=TRUE)) #source data
    }

    D   			<- NA
    cd.mat 			<- NA
    subcd.vec 		<- NA
    subcd.samples 	<- NA
    
	if(Concentrations == 'Optional File') { 
        noconc.flag	= 1
        file.flag	<- paste(file.flag,"noconc",sep='') 
        
    } else {
        D           <- try((utils::read.table(Concentrations,sep='\t',header=TRUE)), silent=TRUE) #source concentration data
        if(dim(D)[2] == 1) {
            D           <- try(utils::read.table(Concentrations,sep=',',header=TRUE), silent=TRUE) #source concentration data
        }
        if(dim(D)[2] == 1) {
            D           <- try(utils::read.table(Concentrations,sep=';',header=TRUE), silent=TRUE) #source concentration data
        }
        D <- D[order(D[,num.iso+1]),]

    }

    Z <- NA
    if(Measurement.Error == 'Optional File') { 
        nome.flag = 1
        file.flag   <- paste(file.flag,"nome",sep='') 
    } else {
        Z           <- try(as.matrix(utils::read.table(Measurement.Error,sep='\t',header=TRUE)), silent=TRUE) #file for measurement error
        if(class(Z) == 'try-error') {stop("Measurement Error file not found")}
        if(dim(Z)[2] == 1) { 
            Z           <- as.matrix(utils::read.table(Measurement.Error, sep=',', header=TRUE)) #file for measurement error
        }
        if(dim(Z)[2] == 1) { 
            Z           <- as.matrix(utils::read.table(Measurement.Error, sep=';', header=TRUE)) #file for measurement error
        }        
    }

    discrim.sd <- NA
    if(Discrimination.Error == 'Optional File') {
        nodiscrim.flag = 1
        #file.flag   <- paste(file.flag,"nodiscrim",sep='')
    } else {
        discrim.sd      <- try((utils::read.table(Discrimination.Error,sep='\t',header=TRUE)), silent=TRUE) #file with the standard deviation of discrimination in sources
        if(class(discrim.sd) == 'try-error') {stop("Discrimination Error file not found")}
        if(dim(discrim.sd)[2] == 1) {
            discrim.sd      <- (utils::read.table(Discrimination.Error,sep=',',header=TRUE)) #file with the standard deviation of discrimination in sources
        }
        if(dim(discrim.sd)[2] == 1) {
            discrim.sd      <- (utils::read.table(Discrimination.Error,sep=';',header=TRUE)) #file with the standard deviation of discrimination in sources
        }   
        discrim.sd <- discrim.sd[order(discrim.sd[,num.iso+1]),]
    }

    nodigest.flag = 1

    options(warn=-1)
    #extract useful info from the read in data to pass to JAGS
	num.sources <- nlevels(as.factor(sources[,num.iso+1]))
	num.groups 		<- nlevels(as.factor(X[,num.iso+1]))

	groupnum.mat 	<- matrix(NA,num.groups,2)
	for(i in 1:num.groups) {
		groupnum.mat[i,1] <- min(which(X[,num.iso+1]==i))
		groupnum.mat[i,2] <- max(which(X[,num.iso+1]==i))
	}

	#determine the proper model to run
	if(file.flag == "") { 
		if(num.groups ==1) { curr.model <- IsotopeRfull} else {curr.model <- IsotopeRfullgroup}
	}
    if(file.flag == "noconc") { 
		if(num.groups ==1) { curr.model <- IsotopeRnoconc } else { curr.model <- IsotopeRnoconcgroup}
	}
    if(file.flag == "noconcnome") { 
		if(num.groups == 1) {curr.model <- IsotopeRnoconcnome} else {curr.model <- IsotopeRnoconcnomegroup}
	}
    if(file.flag == "nome") { 
		if(num.groups == 1) { curr.model <- IsotopeRnome } else {curr.model <- IsotopeRnomegroup }
	}

    #prior diet proportions
    alpha <- rep(1,num.sources)/num.sources

    #center Z (observation error) around 0
    Nz <- NA
    if(!nome.flag) {
        Z[,1] <- Z[,1] - mean(Z[,1])
        Z[,2] <- Z[,2] - mean(Z[,2])
        Nz <- dim(Z)[1]
    }
    
    #prior parameters for concentrations
    dmu.prior.tau 	<- diag(num.iso)*1/100
    dmu.prior.mu 	<- matrix(50, num.iso, num.sources)

    #population prior
    alpha.clr <- log(alpha/prod(alpha)^(1/length(alpha))) #transform onto CLR scale

    #measurement error covariance matrix and mean
    tauZ    <- diag(num.iso)
    muz     <- rep(0,num.iso)
    
    #puts sources into an array- makes it easier to use in JAGS (this is just processing stuff and not very important to understand)   
    source.ss <- vector('numeric', num.sources)
    counter <- 1

	for(i in levels(sources[,num.iso+2])) {    
	  source.ss[counter] <- length(which(sources[,num.iso+2] == i))
      counter <- counter+1        
    }    
    names(source.ss) <- levels(sources$source)
    
	##get array indices for the different subsources
	subsources <- sources[,num.iso+2]
	subsource.vec <- nlevels(as.factor(sources[,num.iso+1]))
	counter <- 1
	for(i in levels(as.factor(sources[,num.iso+1]))) {
		subsource.vec[counter] <- nlevels(as.factor(sources[which(sources[,num.iso+1] == i), num.iso+2]))
		counter <- counter+1
	}

	subsource.samples <-array(NA, c(num.sources, max(subsource.vec), 2))
	source.counter <- 1		
	for(i in levels(as.factor(sources[,num.iso+1]))) {
 		curr.source.index <- which(sources[,num.iso+1] == i)
		subsource.counter <- 1
		for(j in levels(as.factor(sources[curr.source.index, num.iso+2]))) {
			curr.subsource.index <- which(sources[curr.source.index, num.iso+2] == j)
			subsource.samples[source.counter, subsource.counter, ] <- 	curr.source.index[1] + curr.subsource.index[c(1, length(curr.subsource.index))] - 1				
			
			subsource.counter <- subsource.counter +1
			
		}
		source.counter <- source.counter + 1
	}
    counter <- 1
    source.mat <- as.matrix(sources[,1:num.iso])
    
    #prior parameters for sources  
    mu.prior.mu <- rep(0,num.iso)#c(0, 0) #apply(sources[,1:num.iso], 2, mean)
    mu.prior.cov <- solve(diag(num.iso)*100) #source mean prior covariance matrix
    cd.array <- NA

    if(!noconc.flag) {
        #bu ild array of isotope concentrations and apply digestability
        cd.mat <- as.matrix(D[,1:num.iso]) #+ rnorm(length(as.matrix(D[,1:num.iso])), 0, 0.01)
         ##get array indices for the different subsource concentrations
		subcd <- D[,num.iso+2]
		subcd.vec <- nlevels(as.factor(D[,num.iso +1]))
		counter <- 1
		for(i in levels(as.factor(D[,num.iso+1]))) {
			subcd.vec[counter] <- nlevels(as.factor(D[which(D[,num.iso+1] == i), num.iso+2]))
			counter <- counter+1
		}

		subcd.samples <-array(NA, c(num.sources, max(subcd.vec), 2))
		source.counter <- 1
		for(i in levels(as.factor(D[,num.iso+1]))) {
 			curr.source.index <- which(D[,num.iso+1] == i)

			subsource.counter <- 1
			for(j in levels(as.factor(D[curr.source.index, num.iso+2]))) {
				curr.subsource.index <- which(D[curr.source.index, num.iso+2] == j)
				
				subcd.samples[source.counter, subsource.counter, ] <- 	curr.source.index[1] + curr.subsource.index[c(1, length(curr.subsource.index))] - 1
				subsource.counter <- subsource.counter +1
			}
			source.counter <- source.counter + 1
		}
    }

  #rescale sources by discrimination variation
    if(!nodiscrim.flag) {
	  for(i in levels(as.factor(sources[,num.iso+1]))) {
        discrim.entry = which(i == discrim.sd[,num.iso+1])
		if(length(discrim.entry) == 0) {stop('No match found between discrimination variation and sources')}
		for(h in 1:num.iso) { 
		  source.mat[which(sources[,num.iso+1] == i),h] = (sources[which(sources[,num.iso+1] == i),h]-mean(sources[which(sources[,num.iso+1] == i),h]))*(1 + discrim.sd[discrim.entry,h]/stats::sd(sources[which(sources[,num.iso+1] == i),h])) + mean(sources[which(sources[,num.iso+1] == i),h])
		}
	  }
    }

    #gets number of individual observations
	dim.x	<- dim(X)
    num.inds 	<- nlevels(as.factor(X[,dim.x[2]]))
	ind.levels	<- as.factor(X[,dim.x[2]])
    ind.counts 	<- vector('numeric',num.inds)

	counter	<- 1
    for(i in ind.levels) {
        ind.counts[counter] <- length(which(i == X[,dim.x[2]]))
		counter				<- counter+1
    }
	
	ind.array 	<- array(NA,c(num.iso,num.inds,max(ind.counts)))

	counter	<- 1
    for(i in levels(ind.levels)) {
        currObs = which(X[,dim.x[2]]==i)
        ind.array[1:num.iso,counter,1:length(currObs)] 	<- t(X[which(X[,dim.x[2]]==i), 1:num.iso])
		counter							<- counter + 1
    }
	rho.flag	<- ifelse(num.iso ==2, 1, 0)

   #individual observation id's
   N <- num.inds
   jags.dump <- list(muz=muz, ind.counts=ind.counts, ind.array=ind.array, N=N, num.sources=num.sources, num.iso=num.iso, Z=Z, dmu.prior.mu=dmu.prior.mu, Nz=Nz, mu.prior.mu=mu.prior.mu, mu.prior.cov=mu.prior.cov, dmu.prior.tau=dmu.prior.tau, alpha.clr=alpha.clr,  subsource.vec=subsource.vec, subsource.samples=subsource.samples, source.mat=source.mat, cd.mat = cd.mat, subcd.vec=subcd.vec, subcd.samples=subcd.samples, num.groups=num.groups, groupnum.mat=groupnum.mat, rho.flag=rho.flag)
    
   if(noconc.flag) {
        jags.rem    <- which( names(jags.dump) == 'cd.mat' | names(jags.dump) == 'dmu.prior.mu' | names(jags.dump) == 'dmu.prior.tau' |  names(jags.dump) == 'subcd.samples' |  names(jags.dump)  == 'subcd.vec')
        jags.dump <- jags.dump[-jags.rem]
        
       jags.rem <- which(jags.params == "mu.conc" | jags.params=="sd.conc") 
	   jags.params <- jags.params[-jags.rem]
	
    }

    if(nome.flag) {
        jags.rem    <- which(names(jags.dump)  == 'Z' | names(jags.dump) == 'muz' | names(jags.dump) == 'tauZ' | names(jags.dump)== 'Nz')
        jags.dump <- jags.dump[-jags.rem]
        
       jags.rem <- which(jags.params == "sd.me") 
	   jags.params <- jags.params[-jags.rem]
    }
    if(num.groups <= 1) {
		jags.rem    <- which(names(jags.dump)  == 'groupnum.mat' | names(jags.dump)  == 'num.groups')
        jags.dump <- jags.dump[-jags.rem]        
	} else {
		jags.params <- c(jags.params,"p.group") 
	}
	if(num.iso != 2) {
		jags.rem 		<- which(jags.params == "rho.mat")
		jags.params	<- jags.params[-jags.rem]
	}
	
	#function used to initialize parameters
    jags.inits <- list( dmu.prior.mu=dmu.prior.mu, mu.prior.mu=mu.prior.mu, p.transform=stats::runif(num.sources), region.sig=0.5, ind.sig=0.5, p.ind = matrix(stats::runif(N*num.sources), N, num.sources) )

	if(run.parallel) {parallel.state <- "parallel"} else { 
		parallel.state <- "interruptible"
		jags.params <- c(jags.params, "dic", 'deviance', 'pd', 'ped') #dic can only be run in when nonparallel calculations are used.
	}
	if(noconc.flag) {jags.adapt=1e3} else{ jags.adapt=1e3}
	jags.out <- runjags::run.jags(model=curr.model, monitor=jags.params, data=jags.dump, adapt=jags.adapt, n.chains=mcmc.chains, burnin=mcmc.burn, sample=(mcmc.chainLength-mcmc.burn), thin=mcmc.thin, psrf.target=1.5, plots=FALSE, silent.jags=FALSE, method=parallel.state)
	if(mcmc.chains > 1) {	
	  r.est <- jags.out$psrf$psrf[,1] 
	  jags.output.mat <- cbind(jags.out$summary$statistics[,1:2], jags.out$summary$quantiles, r.est)
	} else {
	  jags.output.mat <- cbind(jags.out$summary$statistics[,1:2], jags.out$summary$quantiles)	
	}
    
    print(jags.output.mat)
    if(!run.parallel & mcmc.chains>1) {
	  print(jags.out$dic)
	}
    save(jags.out, jags.output.mat, X, sources, nome.flag, num.sources, num.iso, mcmc.chains, N, num.groups, file=output.name)
    utils::write.table(jags.output.mat, file=paste(strsplit(output.name, ".Rdata")[[1]],'.txt',sep=''))
	if(!run.parallel & mcmc.chains>1) {
	  sink(file=paste(strsplit(output.name, ".Rdata")[[1]],'.txt',sep=''), append=TRUE)
	  print(jags.out$dic)
      #print(jags.out$pd)
	  sink()
	}
	
	xlab = switch(dimnames(X)[[2]][1], "delta.13C"=expression(paste(delta^13, 'C')), "d.13C"=expression(paste(delta^13, 'C')), "delta13C"=expression(paste(delta^13, 'C')), "d13C"=expression(paste(delta^13, 'C')), "c13"=expression(paste(delta^13, 'C')), "13c"=expression(paste(delta^13, 'C')), "C13"=expression(paste(delta^13, 'C')), "13C"=expression(paste(delta^13, 'C')),
		"delta.15N"=expression(paste(delta^15, 'N')), "d.15N"=expression(paste(delta^15, 'N')), "delta15N"=expression(paste(delta^15, 'N')), "d15N"=expression(paste(delta^15, 'N')), "n15"=expression(paste(delta^15, 'N')), "15n"=expression(paste(delta^15, 'N')), "N15"=expression(paste(delta^15, 'N')), "15N"=expression(paste(delta^15, 'N')),
		"delta.34S"=expression(paste(delta^34, 'S')), "d.34S"=expression(paste(delta^34, 'S')), "delta13C"=expression(paste(delta^34, 'S')), "d34S"=expression(paste(delta^34, 'S')), "s34"=expression(paste(delta^34, 'S')), "34s"=expression(paste(delta^34, 'S')), "S34"=expression(paste(delta^34, 'S')), "34S"=expression(paste(delta^34, 'S')), dimnames(X)[[2]][1])
	ylab = switch(dimnames(X)[[2]][2], "delta.13C"=expression(paste(delta^13, 'C')), "d.13C"=expression(paste(delta^13, 'C')), "delta13C"=expression(paste(delta^13, 'C')), "d13C"=expression(paste(delta^13, 'C')), "c13"=expression(paste(delta^13, 'C')), "13c"=expression(paste(delta^13, 'C')), "C13"=expression(paste(delta^13, 'C')), "13C"=expression(paste(delta^13, 'C')),
		"delta.15N"=expression(paste(delta^15, 'N')), "d.15N"=expression(paste(delta^15, 'N')), "delta15N"=expression(paste(delta^15, 'N')), "d15N"=expression(paste(delta^15, 'N')), "n15"=expression(paste(delta^15, 'N')), "15n"=expression(paste(delta^15, 'N')), "N15"=expression(paste(delta^15, 'N')), "15N"=expression(paste(delta^15, 'N')),
		"delta.34S"=expression(paste(delta^34, 'S')), "d.34S"=expression(paste(delta^34, 'S')), "delta13C"=expression(paste(delta^34, 'S')), "d34S"=expression(paste(delta^34, 'S')), "s34"=expression(paste(delta^34, 'S')), "34s"=expression(paste(delta^34, 'S')), "S34"=expression(paste(delta^34, 'S')), "34S"=expression(paste(delta^34, 'S')), dimnames(X)[[2]][2])
		
    if(plot.observations) {
		
		if(num.sources >= 3 & num.iso==2) {grDevices::dev.new(); Tri.plots(jags.out, X, sources=sources, plot.ind.flag=TRUE, me.flag=!nome.flag, color.plots=color.plots, xlab=xlab, ylab=ylab) } 
		if(num.iso == 3) { requireNamespace("rgl", quietly = TRUE); rgl::open3d(); RGL.plots(jags.out, X=X, sources=sources, plot.mix=FALSE, color.plots=color.plots) }
		if(num.sources == 2 & num.iso==2) {grDevices::dev.new(); Bi.plots(jags.out, X, sources=sources, plot.ind.flag=TRUE, me.flag=!nome.flag, color.plots=color.plots, xlab=xlab, ylab=ylab)} 
		
    }
    options(warn=0)

    if(plot.mixing.estimates) {
        if(num.sources >= 3 & num.iso==2) {grDevices::dev.new(); Tri.plots(jags.1=jags.out, X=X, sources=sources, plot.mix=TRUE, me.flag=!nome.flag, color.plots=color.plots, xlab=xlab, ylab=ylab) } 
        if(num.iso == 3) { requireNamespace("rgl", quietly = TRUE); rgl::open3d(); RGL.plots(jags.out, X=X, sources=sources, plot.mix=TRUE, color.plots=color.plots) } 
        if(num.sources == 2 & num.iso==2) {grDevices::dev.new(); Bi.plots(jags.out, X,  sources=sources, plot.mix=TRUE, me.flag=!nome.flag, color.plots=color.plots, xlab=xlab, ylab=ylab) } 
    }
    
    if(plot.dietary.source.contributions) {
        grDevices::dev.new()
        curves.plot(jags.1=jags.out, num.sources=num.sources, num.chains=mcmc.chains, color=color.plots, individuals=N,  xlab.vec=levels(as.factor(sources[,num.iso+1])), num.groups=num.groups)
    }

}

##load previously run data and outputs graphs
load.prev.func <- function(file.name="SampleOutput.Rdata", plot.observations=TRUE, plot.mixing.estimates=TRUE, plot.dietary.source.contributions=TRUE, color.plots=TRUE) {

	jags.out=NA
	sources=NA
	nome.flag=NA
	num.sources=NA
	mcmc.chains=NA
	N=NA
	X=NA
	num.iso=NA
	num.groups=NA
	jags.output.mat=NA
	load(file=file.name)

	if(!exists("jags.out") | !exists("X") | !exists("sources") | !exists("nome.flag") | !exists("num.iso") | !exists("num.sources") | !exists("mcmc.chains") | !exists("N") | !exists("num.groups") | !exists("jags.output.mat"))
	{stop(".Rdata file error")}
	xlab = switch(dimnames(X)[[2]][1], "delta.13C"=expression(paste(delta^13, 'C')), "d.13C"=expression(paste(delta^13, 'C')), "delta13C"=expression(paste(delta^13, 'C')), "d13C"=expression(paste(delta^13, 'C')), "c13"=expression(paste(delta^13, 'C')), "13c"=expression(paste(delta^13, 'C')), "C13"=expression(paste(delta^13, 'C')), "13C"=expression(paste(delta^13, 'C')),
			"delta.15N"=expression(paste(delta^15, 'N')), "d.15N"=expression(paste(delta^15, 'N')), "delta15N"=expression(paste(delta^15, 'N')), "d15N"=expression(paste(delta^15, 'N')), "n15"=expression(paste(delta^15, 'N')), "15n"=expression(paste(delta^15, 'N')), "N15"=expression(paste(delta^15, 'N')), "15N"=expression(paste(delta^15, 'N')),
			"delta.34S"=expression(paste(delta^34, 'S')), "d.34S"=expression(paste(delta^34, 'S')), "delta13C"=expression(paste(delta^34, 'S')), "d34S"=expression(paste(delta^34, 'S')), "s34"=expression(paste(delta^34, 'S')), "34s"=expression(paste(delta^34, 'S')), "S34"=expression(paste(delta^34, 'S')), "34S"=expression(paste(delta^34, 'S')), dimnames(X)[[2]][1])
	ylab = switch(dimnames(X)[[2]][2], "delta.13C"=expression(paste(delta^13, 'C')), "d.13C"=expression(paste(delta^13, 'C')), "delta13C"=expression(paste(delta^13, 'C')), "d13C"=expression(paste(delta^13, 'C')), "c13"=expression(paste(delta^13, 'C')), "13c"=expression(paste(delta^13, 'C')), "C13"=expression(paste(delta^13, 'C')), "13C"=expression(paste(delta^13, 'C')),
			"delta.15N"=expression(paste(delta^15, 'N')), "d.15N"=expression(paste(delta^15, 'N')), "delta15N"=expression(paste(delta^15, 'N')), "d15N"=expression(paste(delta^15, 'N')), "n15"=expression(paste(delta^15, 'N')), "15n"=expression(paste(delta^15, 'N')), "N15"=expression(paste(delta^15, 'N')), "15N"=expression(paste(delta^15, 'N')),
			"delta.34S"=expression(paste(delta^34, 'S')), "d.34S"=expression(paste(delta^34, 'S')), "delta13C"=expression(paste(delta^34, 'S')), "d34S"=expression(paste(delta^34, 'S')), "s34"=expression(paste(delta^34, 'S')), "34s"=expression(paste(delta^34, 'S')), "S34"=expression(paste(delta^34, 'S')), "34S"=expression(paste(delta^34, 'S')), dimnames(X)[[2]][2])
			
	if(plot.observations) {
			
			if(num.sources >= 3 & num.iso==2) {dev.new(); Tri.plots(jags.out, X, sources=sources, plot.ind.flag=TRUE, me.flag=!nome.flag, color.plots=color.plots, xlab=xlab, ylab=ylab) } else {
			if(num.iso == 3) { requireNamespace("rgl", quietly = TRUE); rgl::open3d(); RGL.plots(jags.out, X=X, sources=sources, plot.mix=FALSE, color.plots=color.plots) } else {
			if(num.sources == 2 & num.iso==2) { dev.new(); Bi.plots(jags.out, X, sources=sources, plot.ind.flag=TRUE, me.flag=!nome.flag, color.plots=color.plots, xlab=xlab, ylab=ylab)} else {warning("No observation plot available for this number of isotopes", call.=FALSE) }
		}
		}
		# 		  if(num.iso==2) {dev.new(); Tri.plots(jags.out, X, sources=sources, plot.ind=TRUE, me.flag=!nome.flag, color.plots=color.plots) } 
    }

    if(plot.mixing.estimates) {
        if(num.iso == 2 & num.sources >= 3) {dev.new(); Tri.plots(jags.1=jags.out, X=X, sources=sources, plot.mix=TRUE, me.flag=!nome.flag, xlab=xlab, ylab=ylab) } else {
        if(num.iso == 3) { requireNamespace("rgl", quietly = TRUE); rgl::open3d(); RGL.plots(jags.out, X=X, sources=sources, plot.mix=TRUE, color.plots=color.plots) } else {
        if(num.iso == 2 & num.sources==2) {dev.new(); Bi.plots(jags.out, X,  sources=sources, plot.mix=TRUE, me.flag=!nome.flag, xlab=xlab, ylab=ylab) } else {warning("Warning: No mixing plot available for this isotope/source combination", call.=FALSE)}
        }
        }    
    }
    
 
    if(plot.dietary.source.contributions) {
        dev.new()
        curves.plot(jags.1=jags.out, num.sources=num.sources, num.chains=mcmc.chains, color=color.plots, individuals=N,  xlab.vec=levels(as.factor(sources[,num.iso+1])), num.groups=num.groups)
    }
    print(jags.output.mat)
	return(NULL)
}

IsotopeR    <- function() {
    if(interactive()) {
		fgui::fguiWindow( basicMenu=FALSE, title="IsotopeR 0.5", text="Please choose an option from the Analysis menu." )
		#win <- mgui(IsoWrapper, argFilter = list(Mixtures="{{} {.csv}}", Sources="{{} {.csv}}", Concentrations="{{} {.csv}}", Measurement.Error="{{} {.csv}}",   Discrimination.Error="{{} {.csv}}", Digestibility.Factor="{{} {.csv}}"), argText = list(mcmc.chains="number of chains", mcmc.burn="MCMC burnin", mcmc.chainLength="MCMC runs", mcmc.thin="thinning rate", output.name="Output file"), argOption = list(run.parallel=c("TRUE", "FALSE"), plot.observations=c("TRUE","FALSE"), plot.dietary.source.contributions= c("TRUE", "FALSE"), plot.mixing.estimates=c('TRUE','FALSE'), color.plots=c("TRUE", "FALSE")), closeOnExec=TRUE, title=c("Analysis","New Run"), exec="Run IsotopeR", helpsFunc="IsoWrapper", output=NULL)
		win <- fgui::mgui(IsoWrapper, argFilter = list(Mixtures="{{} {.csv}}", Sources="{{} {.csv}}", Concentrations="{{} {.csv}}", Discrimination.Error="{{} {.csv}}", Measurement.Error="{{} {.csv}}"), argText = list(mcmc.chains="number of chains", mcmc.burn="MCMC burnin", mcmc.chainLength="MCMC runs", mcmc.thin="thinning rate", output.name="Output file"), argOption = list(run.parallel=c("TRUE", "FALSE"), plot.observations=c("TRUE","FALSE"), plot.dietary.source.contributions= c("TRUE", "FALSE"), plot.mixing.estimates=c('TRUE','FALSE'), color.plots=c("TRUE", "FALSE")), closeOnExec=TRUE, title=c("Analysis","New Run"), exec="Run IsotopeR", helpsFunc="IsoWrapper", output=NULL)
        
		win <- fgui::mgui(load.prev.func, title=c("Analysis","Load Previous Run"), argFilter= list(file.name="{{} {.Rdata}}"), argOption = list(plot.observations=c("TRUE","FALSE"), plot.dietary.source.contributions= c("TRUE", "FALSE"), plot.mixing.estimates=c('TRUE','FALSE'), color.plots=c("TRUE", "FALSE")), closeOnExec=TRUE , exec="Plot", output=NULL, helpsFunc="IsoWrapper")
    }


}

