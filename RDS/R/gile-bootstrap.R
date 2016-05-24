
SSBS.estimates <- function(rds.data,trait.variable,
		number.of.bootstrap.samples=NULL,
		number.ss.samples.per.iteration=NULL,
		refs.to.trait=NULL,
		N=NULL, 
		confidence.level=NULL,
		fast=TRUE, useC=FALSE,
		control=control.rds.estimates(),
		continuous=NA, is.cts=FALSE, is.quantile=FALSE,
		weight.type="Gile's SS", verbose=TRUE)
{
	
	if(!(trait.variable %in% names(rds.data))){
		stop(sprintf("No variable called %s appears in the data.",
						trait.variable))
	}
	network.size <- attr(rds.data, "network.size.variable")
	
	if(is.null(N)){
		N <- attr(rds.data, "population.size.mid")
		if(is.null(N)){
			N <- ceiling(nrow(rds.data)/0.04)
			warning(paste("Parameter N missing, with no default for this data set. Using N =",N,"\n"))
		}
	}
	
	if(is.null(confidence.level)){
		confidence.level <- 0.95
	}
	if(is.null(number.of.bootstrap.samples)){ 
		number.of.bootstrap.samples <- 500
	}
	if(is.null(number.ss.samples.per.iteration)){
		number.ss.samples.per.iteration <- 500
	}
	if(weight.type=="Good-Fellows"){useC <- FALSE}
	
	trait.min.value <- NULL
	tv <- rds.data[[trait.variable]]
	#rds.data[[trait.variable]] <- tv     ###??? why? Are elements of trait.variable all numbers?
	remvalues <- !is.na(rds.data[[trait.variable]])
	if(sum(remvalues) < nrow(rds.data)){
		warning(paste(nrow(rds.data)-sum(remvalues),"of",nrow(rds.data),
						"trait variable values were missing and were removed."), call. = FALSE)
		if(is.numeric(rds.data[[trait.variable]])){
			trait.min.value <- min(rds.data[[trait.variable]],na.rm=TRUE)-1
			rds.data[[trait.variable]][!remvalues] <- trait.min.value
		}else{
			a=as.vector(rds.data[[trait.variable]])
			a[!remvalues] <- "NA"
			outcome <- factor(a,exclude=NULL)
# 			Make sure the factor labels are alphabetic!
			outcome=factor(outcome,levels=levels(outcome)[order(levels(outcome))],exclude=NULL)
			rds.data[[trait.variable]] <- outcome
			trait.min.value <- "NA"
		}
	}
	

	recruiter.id <- get.rid(rds.data)#as.character(rds.data[[attr(rds.data,"recruiter.id")]])
	seed.rid <- get.seed.rid(rds.data)
	outcome <- factor(rds.data[[trait.variable]],exclude=NULL)
# 	Make sure the factor labels are alphabetic!
	outcome=factor(outcome,levels=levels(outcome)[order(levels(outcome))])
	outclasses <- levels(outcome)
	deg <- rds.data[[network.size]]
	deg <- as.numeric(as.vector(deg))
	g <- length(outclasses)
	
	if(is.null(refs.to.trait)){
		# matrix of numbers of referrals to diseased or non-diseased nodes
		refs.to.trait<-matrix(0,nrow=nrow(rds.data),ncol=g)
		for(j in 1:g){
			# b is the number of referrals for each recruiter to those of group j
			b <- tapply(outcome==outclasses[j],recruiter.id,sum,na.rm=TRUE)[-1]
			# match the names of the recruiters to the id
			#d=match(names(b),as.character(rds.data[[attr(rds.data,"id")]]))
			d <- match(names(b),get.id(rds.data))
			# match the names of the recruiters to the id
			refs.to.trait[d[!is.na(d)],j] <- b[!is.na(d)]
		}
	}else{
		if(is.character(refs.to.trait)){
			TT<-matrix(0,nrow=nrow(rds.data),ncol=g)
			for(j in 1:g){
				TT[,j] = as.vector(rds.data[[refs.to.trait[j]]])
			}
			refs.to.trait<-TT
			rm(TT)
		}
	}
	
	n0 <- sum(recruiter.id==seed.rid, na.rm=TRUE)
	recID <- match(recruiter.id,get.id(rds.data))

	
	###########################################################################################
# Check for missing values and warn the user if any are removed.   This should really taken
# care of elsewhere.  NB: It is also worth considering the semantics of the message 
# "values were missing and were removed".    
	
	remvalues <- !is.na(as.vector(outcome))
	remvalues[as.vector(outcome)=="NA"] <- FALSE
	if(sum(remvalues) < nrow(rds.data)){
		warning(paste(nrow(rds.data)-sum(remvalues),"of",nrow(rds.data),
						"values were missing and were removed."))
		deg <- deg[remvalues]
		outcome <- factor(as.vector(outcome[remvalues]),exclude=NULL)
# 		Make sure the factor labels are alphabetic!
		outcome=factor(outcome,levels=levels(outcome)[order(levels(outcome))],exclude=NULL)
		outclasses <- levels(outcome)
		g <- length(outclasses)
		refs.to.trait <- refs.to.trait[remvalues,]
		recruiter.id <- recruiter.id[remvalues]
		recID <- recID[remvalues]
		n0 <- sum(recruiter.id==seed.rid, na.rm=TRUE)
	}
	
	if(length(unique(outcome))>1){
		result<-sppsboot4ppsreal(degs=deg,dis=outcome,n=N, n0=n0,
				refs.to.trait=refs.to.trait,
				nit=number.of.bootstrap.samples,
				number.ss.samples.per.iteration=number.ss.samples.per.iteration,
				fast=fast,useC=useC,
				weight.type=weight.type,
				recID=recID,
				control=control,
				continuous=continuous, is.cts=is.cts, is.quantile=is.quantile,
				verbose=verbose) 
		if(result$is.cts){
		 names(result$point_estimate) <-  trait.variable
		 names(result$se_estimate) <-  trait.variable
		}else{
		 names(result$point_estimate) <- outclasses
		 names(result$se_estimate) <- outclasses
		 colnames(result$bsests) <- outclasses
		}
	}else{
		result<-list(point_estimate=c(0,1),
				bsests=cbind(rep(1,number.of.bootstrap.samples),
						rep(1,number.of.bootstrap.samples)),
			        is.quantile=is.quantile, is.cts=is.cts,
			        mu=NULL, sigma2=NULL,
				se_estimate=c(0,0))
		names(result$point_estimate) <- c("0",outclasses)
		names(result$se_estimate) <- c("0",outclasses)
		colnames(result$bsests) <- c("0",outclasses)
	}
	
	observed.estimate <- result$point_estimate
	crit <- qnorm((1-confidence.level)/2,lower.tail=FALSE)
	a <- cbind(observed.estimate,
			observed.estimate-crit*apply(result$bsests, 2, sd,na.rm=TRUE),
			observed.estimate+crit*apply(result$bsests, 2, sd,na.rm=TRUE))
	colnames(a)[1] <- trait.variable
	
	attr(a,"bsresult") <- result
	attr(a,"is.cts") <- result$is.cts
	attr(a,"is.quantile") <- result$is.quantile
	attr(a,"mu") <- result$mu
	attr(a,"sigma2") <- result$sigma2
	return(a)
}



#was vh_est
# computes horvitz-thompson estimator
vh.est<-function(degs,dis,wstart=0,wsample=NULL){
	degs[degs==0]<-1 #added 0707
	if(wstart>0){
		degsnew<-degs[wsample>=wstart]
		disnew<-dis[wsample>=wstart]		
	}else{
		degsnew<-degs
		disnew<-dis
	}
	if(is.factor(dis)){
		num <- rep(0,length=length(levels(disnew)))
		a <- tapply(1/degsnew,as.numeric(disnew),sum)
		num[as.numeric(names(a))] <- a
	}else{
		num<-sum(disnew/degsnew)
	}
	den<-sum(1/degsnew)
	num/den
}
#was spps_est_keep
spps.est.keep<-function(degs,dis,nguess,wstart=0,wsample=NULL,number.ss.iterations=5,number.ss.samples.per.iteration=2000,SS.infinity=0.04,
		   weight.type="Gile's SS",recID=NULL){
	degs[degs==0]<-1 
	mapping<-getestCstacked(degs,n=nguess,nit=number.ss.iterations,nsampsamp=number.ss.samples.per.iteration,trace=FALSE,SS.infinity=SS.infinity)
	weights=stats::approx(x=mapping$classes,y=1/mapping$probs,xout=degs,rule=2)$y
	pis <- 1/weights
#	if(weight.type == "Good-Fellows"){
#	  est=gf.est(deg=degs,outcome=dis,N=nguess,recID=recID,SS.infinity=SS.infinity)
#	}else{
	 if(is.factor(dis)){
		num <- rep(0,length=length(levels(dis)))
		a <- tapply(weights,as.numeric(dis),sum)
		num[as.numeric(names(a))] <- a
	 }else{
		num<-sum(dis*weights)
	 }
	 est<-num/sum(weights)
#	}
	
	# Note wts are actually the pis here
	list(samplewts=pis,
			est=est, 
			classes=mapping$classes, pvec=mapping$probs,
			mapping=mapping) 
}

#was ssps_est
spps.est<-function(degs,dis,nguess,wstart=0,wsample=NULL,nsampsamp=500,hajek=TRUE,mapping=NULL,SS.infinity=0.04,
		   weight.type="Gile's SS",recID=NULL){
	  degs[degs==0]<-1 #added 070708
	  if(is.null(mapping)){
		mapping<-getestCstacked(degs,n=nguess,nit=3,nsampsamp=nsampsamp,trace=FALSE,hajek=hajek,SS.infinity=SS.infinity)
	  }
	  pis=1/stats::approx(x=mapping$classes,y=1/mapping$probs,xout=degs,rule=2)$y
#	  if(weight.type == "Good-Fellows"){
#	    list(est=gf.est(deg=degs,outcome=dis,N=nguess,recID=recID,SS.infinity=SS.infinity),pis=pis)
#	  }else{
	    list(est=vh.est(pis,dis,wstart=wstart,wsample=wsample),pis=pis)
#	  }
}

######################################
#####  Main Function #################
######################################

# compute the bootstrap.  Returns the point estimate and bootstrap samples. We then take the sd of the bootstrap estimates as the estimate of the se (in simulations, this performs better than the quantiles of the bootstrap estimates)

#############  Arguments:
# degs = degrees of sample units
# n = population size,                         
# n0 = number of seeds, 
# dis is a nsamp vector with elements 1,...,g 
# refs.to.trait is a nsamp*g matrix, where refs.to.trait[i,j]=number of referrals from respondent i to recruits with z=j.
sppsboot4ppsreal<-function(degs,dis,n,n0,refs.to.trait,nit,                                   # ignore fullcup
		fullcup=rep(TRUE,length(degs)),
		number.ss.samples.per.iteration=500,
		fast=TRUE, useC=FALSE,
		weight.type="Gile's SS",
		recID=NULL,
		control=control.rds.estimates(),
		continuous=NA, is.cts=FALSE, is.quantile=FALSE,
		verbose=TRUE){ 
	fdis <- dis
	dis <- as.numeric(dis)
	outclasses <- levels(fdis)
	g=length(outclasses)
	nsamp<-length(degs)
	#
	# The next section pools small categories
	if(is.cts){
	        num.dis <- as.numeric(as.character(fdis))
	        oo <- order(order(as.numeric(outclasses)))
		olevels <- outclasses[order(oo)]
		nc0 <- as.character(oo[dis])
		nc <- mapply(gsub, list('\\d+'), sprintf('%.4d', as.numeric(regmatches(nc0, gregexpr('\\d+', nc0)))), nc0)
		tnc <- nc
		g <- sort(unique(tnc))
		tg <- g
		tt <- table(nc)
		ntt <- names(tt)
		sa <- tt[-1]+tt[-length(tt)]
		while(any(sa<10)){
		  i <- which.max(sa<10)
		  nc[nc==(ntt[i+1])] <- ntt[i]
	           g[ g==(ntt[i+1])] <- ntt[i]
		  tt <- table(nc)
		  ntt <- names(tt)
		  sa <- tt[-1]+tt[-length(tt)]
		}
	        mnc <- match(nc,tnc)
		gn=as.numeric(nc0)[match(g, tnc)]
		gm=match(gn,sort(unique(gn)))
		nc=gm[as.numeric(nc0)[mnc]]
		dis.agg <- tapply(num.dis,nc,mean)
#
		fdis.orig <- fdis
		dis.orig <- dis
	        outclasses <- as.character(dis.agg)
	        fdis <- factor(outclasses[nc],labels=outclasses)
	        dis <- nc
		g <- length(outclasses)
		refs.to.trait.orig <- refs.to.trait[,order(oo)]
		refs.to.trait <- matrix(0,ncol=g,nrow=nsamp)
	        for(i in 1:nsamp){
		  refs.to.trait[i,] <- tapply(refs.to.trait.orig[i,],gm,sum)
		}
	}	

	nrefs<-apply(refs.to.trait,1,sum)[which(fullcup)]
# nrefs[i] is the number of referrals for respondent i (i.e., the number of returned coupons)
	if(n*control$SS.infinity > nsamp){n <- round(nsamp/control$SS.infinity)}  
	sek<-spps.est.keep(degs,fdis,n,number.ss.samples.per.iteration=number.ss.samples.per.iteration,SS.infinity=control$SS.infinity,
		   weight.type=weight.type,recID=recID)
	# Note wts are actually the pis here
	wts<-sek$samplewts
	est<-sek$est         
	pvec<-sek$pvec
	if(fast){
		data.mapping<-sek$mapping
	}else{
		data.mapping<-NULL
	}
	classes<-sort(unique(degs))
	
	K=max(degs)
	classesbotha<-classes # Create unique labels for each (deg, dis) class
	for(i in 2:g)
		classesbotha<-c(classesbotha,classes+(i-1)*K)
	#g times the length of the vector of degree classes, to allow for degree-(g)disease categories. (We can assume the disease has g status ^_^)
	
	samplecountsa<-c(table(degs,dis))
	wtsbotha<-rep(pvec,g)
	
	CC=matrix(0,g,g)
	for(i in 1:nsamp){
		for(j in 1:g){
			CC[dis[i],j]=CC[dis[i],j]+refs.to.trait[i,j]
		}
	}
	# CC is a g*g matrix, where CC[i,j]=number of referrals from a recruiter with z=i to a recruit with z=j.
	
	deg=numeric()
	for(i in 1:g){
		deg[i]<-vh.est(wts[dis==i],degs[dis==i])
	}
	# deg[i] is the estimated population mean degree of nodes with z=i
	
	R = sweep(CC,1,apply(CC,1,sum),"/")
	R[is.na(R)] <- 0
	# R is a g*g matrix, where R[i,j]=observed proportional rates of referral of nodes with z=i for z=j.
	
	tiemtx = sweep(R,1,deg*est*n,"*")
	# est*n is the vector of counts of nodes with z=i in the population
	# deg*est*n is the vector of counts of the number ot ties
        #!!!! tiesnodes with z=i in the population
	# tiemtx is a matrix with the estimated number of ties (not referrals)
	# from class i to class j (in the population)
	
#	classesboth<-classesbotha[samplecountsa>0]
#	samplecounts<-samplecountsa[samplecountsa>0]   
#	wtsboth<-wtsbotha[samplecountsa>0]
	classesboth<-classesbotha
	samplecounts<-samplecountsa   
	wtsboth<-wtsbotha
	
	manynewests<-matrix(0,nrow=nit,ncol=g)
	manynewnms<-rep(0,nit)
	manysamp<-vector("list",length=nit)
	
	if(weight.type %in% c("Gile's SS","Good-Fellows")){
		effn <- n
	}else{
		effn <- 10^8
	}
	
# moved
	newprops<-probtodist(classesboth,samplecounts, wtsboth, n)$props
	props2<-newprops/sum(newprops)
	nbyclass<-round(n*props2)  
#	nbyclass[nbyclass==0]<-1
	offby<-n-sum(nbyclass)
	if(is.na(offby)){print("Error in getincl: offby"); } 
 	tempties <- (tiemtx+t(tiemtx))/2
#	symetrize the ties so that
	# tempties is a matrix with the estimated number of ties (not referrals)
	# from class i to class j (in the population). It is based on extrapolating
	# the number of referrals 
# moved
	pis=1/stats::approx(x=data.mapping$classes,y=1/data.mapping$probs,xout=rep(classes,g),rule=2)$y
	data<-list(n=n, effn=effn,
			classesboth=classesboth,
			n0=n0, nsamp=nsamp, nrefs=nrefs,
			mapping=data.mapping, K=K, g=g, pis=pis, nit=nit,
			control=control,
			number.ss.samples.per.iteration=number.ss.samples.per.iteration,
			props2=props2,offby=offby,nbyclass=nbyclass,tempties=tempties,
			weight.type=weight.type
                  )
	
	if(useC){
	  Cret <- .C("bsC",
			nbyclass=as.integer(data$nbyclass),
			classesboth=as.integer(data$classesboth),
			nrefs=as.integer(data$nrefs),
			props2=as.double(data$props2),
			tempties=as.double(data$tempties),
			pis=as.double(data$pis),
			est=as.double(rep(0,data$g*data$nit)),
			nm=as.double(rep(0,data$nit)),
			numsamp=as.integer(data$nit),
			offby=as.integer(data$offby),
			K=as.integer(data$K),
			nc=as.integer(length(data$nbyclass)),
			g=as.integer(data$g),
			N=as.integer(data$n),
			n=as.integer(data$nsamp),
			n0=as.integer(data$n0),
			PACKAGE="RDS")
          manynewests <- matrix(Cret$est,ncol=g,byrow=TRUE)
          colnames(manynewests) <- outclasses
          manynewnms <- Cret$nm
	}else{
	  bsfn <- function(i, data){
                offby <- data$offby
                nbyclass <- data$nbyclass
		
		if(offby>0){
			for(ii in 1:offby){
				dif<-data$props2-nbyclass/data$n
				dif[dif<0]<-0
				tochange<-sample(1:length(dif),1,prob=dif)
				nbyclass[tochange]<-nbyclass[tochange]+1
			}
		}else{if(offby<0){
				for(ii in 1:abs(offby)){                
					dif<-nbyclass/data$n - data$props2
					dif[dif<0]<-0
					dif[nbyclass==1]<-0 
					if(sum(dif==0)){
						dif<-1-(data$props2-nbyclass/data$n)
						dif[nbyclass==1]<-0
					        dif[dif<0]<-0
					}       
					tochange<-sample(1:length(dif),1,prob=dif)
					while(any(nbyclass[tochange]<1)){
					 tochange<-sample(1:length(dif),1,prob=dif)
						print(nbyclass)
						print(tochange)
						print(dif)
					}
					nbyclass[tochange]<-nbyclass[tochange]-1
				}
			}}   
		
#		popclass<-rep(data$classesboth,times=nbyclass)    # Create a population of labels (for each (deg, dis) class)
		popclass<-data$classesboth
		nc <- length(popclass)
		K <- data$K
		idis<-((popclass-1)%/%K)+1  
		ideg<-round(popclass-K*(popclass%/%K))
		ideg[ideg==0]<-K
		is.seed <- rep(FALSE,nsamp)
		recID <- rep(NA,data$nsamp)
		dissample <- rep(0,data$nsamp)
		degsample <- dissample
		
#		newdis<-data$newdis
#		newdeg<-data$newdeg
 		tempties<-data$tempties
#		nsample0<-sample.int(data$n,size=data$n0,prob=abs(newdeg),replace=FALSE)
		pclass=ideg*nbyclass
                for(i in 1:data$n0){
		 j<-sample.int(nc,size=1,prob=ideg*nbyclass)
		 pclass[j]<-pclass[j]-ideg[j]
                 degsample[i] <- ideg[j]
                 dissample[i] <- idis[j]
                 is.seed[i] <- TRUE
		}
#		todisnew=matrix(0,nrow=data$nsamp,ncol=g)
		
#		tempties <- (data$tiemtx+t(data$tiemtx))/2
 		
# 		nprop <- tabulate(newdis,nbins=g)
 		nprop <- tapply(nbyclass,idis)
 		nprop <- nprop / sum(nprop)
 		for(k in 1:g){
 			if(all(tempties[k,]==0)){
 				tempties[k,] <- nprop
 			}
 			if(all(tempties[,k]==0)){
 				tempties[,k] <- nprop
 			}
 		}
		
		nrefs<-data$nrefs 
                crefs <- which(nrefs>0)
                ncrefs <- length(crefs)
		crefs.rs<-nrefs[crefs[sample.int(ncrefs, nsamp, replace=TRUE)]]
		icrefs<-1
		numberfrom<-crefs.rs[icrefs] 

#		nsample<-c(nsample0,rep(0,(data$nsamp-data$n0)))    ### Vector W !!!
		activnode<-1
		countrefs<-0
#		temptiess <- apply(tempties,1,sum)
		for(mm in (data$n0+1):data$nsamp){                 ### loop begins!!!  mm is not m in the paper!
			thisdis = dissample[activnode]
#			thisdis = newdis[nsample[activnode]]
#			p <- tempties[thisdis,] / temptiess[thisdis]
#        ## p is the vector of probability

		        nextdis<-sample.int(g, 1, replace=FALSE, prob=tempties[thisdis,])
			
#       numsfrom are the indices of those which have not been sampled
#       who have the same disease status as "nextdis".
#			nprob<-pnewdeg*(newdis==nextdis)
			nprob<-pclass*(idis==nextdis)
			
			if(all(nprob==0)){
			       nprob<-pclass
			       warning(paste("Ran out of",nextdis))
			}
			if(sum(nprob>0)==1){
				nextresp<-which(nprob > 0)
			}else{     
				# Draw a node from those with newdis left proportional to degree
		                nextresp<-sample.int(length(nprob), 1, replace=FALSE, prob=nprob)
			}
			
#			totaltmp<-sum(newdeg[numsfrom])  # the sum of degrees left with newdis
			totaltmp<-sum(nprob)  # the sum of degrees left with newdis
#                       pnewdeg[nextresp] <- 0
		 	pclass[nextresp]<-pclass[nextresp]-ideg[nextresp]
			
			nextdis <- idis[nextresp]
#			todisnew[activnode,nextdis]=todisnew[activnode,nextdis]+1
			tdeg<-ideg[nextresp] # the degree of the sample
			tempties[,nextdis]<-tempties[,nextdis]*(totaltmp-tdeg)/totaltmp       ###??? Shouldn't tempties be symmetric?

#			nsample[mm]<-nextresp
                 	degsample[mm] <- tdeg
                 	dissample[mm] <- nextdis
                 	recID[mm] <- activnode
			countrefs<-countrefs+1
# numberfrom is the number of recruits to get for the current recruiter
			if((mm<nsamp)&(countrefs==numberfrom)){
				activnode<-activnode+1                     ### move to the next seed (or node)!!! i=i+1
				
				countrefs<-0 
		                icrefs<-icrefs+1
				numberfrom<-crefs.rs[icrefs] 
			}
		}                                            ### loop ends!!!
		
#		degsample<-newdeg[nsample]
#		dissample<-factor(newdis[nsample],
		dissample<-factor(dissample,
				levels=1:length(outclasses),labels=outclasses,exclude=NULL)
#		manysamp<-list(deg=degsample,dis=dissample,refs.to.trait=todisnew)
#THE ORIGINAL ONE:
		manynewests<-spps.est(degsample,
				dissample,
				data$effn,nsampsamp=ceiling(number.ss.samples.per.iteration/4), 
				mapping=data$mapping,
				SS.infinity=control$SS.infinity,
				weight.type=data$weight.type,recID=recID)
  
#		list(samp=manysamp, newests=manynewests)
		list(newest=manynewests$est,
		     newnm=sum(manynewests$pis,na.rm=TRUE)^2/sum(manynewests$pis^2,na.rm=TRUE))
		
	  }
	  if(verbose)  cat(paste('Computation 0% completed ...\n',sep=""))
	  for(i in 1:nit){
		out <- bsfn(i,data)
#		manysamp[[i]] <- out[["samp"]]
		manynewests[i,] <- out[["newest"]]
		manynewnms[i] <- out[["newnm"]]
		
		if(verbose) {
			if(i == trunc(i/(nit/10))*(nit/10)){
				cat(paste((100/10)*trunc(i/(nit/10)),'% completed ...\n',sep=""))
			}
		}
	  }
	}

	if(is.cts){
	  # Note sek$samplewts are actually the pis here
	  wts.agg <- tapply(1/sek$samplewts,dis,sum)
	  fn.mu <- function(x,dis){wtd.mean(dis,x,na.rm=TRUE,normwt=TRUE)}
	  if(is.quantile){ 
	    fn <- function(x,dis){wtd.quantile(dis,x,probs=continuous,na.rm=TRUE,normwt=TRUE)}
	    list(point_estimate=fn(wts.agg,dis.agg),
	       bsests=matrix(apply(manynewests,1,fn,dis.agg),ncol=1),bsnm=manynewnms,
	       se_estimate=stats::sd(apply(manynewests,1,fn,dis.agg)),is.cts=is.cts,is.quantile=is.quantile,
	       mu=fn.mu(wts.agg,dis.agg),
	       sigma2=fn.mu(wts.agg,dis.agg^2)-fn.mu(wts.agg,dis.agg)^2)
	  }else{
	    list(point_estimate=fn.mu(wts.agg,dis.agg),
	       bsests=matrix(apply(manynewests,1,fn.mu,dis.agg),ncol=1),bsnm=manynewnms,
	       se_estimate=stats::sd(apply(manynewests,1,fn.mu,dis.agg)),is.cts=is.cts,is.quantile=is.quantile,
	       mu=fn.mu(wts.agg,dis.agg),
	       sigma2=fn.mu(wts.agg,dis.agg^2)-fn.mu(wts.agg,dis.agg)^2)
	  }
	}else{
	  list(point_estimate=sek$est,bsests=manynewests,bsnm=manynewnms,
		se_estimate=apply(manynewests,2,sd),is.cts=is.cts,is.quantile=is.quantile)
	}
}
