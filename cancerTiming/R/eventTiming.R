#for now, history is given as matrix 
#(1/ncopies, ..., 1); Must be length equal to the number of copies of chromosome for the history
#later will have function that takes encoding of history and creates this vector
#if exactAllele=TRUE, then assume that x/m is exactly the true allele frequency P -- not implemented yet
#... corresponds to arguments passed to the internal fitting functions (e.g. Bayes method of approximation)
###Should code some EM steps to get good init point...
##Version 1.0.2 add ids for mutations; made returnAssignments and returnData the same thing -- one data.frame
eventTiming<-function(x, m, history, totalCopy, method=c("fullMLE","partialMLE","Bayes"),type=c("gain","CNLOH"),
seqError=0, bootstrapCI=NULL, B=if(method=="Bayes") 10000 else 500,CILevel=0.95, normCont=0,verbose=TRUE,
returnAssignments=FALSE,coverageCutoff=1,minMutations=10,init=NULL,maxiter=100, tol=0.0001,mutationId=1:length(x),...){
	method<-match.arg(method)
	doCI<-!is.null(bootstrapCI)
	type<-match.arg(type)

	nMuts<-length(x)
	if(length(m)!=nMuts) stop("x and m must be of same length")
	if(length(mutationId)!=length(unique(mutationId)) & verbose) warning("mutationId not unique values; using index as id instead.")
	if(length(mutationId)!=nMuts & verbose) warning("mutationId an invalid length; using index as id instead.")
	
	if(returnAssignments) returnData<-TRUE 
	else returnData<-FALSE
	
	A<-history #this will be a function later
	nCopies<-totalCopy
	if(is.null(dim(A))) stop("'history' should be a matrix of size (nEvents +1) x (nEvents +1)")
	if(ncol(A)!=nrow(A)) stop("'history' should be a square matrix (hint: do not include the allele frequency 1 except for CNLOH events)")
	K<-ncol(A)-1
	if(length(normCont)!=1 || normCont>1 || normCont<0) stop("normCont must be given and should be single value between 0 and 1.")

	#get final allele frequencies
	possAlleles<-allAF(nCopies,normCont=normCont,type="mutation")[[1]] 
	if(type=="gain") possAlleles<-head(possAlleles,-1) #can't have the allele frequency 1
	if(nrow(A)!=length(possAlleles)) stop(length(possAlleles),"possible alleles for this scenario, but history only has",nrow(A),"rows.")
	possAlleles<-errorAF(possAlleles,seqError=seqError)

	if(coverageCutoff<1) stop("coverageCutoff must be at least 1")
	nBelowCutoff<-length(which(m<coverageCutoff))
	nZero<-length(which(m==0))
	if(any(m<coverageCutoff)){		
		if(verbose) warning(nBelowCutoff, "mutations present below cutoff",coverageCutoff,", will be ignored")
		if(verbose & nZero>0) warning(nZero, "mutations have no reads at all")
		x<-x[which(m>=coverageCutoff)]
		mutationId<-mutationId[which(m>=coverageCutoff)]
		m<-m[which(m>=coverageCutoff)]
	}
	summaryTable<-cbind(nMuts,length(m),nZero,nBelowCutoff)
	colnames(summaryTable)<-c("Total Mutations","Qualify for fitting","Zero Read Coverage","Coverage Below Cutoff")
	summaryTable<-t(summaryTable)
	colnames(summaryTable)<-"Number"
	call<-list(alleleSet=possAlleles,
		history=history, 
		totalCopy=totalCopy,
		type=type,
		method=method,
		seqError=seqError,
		normCont=normCont,
		coverageCutoff=coverageCutoff,
		minMutations=minMutations,
		B=B,
		init=init,
		maxiter=maxiter, 
		tol=tol)		
	if(returnData) inputData=data.frame(mutationId=mutationId,x=x,m=m)
	
	failReason<-NULL
	success<-TRUE
	if(length(m)<minMutations){
		failReason<-paste(failReason,"not enough mutations above cutoff that satisfy coverage criteria",sep=",")
		if(verbose) warning("less than ",minMutations," mutated locations meet coverage criteria; no estimate of pi will be calculated")
		success<-FALSE
	}
	# if(all(x==m)){
	# 	failReason<-paste(failReason,"all x equal to m",sep=",")
	# 	success<-FALSE
	# }
	dummyOutput<-list("pi"=rep(NA,length=ncol(A)),alleleSet=possAlleles,optimDetails=list(nIter=NA,initial=NA,pathOfQ=NA),perLocationProb=NA,q=NA)
	names(dummyOutput$pi)<-paste("Stage",0:K,sep="")
	ci<-matrix(rep(NA,length=2*length(dummyOutput$pi)),ncol=2)
	row.names(ci)<-names(dummyOutput$pi)
	dummyOutput$piCI<-ci
	if(!is.null(failReason)){
		output<-dummyOutput
	}
	else{	
		rankA<-qr(A)$rank
		piId <- rankA==ncol(A)
		piEst<-rep(NA,ncol(A))
		if(piId){
			if(method%in%c("fullMLE","Bayes")){		
				### Estimate the q vector:
				output<-try(.estimateQ(x,m,possAlleles,history=history,init=init,maxiter=maxiter, tol=tol))
				if(output$optimDetails$nIter==maxiter){
					failReason<-paste(failReason,"hit maximum number of iterations",sep=",")
					success<-FALSE
				}
			}
			if(method=="partialMLE" ){
				mleAlleles<-mleAF(x,m,totalCopy=nCopies,seqError=seqError,normCont=normCont,maxCopy=if(type=="gain") totalCopy-1 else totalCopy)
				#q<-mleAlleles$alleleSet$Freq/sum(mleAlleles$alleleSet$Freq)
				output<-try(.estimateQ(x,m,possAlleles,alleleFreq=mleAlleles$alleleSet$Freq,history=history,init=init,maxiter=maxiter, tol=tol))				
#				output<-list(q=q,perLocationProb=mleAlleles$perLocationProb,alleleSet=possAlleles,optimDetails=list(nIter=NA,initial=NA,pathOfQ=NA))
			}
			if(!inherits(output, "try-error")){	
				piEst<-solve(A)%*%output$q
				piEst<-as.vector(piEst)
				piEst<-piEst/sum(piEst)
				if(method=="Bayes"){
					initBayes<-piEst
					#initBayes<-seq(1,length(possAlleles))/length(possAlleles)
					bayesout<-.bayesEstimate(x,m,alleleSet=possAlleles,history=A,init=initBayes,nsim=B,CILevel=CILevel,...)
					piEst<-bayesout$pi
					q<-A%*%piEst
					q<-q/sum(q)
					output<-list(piCI=bayesout$piCI,q=q,perLocationProb=NA,alleleSet=possAlleles,optimDetails=list(nIter=NA,initial=rbind(initBayes,bayesout$mode),pathOfQ=NA))					
				}
			}
			else{
				if(method=="Bayes"){
					initBayes<-seq(1,length(possAlleles))/length(possAlleles)
					bayesout<-.bayesEstimate(x,m,alleleSet=possAlleles,history=A,init=initBayes,nsim=B,...)
					piEst<-bayesout$pi
					q<-A%*%piEst
					q<-q/sum(q)
					output<-list(piCI=bayesout$piCI,q=q,perLocationProb=NA,alleleSet=possAlleles,optimDetails=list(nIter=NA,initial=rbind(initBayes,bayesout$mode),pathOfQ=NA))					
					
				}
				else{
					output<-dummyOutput
					failReason<-paste(failReason,"error returned in estimating q",sep=",")
					success<-FALSE	
				}
				
			}
			if(method=="Bayes") row.names(output$optimDetails$initial)<-c("mle","postMode")
		}
		else{
			output<-dummyOutput
			failReason<-paste(failReason,paste("history matrix not invertible, rank=",rankA,sep=""),sep=",")
			success<-FALSE
			if(verbose) warning(paste("history matrix not invertible, rank=",rankA,sep=""))
		}
		names(piEst)<-paste("Stage",0:K,sep="")
		output$pi<-piEst
		
		if(!returnAssignments) output<-output[-grep("perLocationProb",names(output))]
		else{ output[["perLocationProb"]]<-data.frame(inputData,output[["perLocationProb"]],check.names=FALSE)}
		
		#do bootstrap CI, or format the Bayesian ones; makes sure the elements returned are in right order
		if(doCI & method!="Bayes"){
			if(!any(is.na(output$pi))){
				bootsample<-bootstrapEventTiming(call=call,B=B,x=x,m=m,type=bootstrapCI,pi=output$pi)
				ci<-t(apply(bootsample,2,quantile,c((1-CILevel)/2,1-(1-CILevel)/2)))
			}
			else {
				ci<-matrix(rep(NA,length=2*length(output$pi)),ncol=2)
				row.names(ci)<-names(output$pi)
			}
			#colnames(ci)<-c("lower","upper")
			output$piCI<-ci
			output<-output[c("pi","piCI","q","perLocationProb","optimDetails")]
		}
		else{
			if(method!="Bayes") output<-output[c("pi","q","perLocationProb","optimDetails")]
			else{
				if(!any(is.na(output$pi))){ row.names(output$piCI)<-names(output$pi)}
				output<-output[c("pi","piCI","q","perLocationProb","optimDetails")]
			}
		}
	}
	return(c(output,list(summaryTable=summaryTable,success=success,failReason=if(!is.null(failReason)) gsub("^,","",failReason) else failReason,call=call)))
}

#simple EM algorithm to calculate \hat{q}, the proportions of the multinomial, constrained
#also gives probabilities of assignments to different alleles.
#if alleleFreq not null, then means know the P[i], and alleleFreq is the number of P[i] equal to each allele
#xGreaterZero determines whether the likelihood will account for the fact that only observe X[i]>0
#useGradient determines whether to use the gradient I calculated, or set gr=NULL
.estimateQ<-function(x,m,alleleSet,alleleFreq=NULL,history,init=NULL,maxiter=100, tol=0.0001,xGreaterZero=TRUE, useGradient=TRUE){
	N<-length(x)
	possAlleles<-alleleSet
	lengthOfTheta<-length(possAlleles)-1
	if(nrow(history)!=length(possAlleles)) stop("length of alleleSet must equal the number of rows of history.")
	##check if it is one of the easy ones:
	Avec<-as.vector(history)
	easyType<-"no"
	if(identical(Avec,c(0, 1, 2, 0))) easyType<-"CNLOH"
	if(identical(Avec,c(1 ,1, 3, 0))) easyType<-"SingleGain"
	
	Ainv<-solve(history)
	if(is.null(init)){
		init<-history%*%rep(1/length(possAlleles),length=length(possAlleles))	
		init<-as.vector(init)/sum(init)
#		init<-c(.2,.8) #just for testing
	} 
	##Break it into 2 because useful to have a function that assigns most likely allele to each
	EStep<-function(q){ ####Calculate P(Pi=a | Xi,q) each i and a
		#P(Pi=a | Xi)=P(Xi|Pi=a)P(Pi=a)/sum(P(Xi|Pi=a)P(Pi=a))
		
		#calculate the log of the probability (numerator)
		logPPerAllele<-mapply(possAlleles,q,FUN=function(a,qa){
			dbinom(x=x, size=m, prob=a,log=TRUE)+log(qa) 
		}) #returns a N x numb.Alleles Matrix; needs to be normalized
		
		#-----
		#calculate log of sum of probabilities (denominator)
		#-----
		#note if individual probabilities very small, then denominator sums to zero, so divide by zero -- problem!
		#use log sum exponential trick:
		# log(e[theta1]+e[theta2]+...)=m+log(e[theta1-m]+e[theta2-m]+...), where m=max(theta[i])
		#here theta1=log[P(Xi|Pi=a)]+log[P(Pi=a)]
		#-----

		rowMax<-apply(logPPerAllele,1,max) #max of theta[i]
		logPPerAlleleMinusMax<-sweep(logPPerAllele,1,rowMax,"-") #subtract off m from each row
		logdenom<-rowMax+log(apply(exp(logPPerAlleleMinusMax),1,sum))
		
		#-----
		#calculate log of sum of probabilities (denominator)
		#-----
		#now overall prob=exp[num]/exp[denom]=exp[num-denom]
		#-----
		logPPerAllele<-sweep(logPPerAllele,1,logdenom,"-")
		
	
 		# #-----
 		# #Old simplistic code:
 		# #-----
 		#     	PPerAllele<-mapply(possAlleles,q,FUN=function(a,qa){
 		# 	dbinom(x=x, size=m, prob=a,log=FALSE)*(qa) 
 		# }) #returns a N x numb.Alleles Matrix; needs to be normalized
 		# PPerAllele<- prop.table(PPerAllele,1) #divide by the sum of each row
 		
		#compare the two:
		# cbind(exp(logPPerAllele),PPerAllele)
 		
		
		return(exp(logPPerAllele))
	}
	
	###Constraints for the q
	#The feasible region is defined by ui %*% theta - ci >= 0. 
	if(easyType=="no"){
		uiL1<-diag(rep(-1,lengthOfTheta),nrow=lengthOfTheta,ncol=lengthOfTheta)
		ciL1<-rep(1,lengthOfTheta)
		uiGr0<-diag(rep(1,lengthOfTheta),nrow=lengthOfTheta,ncol=lengthOfTheta)
		ciGr0<-rep(0,lengthOfTheta)
		uiA<-Ainv%*%(rbind(diag(lengthOfTheta),rep(-1,lengthOfTheta)))
		ciA<- - Ainv%*%c(rep(0,lengthOfTheta),1)
		ui<-rbind(uiL1,uiGr0,uiA)
		ci<- c(-ciL1,-ciGr0,ciA)		
	}
	
	#make matrix of (1-a[j])^m[i]
	MStep<-function(PMat=NULL,qinit){ #estimate of q constrained so that solve(history)%*%q is all positive
		if(is.null(alleleFreq)) Na<-colSums(PMat) #from E-Step, the estimated values from the multinomial...basically just the maximum likelihood estimation if knew Pi exactly
		else Na<-alleleFreq
		if(easyType %in% c("CNLOH","SingleGain")){
			#print("easytype")
			if(easyType=="CNLOH"){ 
				a<-0
				b<-1
			}
			if(easyType=="SingleGain"){
				a<-1/2
				b<-1
			}
			qout<-Na[1]/sum(Na)
			if(qout>b) qout<-b
			if(qout<a) qout<-a
			qout<-c(qout,1-qout)
		}
		else{
			f<-function(q){
				qS<-1-sum(q)
				#-dmultinom(Na, prob=c(q,qS), log = TRUE) #just to check; should be equivalent, but faster to not calculated normalizing constant
				ll<- sum(Na*log(c(q,qS)))
				if(xGreaterZero){
					#missZero[i]= 1- sum_j{(1-a[j])^m[i] q[j]}
					missZero<-1-rowSums(mapply(possAlleles,c(q,qS),FUN=function(a,qa){
						(1-a)^m*qa 
					})#returns a N x numb.Alleles Matrix; needs to be normalized
					) #sum over all j, so missZero is vector of length i single value for each i
					ll<-ll-sum(missZero)
				} 
				return( - ll) #because minimizing
			}

			#used in the gradient
			#calculate (1-a[S])^{m[i]} - (1-a[j])^{m[i]} 
			#returns a N x S-1 Matrix; needs to be normalized
			aS<-tail(possAlleles,1)
			# #bug reported
			# num<-mapply(head(possAlleles,-1),q,FUN=function(a,qa){
			# 	(1-aS)^m - (1-a)^m
			# })
			#fix with
			num<-sapply(head(possAlleles,-1),FUN=function(a,qa){
				(1-aS)^m - (1-a)^m
			})
			gr<-function(q){
				qS<-1-sum(q)
				NS<-tail(Na,1)
				g<-(head(Na,-1)/q - NS/qS) #was - (head(Na,-1)/q - NS/qS*q) #which is error; why got wrong answer.
				if(xGreaterZero){
					#missZero[i]= 1- sum_j{(1-a[j])^m[i] q[j]}
					missZero<-1-rowSums(mapply(possAlleles,c(q,qS),FUN=function(a,qa){
						(1-a)^m*qa 
					})#returns a N x numb.Alleles Matrix; needs to be normalized
					) #sum over all j, so missZero is vector of length i single value for each i
					g<-g-colSums(sweep(num,1,missZero,"/"))
				}
				return(-g) #because minimizing
			}
			theta<-head(qinit,-1)
			if(lengthOfTheta==1){
				upper<-min((ci/ui)[ui<0])
				lower<-max((ci/ui)[ui>0])
				out<-optimize(f=f,interval=c(lower,upper),tol=.Machine$double.eps)
				qout<-c(out$minimum,1-out$minimum)
			}
			else{
				out<-constrOptim(theta=theta,f = f, grad=if(useGradient) gr else NULL, ui=ui, ci=ci) #when use gr, get wrong answer -- returns me to init, doesn't iterate, etc. WHY????
				qout<-c(out$par,1-sum(out$par))
				}
		}
		names(qout)<-names(possAlleles)
		return(qout)	
	}
	#print("gr uses - (head(Na,-1)/q - NS/qS) ")
	TotalStep<-function(q){
	 	MStep(EStep(q),qinit=q)
	}
	if(is.null(alleleFreq)){
		qOld <- init;
		qNew<-TotalStep(qOld)
		allQ<-cbind(qOld,qNew)
		nIter<-1
		#convMessages<-c(firstOut$convergence)
		while(sum(abs(qNew - qOld)) >= tol & nIter <= maxiter){
			qOld<-qNew
			pMat<-EStep(qOld)
			qNew<-MStep(pMat,qinit=qOld)
			#qNew<-c(mstepOut$par,1-sum(mstepOut$par))
			#convMessages<-c(convMessages,mstepOut$convergence)
			nIter <- nIter + 1
			allQ<-cbind(allQ,qNew)
			#if(nIter %% 20 ==0) print(sprintf("finished iteration %s", nIter))
		}
		pMat<-EStep(qNew)
	}
	else{
		qNew<-MStep(qinit=init)
		pMat<-NA
		nIter<-NA
		allQ<-NA
	}
	names(qNew)<-names(possAlleles)
	return(list(q=qNew,perLocationProb=pMat,alleleSet=possAlleles,optimDetails=list(nIter=nIter,initial=init,pathOfQ=allQ)))
}

# > x$pi
#    Stage0    Stage1 
# 0.1096573 0.8903427 
# > x$q
#        1/2        2/2 
# 0.96206078 0.03793922 

#the Learn Bayes does weird things with the error
#I copy them here to try to fix it; If I comment out the warn options, then will get wanring about optim not appropriate for 1-dim.
.laplace<-function (logpost, mode, ...) 
{
	currWarn<-options()$warn
	#print(currWarn)
    options(warn = -1)
    fit = optim(mode, logpost, gr = NULL, ..., hessian = TRUE, 
        control = list(fnscale = -1))
    #options(warn = 0)
	options(warn = currWarn)
    mode = fit$par
    h = -solve(fit$hessian)
    p = length(mode)
    int = p/2 * log(2 * pi) + 0.5 * log(det(h)) + logpost(mode, ...)
    stuff = list(mode = mode, var = h, int = int, converge = fit$convergence == 0)
    return(stuff)
}
.sir<-function (logf, tpar, n) 
{
    k = length(tpar$m)
    theta = LearnBayes::rmt(n, mean = c(tpar$m), S = tpar$var, df = tpar$df)
    lf = matrix(0, c(dim(theta)[1], 1)) #really just n x 1 bc dim(theta)[1]=n
    for (j in 1:dim(theta)[1]) lf[j] = logf(theta[j, ])
    lp = LearnBayes::dmt(theta, mean = c(tpar$m), S = tpar$var, df = tpar$df, 
        log = TRUE)
    md = max(lf - lp)
    wt = exp(lf - lp - md)
    probs = wt/sum(wt)
    indices = sample(1:n, size = n, prob = probs, replace = TRUE)
    if (k > 1) 
        theta = theta[indices, ]
    else theta = theta[indices]
    return(theta)
}
.rejectsampling<-function (logf, tpar, dmax, n) 
{
    d = length(tpar$m)
    theta = LearnBayes::rmt(n, mean = c(tpar$m), S = tpar$var, df = tpar$df)
    lf = matrix(0, c(dim(theta)[1], 1))
    for (j in 1:dim(theta)[1]) lf[j] = logf(theta[j, ])
    lg = LearnBayes::dmt(theta, mean = c(tpar$m), S = tpar$var, df = tpar$df, 
        log = TRUE)
    if (d == 1) {
        prob = exp(c(lf) - lg - dmax)
        return(theta[runif(n) < prob])
    }
    else {
        prob = exp(lf - lg - dmax)
        return(theta[runif(n) < prob, ])
    }
}

# data(mutData)
# ACNLOH<-matrix(c(1,3,1,0),ncol=2,nrow=2,byrow=TRUE)
# onlyMuts<-subset(mutData,is.na(rsID) & position <= 1.8E7)
# onlyMuts$t_depth<-onlyMuts$t_ref_count+onlyMuts$t_alt_count
# x<-eventTiming(x=onlyMuts$t_alt_count,m=onlyMuts$t_depth,history=ACNLOH,totalCopy=2,type="CNLOH",normCont=0.22)

#outnew<-apply(combSeq,1,function(th){.logPosterior(th,x=xout$nMutAllele,m=xout$nReads,alleleSet= eventOut$call$alleleSet,history=AMat,jacob="new")})
# > length(outnew)
# [1] 10000
# > setwd('/Users/epurdom/Documents/Sequencing/CancerTiming/Simulations')
# > outold<-sapply(xseq,function(th){.logPosterior(th,x=xout$nMutAllele,m=xout$nReads,alleleSet= eventOut$call$alleleSet,history=AMat,jacob="old")})

# xseq<-seq(-10,10,length=100)
# outnew<-sapply(xseq,function(th){.logPosterior(th,x=onlyMuts$t_alt_count,m=onlyMuts$t_depth,alleleSet= x$call$alleleSet,history=ACNLOH,jacob="new")})
# outold<-sapply(xseq,function(th){.logPosterior(th,x=onlyMuts$t_alt_count,m=onlyMuts$t_depth,alleleSet= x$call$alleleSet,history=ACNLOH,jacob="old")})


#run into problems when get close to the boundary -- 
#conversion between theta and pi is tenuous
#replace log(c+sum(a*exp(theta))) with the following
.mylogSumExp<-function(x,c=1,a=1){
	if(length(a)==1) a<-rep(a,length(x))
	# if(c==0 || any(exp(x)>1e10)){ 
	# 	#ignore the +c term in this case and just calculate log( sum_i a_i exp(x_i))
	# 	#then log( sum_i a_iexp(x_i))=m + log(sum_i a_i exp(x_i-m)), where m = max(x_i) http://lingpipe-blog.com/2009/06/25/log-sum-of-exponentials/
	# }
	#else	log(c+sum(a*exp(x))) #this can always be written as c=0 if add term to the x...
	if(c!=0){
		x<-c(x,log(c))
		a<-c(a,1)
	}
	m<-max(x)
	return(m+log(sum(a*exp(x-m))))			

}
#always return on the theta scale, but can calculate either log density of theta or pi, depending on choice of value of 'parameter'
.logPosterior<-function(theta,x, m, alleleSet, history, alpha = 1,parameter=c("theta","pi")){
	parameter<-match.arg(parameter)
	K<-length(theta) 	
	if(length(alleleSet)!=K+1) stop("invalid values for theta or alleleSet")


	##############
	##Parts of calculating the posterior
    pOfPi <- function(theta) {
#        (alpha - 1) * (sum(theta - log(1 + sum(exp(theta))))) #change 9/26, but should be equivalent...
        (alpha - 1) * (sum(theta) - (K+1)*.mylogSumExp(theta) ) 
    }
 	pOfX<-function(theta){ #changed 9/30/2013 so that can handle large theta:
  		PPerX<-mapply(x,m,FUN=function(dat,size){
  			p<-dbinom(x=dat, size=size, prob=alleleSet) #vector of probabilities (per allele) of seeing data [i] with that allele
  			#p'A
  			Avec<-apply(history,2,function(x){sum(x*p)}) #avec[j]=sum_a P_a(X) A[a,j] for each stage j
  			.mylogSumExp(c(theta,log(1)),c=0,a=Avec) #	log( sum_j  avec[j]exp(theta[j]) + avec[K] )
  		}) 
  		if(length(PPerX)!=length(x)) stop("coding error")

  		#Nlog(sum(S[j]*v[j]))
  		 S<-2+0:(K) # this is sequential amplification; could adjust here to be for more complicated
  		denomTerm<-length(x)*.mylogSumExp(c(theta,log(1)),c=0,a=S) #log( sum_j (2+j)*exp(theta_j) + (2+K))

   		return(sum(PPerX)-denomTerm)
   	}
 	chgOfVar<-function(theta){ #new version 9/26/2013; doesn't affect unless dim(theta)>1
 		K<-length(theta)
 		#sum(theta)-(K+1)*log(1+sum(exp(theta)))
 		sum(theta)-(K+1)*.mylogSumExp(theta)
 	}
	if(parameter=="theta") return(pOfPi(theta) + chgOfVar(theta) + pOfX(theta))
	else return(pOfPi(theta) + pOfX(theta))
}


.bayesEstimate<-function (x, m, alleleSet, history, alpha = 1, bayesApproxMethod = c("sir", "inv", "rej","exact"), init, tdf = 4, nsim = 10000, CILevel = 0.95) 
{
	#change of variables so that now theta[j]=log(pi[j-1]/pi[K])
	#assume alpha a single scalar
	K<-length(init)-1
    method <- match.arg(bayesApproxMethod)
	if(method=="inv" & K>1) stop("direct inversion method is only implemented for K=1")
    if (!method %in% c("sir","inv")) 
        stop("Other methods than 'sir' are not yet operational")
   # require(LearnBayes) #copied the functions into the code so I could have better control of them, but still need some helper functions...
    pi2theta <- function(piV) {
        head(log(piV/tail(piV, 1)), -1)
    }
    theta2pi <- function(theta, returnFullPi=TRUE) {
		logPiV<-theta-.mylogSumExp(theta,c=1)
		piV<-exp(logPiV)
        if(returnFullPi) return(c(piV, 1 - sum(piV)))
		else return(piV)
    }
    if (method == "inv") { #if 1-dimensional, can approximate it easily with a grid of the possibilities
		logPiPosterior<-function(theta){.logPosterior(theta,x=x,m=m,alleleSet=alleleSet,history=history,alpha=alpha,parameter="pi")}
		piVec<-seq(1e-5,1-1e-5,length=1000)
		piVec<-cbind(piVec,1-piVec)
		thetaVec<-apply(piVec,1,pi2theta)
		#store the values on the log scale so as to make it numerically stable.
		logpdensity<- sapply(thetaVec,logPiPosterior)
		logNormFac<-.mylogSumExp(logpdensity,c=0,a=1)
		logpdensity<-logpdensity-logNormFac
		logcumSum<-sapply(1:length(logpdensity),function(ii){
			.mylogSumExp(logpdensity[1:ii],c=0)
		})
		piV<-c(piVec[,1],1)
		cumSum<-c(0,exp(logcumSum),1)
		#cdf<-stepfun(piV,cumSum,right=TRUE) #f(x[i])=y[i]; otherwise f(x[i])=y[i-1]; #don't need this, actually
		meanPi0<-exp(.mylogSumExp(logpdensity,a=piVec[,1],c=0))
		upperCred<-uniroot(stepfun(piV,cumSum-(1 - (1 - CILevel)/2),right=TRUE),interval=c(0,1))$root
		lowerCred<-uniroot(stepfun(piV,cumSum-(1-CILevel)/2,right=TRUE),interval=c(0,1))$root
		ciMat<-rbind(c(lowerCred,upperCred),c(1-upperCred,1-lowerCred))
		colnames(ciMat)<-c("2.5%" ,    "97.5%")
		return(list(pi = c(meanPi0,1-meanPi0), piCI = ciMat, modeFit = c(NA,NA), 
	        nsim = nsim)) #return the function that creates the logposterior for checking purposes.
	    
    }
	else{
	    logPosterior<-function(theta){.logPosterior(theta,x=x,m=m,alleleSet=alleleSet,history=history,alpha=alpha)}
		fit <- try(.laplace(logPosterior, pi2theta(init)))
	    if (inherits(fit, "try-error")) {
	        fit <- try(.laplace(logPosterior, pi2theta(rep(1, length(init))/length(init))))
	    }
	    tpar <- list(m = fit$mode, var = 2 * fit$var, df = tdf)
	 	if (method == "sir") {
	        theta.s <- .sir(logPosterior, tpar, nsim)
	    }


	    if (method == "rej") {
	        Tdiff <- function(theta, datapar) {
	            data <- datapar$data
	            fit <- datapar$fit
	            d <- logPosterior(theta, x) - LearnBayes::dmt(theta, mean = fit$mode, 
	                S = 2 * fit$var, df = tdf, log = TRUE)
	            return(d)
	        }
	        maxD <- .laplace(Tdiff, fit$mode, list(data = x, fit = fit))
	        d <- Tdiff(maxD$mode, list(data = x, fit = fit))
	        theta.s <- .rejectsampling(logPosterior, tpar, d, nsim)
	    }
	    if (is.null(dim(theta.s))) {
	        pis <- sapply(theta.s, theta2pi)
	    }
	    else {
	        pis <- do.call(cbind, lapply(1:nrow(theta.s), function(i) {
	            theta2pi(theta.s[i, ])
	        }))
	    }
	    out <- t(apply(pis, 1, quantile, prob = c((1 - CILevel)/2, 
	        0.5, 1 - (1 - CILevel)/2)))
		meanOut<-	apply(pis, 1, mean)
	    return(list(pi = meanOut, piCI = out[, c(1, 3)], modeFit = theta2pi(fit$mode), 
	        nsim = nsim)) #return the function that creates the logposterior for checking purposes.
		
	}
}



