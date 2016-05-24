# B is the number of bootstrap observations you would like to generate
# eventOrdering is the list returned by eventTiming
#optional:
# pi is the estimate of pi you want to use for bootstrapping
# x, m the data to bootstrap
# note that assumes that the $alleleSet from '$call' is the right alleles (i.e. after seqError and norm. Cont. included)
bootstrapEventTiming <- function(eventOrdering, B,  type=c("parametric","nonparametric"),  pi, x, m,call){
	type<-match.arg(type)
	if(missing(x) | missing(m)){
		if(!"perLocationProb"%in% names(eventOrdering)) stop("perLocationProb must be one of the named elements of 'eventOrdering' if 'x' and 'm' are not supplied.")
		x<-eventOrdering$perLocationProb$x
		m<-eventOrdering$perLocationProb$m
	}
	if(missing(pi)){
		if(!"pi"%in%names(eventOrdering) | (any(is.na(eventOrdering$pi)) & type=="parametric")) stop("invalid values of 'pi' in 'eventOrdering' object")
		pi<-eventOrdering$pi
	}
	if(missing(call)){
		if(!"call"%in%names(eventOrdering)) stop("invalid values of 'call' in 'eventOrdering' object")
		call<-eventOrdering$call
	}
	N<-length(m)
	eventFunctionNames<-c("history","totalCopy" ,"method","seqError","type","normCont","coverageCutoff","minMutations","init","maxiter","tol")
	requiredNames<-c("alleleSet",eventFunctionNames)
	if(!all(requiredNames%in%names(call))) stop("Missing call information:",requiredNames[!requiredNames%in%names(call)])
	possAlleles<-call$alleleSet
	A<-call$history
	initial <- call$init;

	if(type=="parametric"){
		if(ncol(A)!=length(pi)) stop("'history' from 'call' does not match length of 'pi'")
		q<-as.vector(A%*%pi) 
		q<-q/sum(q)
		if(length(q)!=length(possAlleles)) stop("'alleleSet' and 'history' from 'call' does not match in size ")
		#deal with any possible problems with numerical inaccuracy
		if(any(q< 0)){
			if(any(q< -1e-10)) stop("Programming error -- negative probabilities")
			else q[q<0]<-0
		} 
		if(any(q> 1)){
			if(any(q> 1+1e-10)) stop("Programming error -- >1 probabilities")
			else q[q>1]<-1
		}
		bootData<-.parametricSimulation(B,m,q,alleleSet=possAlleles,onlyReturnData=TRUE)
	}
	else{
		wh<-matrix(sample(1:N,N*B,replace=TRUE),nrow=N,ncol=B)
		bootData<-lapply(1:B,function(kk){cbind(nMutAllele=x[wh[,kk]],nReads=m[wh[,kk]])})
	}
	piBoot<-lapply(bootData,function(z){do.call(eventTiming,c(list(x=z[,"nMutAllele"],m=z[,"nReads"]),call[eventFunctionNames]))$pi})
	piBoot<-do.call("rbind",piBoot)
	return(piBoot);
}


#assumes alleleSet is the right one (already accounts for sequencing error, normal contamination, etc.)
.parametricSimulation<-function(B,m,q,alleleSet,onlyReturnData=FALSE){
	if(length(q)!=length(alleleSet)) stop("'alleleSet' and 'q' must be of the same length")
	if(is.null(dim(m))){
		nMut<-length(m)
		m<-matrix(m,ncol=B,nrow=nMut,byrow=FALSE)
	}
	else{
		nMut<-nrow(m)
		if(B!=ncol(m)) stop("'B' must be same as number of columns of m")
	}
	if(nMut==0){
		data.frame(AF=rep(NA,nMut),nReads=rep(NA,nMut),nMutAllele=rep(NA,nMut),obsAF=rep(NA,nMut))		
	}
	nAlleles<-length(q)
	#adjust the allele frequencies to account for Sequencing error:
	
	###########Random generation of the allele frequencies: #
	#number of each allele frequency will observe
	nPerCategory<-rmultinom(n=B,size=nMut,prob=q) #a nallele x B matrix of which allele frequency get for each 	
	#repeat each of the alleles according to nPerCategory and then permute them so that each m[i] equally likely to be paired with each allele
	AF<-apply(nPerCategory,2,function(x){sample(rep(alleleSet,times=x))}) #matrix of nMut x B

	#if B=1 of nMut=1, make it into matrix
	if(is.null(dim(AF))){
		AF<-matrix(AF,ncol=B,nrow=nMut,)
	}
	#run binomial for each B and turn into reasonable data.frame
	return(lapply(1:B,function(kk){
		x<-rbinom(nMut,size=m[,kk],prob=AF[,kk])
		if(onlyReturnData) cbind(nMutAllele=x,nReads=m[,kk])
		else data.frame(AF=AF[,kk],nMutAllele=x,nReads=m[,kk],obsAF=x/m[,kk])
	}))
}

readSimulation<-function(B,alleleSet, q, totalCopy,mutRate=NULL,seqError=0,fixedN=FALSE,normCont=0, aveReadCoverage=30,countDistribution=NULL){
	#from mutation rate, how many actual mutations
	if(fixedN & is.null(mutRate)) stop("must give mutRate if fixedN value")
	if(length(mutRate)!=1) stop("mutRate must be a single value")
	if(!fixedN) nMut<-rpois(B,mutRate)  #vector of length B
	else nMut<-mutRate #single value
	contAFSet<-contAF(alleleSet,totalCopy=totalCopy,normCont=normCont,type="mutation")

	contErrorAFSet<-errorAF(contAFSet,seqError=seqError)
	
	singleCall<-function(nMut,B){
		if(nMut>0){
			if(is.null(countDistribution)){
				m<-matrix(rpois(nMut*B,aveReadCoverage),nrow=nMut,ncol=B) #actual read coverage per mutation
			}
			else{ #otherwise, use empirical distribution given by countDistribution (i.e. ignore aveReadCoverage)
				aveReadCoverage<-round(mean(countDistribution),2) #so label is still relevant
				m<-matrix(sample(countDistribution,size=nMut*B,replace=TRUE),nrow=nMut,ncol=B)
			}
			simData<-.parametricSimulation(B=B,m=m,q=q,alleleSet=contErrorAFSet)
			simData<-lapply(simData,function(x){data.frame(tumorAF=alleleSet[match(x$AF,contErrorAFSet)],contAF=contAFSet[match(x$AF,contErrorAFSet)],x,normCont=unname(normCont))})
		}
		else simData<-list(data.frame(tumorAF=NA,contAF=NA,AF=NA,nMutAllele=NA,nReads=NA,obsAF=NA,normCont=NA)) #doesn't matter what normCont is, because 0 rows!
		return(simData)
	}
	if(fixedN) simData<-singleCall(nMut,B)
	else simData<-unlist(lapply(nMut,singleCall,B=1),recursive=FALSE) #have to call them separately...
	
	return(simData)
}

