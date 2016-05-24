#' Probability of being a sampled ancestor of another sampled taxon
#'
#' Uses models from Foote (1996) to calculate the probability 
#'
#' probAnc obtains the probability of sampling a descendant of a 
#' morphotaxon in the fossil record, given the sampling probability 
#' and estimates of origination and extinction rates. These values are 
#' always calculated assuming infinite time for the potential ancestor
#' to produce daughter taxa (assuming it lives that long) and under 
#' homogenous birth, death and sampling rates/probabilities, which is a
#' situation that may be overly ideal relative to many real fossil records.
#' 
#' This can be calculated for either direct descendants, i.e. the probability 
#' of sampling any morphotaxa that arise immediately from the particular 
#' morphotaxon that could be an ancestor, or indirect descendants, i.e. the 
#' probability for any morphotaxon that has the morphotaxon of question as an
#' ancestor, no matter how distant. See the argument \code{analysis} for 
#' details. Mode of differentiation can also be varied
#' for three different models, see the argument \code{mode}.
#'
#' This probability is calculated including the probability that extinction might
#' occur before any descendants are produced. Thus, if \code{p = q}, the probability of 
#' a taxon going extinct before it produces any descendants will be 0.5, which 
#' means that even when sampling is perfect (\code{R = 1}, meaning completeness of 
#' 100%) the probability of a taxon being an ancestor of another sampled taxon
#' can be no higher than 0.5. See Foote (1996) for a graphic depiction of this
#' non-intuitive ceiling. For reasons (probably?) having to do with finite
#' approximations of infinite summations, values close to perfect sampling
#' may have values slightly higher than this ceiling, which is also apparent
#' visually in the figures in Foote (1996). Thus, values higher than 0.5 when p=q
#' should be discounted, and in general when sampling rate is high, results should
#' be treated cautiously as overestimates.

#' @inheritParams SamplingConv

#' @param analysis The type of analysis to be performed, either the probability of sampling direct
#' descendants (\code{"directDesc"}) or of sampling indirect descendants (\code{"indirectDesc"}).

#' @param Mmax The maximum number of direct descendants (M) to sum over in the function, which
#' is ideally meant to be a sum from zero to infinity, like nrep. Unfortunately,
#' \code{(2*M)} is used in a factorial, which means we are limited to a relatively
#' small upper bound on M.

#' @seealso \code{\link{SamplingConv}} 

#' @references 
#' Foote, M. 1996 On the Probability of Ancestors in the Fossil
#' Record. \emph{Paleobiology} \bold{22}(2):141--151.

#' @examples
#' #examples, run at very low nrep for sake of speed (examples need to be fast)
#'
#' #default: probability of sampling a direct descendant
#' probAnc(p = 0.1, q = 0.1, R = 0.5, mode = "budding", analysis="directDesc",nrep=100)
#' 
#' #other modes
#' probAnc(p = 0.1, q = 0.1, R = 0.5, mode = "bifurcating", analysis="directDesc",nrep=100)
#' probAnc(p = 0.1, q = 0.1, R = 0.5, mode = "anagenesis", analysis="directDesc",nrep=100)
#'
#' #probability of having sampled indirect descendants of a taxon
#' probAnc(p = 0.1, q = 0.1, R = 0.5, mode = "budding", analysis="indirectDesc",nrep=100)	#default
#' probAnc(p = 0.1, q = 0.1, R = 0.5, mode = "bifurcating", analysis="indirectDesc",nrep=100)
#' probAnc(p = 0.1, q = 0.1, R = 0.5, mode = "anagenesis", analysis="indirectDesc",nrep=100)
#'

#'@export
probAnc<-function(p,q,R,mode="budding",analysis="directDesc",Mmax=85,nrep=10000){
	#see Foote, 1996	
	#calculates prob of taxa with indirect desc under budding speciation
		#under infinite time, with p=q or p<q
	#When mode=anagenesis, p is taken to be the rate of anagenesis
	#unused equations related to estimating indirect ancestry? DWB 12-09-13
		#Pinf<-function(p,q,Pp){
		#	res<-numeric()
		#	for(M in 1:80){
		#	res[M]<-Qm(p,q,M)*(1-(1-Pp)^M)
		#		}
		#	sum(res)
		#	}
	#test
	if(!any(mode==c("budding","bifurcating","anagenesis"))){
		stop("Mode not designated, must be 'budding', 'bifurcating' or 'anagenesis'")}
	if(mode=="anagenesis"){message("p will be treated as the rate of anagenesis/pseudospeciation")}
	if(!any(analysis==c("directDesc","indirectDesc"))){
		stop("Analysis type not designated, must be 'directDesc' or 'indirectDesc'")}
	if(nrep<0){stop("Nrep is less than zero?")}
	if(analysis=="directDesc"){
		#get completeness
		Pp<-qsProb2Comp(R=R,p=p,q=q,mode=mode)
		#functions dependent on mode
		if(mode=="budding"){
			Pd<-function(p,q,Ti){
				exp(-q*(Ti-1))-exp(-q*Ti)
				}
			PN<-function(p,q,Ti,Ni){
				(exp(-p*Ti)*((p*Ti)^Ni))/factorial(Ni)
				}
			#should approximate 2*(p/q)*Pp when R and p are both <<1
			#approx<-2*(p/q)*Pp
			maxN<-100
			}
		if(mode=="bifurcating"){
			Pd<-function(p,q,Ti){exp(-(p+q)*(Ti-1))-exp(-(p+q)*Ti)}
			PN<-function(p,q,Ti,Ni){
				if(Ni==0){res<-q/(p+q)}
				if(Ni==2){res<-p/(p+q)}
				if(Ni!=2 & Ni!=0){res<-0}
				return(res)
				}
			#approx<-2*(p/(p+q))*Pp
			maxN<-2
			}
		if(mode=="anagenesis"){
			Pd<-function(p,q,Ti){exp(-(p+q)*(Ti-1))-exp(-(p+q)*Ti)}
			PN<-function(p,q,Ti,Ni){
				if(Ni==0){res<-q/(p+q)}
				if(Ni==1){res<-p/(p+q)}
				if(Ni!=1 & Ni!=0){res<-0}
				return(res)
				}
			#approx<-2*(p/(p+q))*Pp
			maxN<-1
			}
		#now get PA, the probability of a direct descendant of a taxon being sampled
		res<-numeric()
		for(t in 1:nrep){
			Nres<-numeric()
			for(N in 0:maxN){
				Nres[N]<-PN(p=p,q=q,Ti=t,Ni=N)*(1-((1-Pp)^N))
				}
			res[t]<-(1-((1-R)^t))*Pd(p=p,q=q,Ti=t)/Pp*sum(Nres)
			}
		}
	#now prob of sampling an indirect descendant over indefinite time 
		#(nrep in the case of anagenesis)
	if(analysis=="indirectDesc"){
		if(mode=="budding" | mode=="bifurcating"){
			if(p>q){stop(
				"Indirect Descendant formulae are unsolved if p>q, see Foote 1996")}
			Qm<-function(p,q,M){
				x<-(4*p*q)/((p+q)^2)
				res<-((p+q)/(2*p))*(factorial(2*M)/((2^(2*M))*factorial(M)^2))*((x^M)/((2*M)-1))
				return(res)
				}
			Pp<-qsProb2Comp(R=R,q=q,mode="budding")	#as per instructions in Foote, 1996, use 'budding' for both
			#get prob using P'		
			res<-numeric()
			for(M in 1:Mmax){	#85 is the maximum number we can do factorial(2*M) too... but doesn't matter mostly
				res[M]<-Qm(p=p,q=q,M=M)*(1-(1-Pp)^M)
				}
			}
		if(mode=="anagenesis"){
			Pp<-qsProb2Comp(R=R,q=q,mode="anagenesis")
			QmStar<-function(p,q,M,T){
				firstTerm<-numeric()
				for(t in 1:(T-1)){
					firstTerm<-(exp(-q*(t-1))-exp(-q*t))*(exp(-p*t)*((p*t)^(M-1)))/factorial(M-1)
					}
				secondTerm<-exp(-q*(T-1))*(exp(-p*T)*((p*T)^(M-1)))/factorial(M-1)
				res<-sum(firstTerm)+secondTerm
				return(res)
				}
			res<-numeric()
			for(T in 1:nrep){
				Tres<-numeric()
				for(M in 1:Mmax){	
					Tres[M]<-QmStar(p=p,q=q,M=M,T=T)*(1-(1-Pp)^M)
					}
				res[T]<-sum(Tres)
				}
			}
		}
	if(any(is.nan(res))){
		message("Input parameters and nrep produce NaN values, which are replaced with zeroes.")
		message("May want to decrease nrep to see if returned estimate holds.")
		res[is.nan(res)]<-0
		}
	res<-sum(res)
	names(res)<-NULL
	if(res>0.5 & p==q){
		message("Treat result with caution: if p = q, then prob of a taxon being an ancestor should be no greater than 0.5.")
		message("Values higher than 0.5 result from limits of finite calculates, particularly with high sampling probabilities.")
		message("See documentation.")
		}
	return(res)
	}
	
