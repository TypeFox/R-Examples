surfaceSimulate <-
function(phy,type="BM",param=0,n_traits=NULL,dat=NULL,vcv=NULL,hansenfit=NULL, shifts=NULL,n_shifts=NULL,n_conv_shifts=NULL, n_regimes=NULL,n_per_regime=NULL,no_nested=TRUE, optima=NULL,  sample_optima=TRUE, optima_distrib=NULL, optima_type="rnorm", sigma_squared=NULL, alpha=NULL, pshift_timefactor=NULL){

if(type%in%c("BM","hansen-paint","hansen-fit")==F)stop("`type` must be `BM`, `hansen-paint` or `hansen-fit`")

ntaxa<-ifelse(class(phy)=="phylo",length(phy$tip.label),phy@nterm)

if(type=="BM"){
	
	if(is.null(n_traits)){
	if(!is.null(vcv)){ 
		n_traits<-dim(vcv)[1]
	}else if(!is.null(dat)){
		n_traits<-dim(dat)[2]
	}else{ 
		n_traits<-1
		}
	}
	if(is.null(vcv)){
		if(is.null(dat)){
			vcv<-diag(n_traits)
			dimnames(vcv)<-list(paste("V",1:n_traits,sep=""),paste("V",1:n_traits,sep=""))

		}else{
			vcv<-ratematrix(phy, dat)
		}	
	}
	
	phytransform<-function(phy,param){
		if(param==0)phy2<-phy
		if(param<0)phy2<-rescale(phy, model="EB", a = param)
		if(param>0)phy2<-rescale(phy, model="OU", alpha = param)
		phy2<-rescale(phy2,"depth", depth = max(branching.times(phy)))
		phy2
		}		

	if(is.null(param))param<-0
	if(length(param)==1|n_traits==1){
		phy2<-phytransform(phy,param=param[1])
		simdat<-data.frame(sim.char(phy2,vcv)[,,1])
	}else{
		if(length(param)!=n_traits)stop("`param` should either be a single value or a vector with number of elements equal to `ntrait`")
		simdat<-sim.char(phytransform(phy,param=param[1]),matrix(vcv[1,1]))[,,1]
		for(i in 2:n_traits){
			simdat<-cbind(simdat,sim.char(phytransform(phy,param=param[i]),matrix(vcv[i,i]))[,,1])
		}
	simdat<-data.frame(simdat)
	}
	names(simdat)<-colnames(vcv)
	optima<-shifts<-tip_regimes<-shifttimes<-hfit<-NULL

}else if(type=="hansen-fit"){

	if(class(hansenfit)=="hansentree"){
		n_traits<-1
	}else if(class(hansenfit[[1]])=="hansentree"){
		n_traits<-length(hansenfit)
	}else{
		stop("for type=`hansen-fit`, object `hansenfit` must either be a fitted hansentree object or a list of such objects (e.g. the `fit` component of a SURFACE output)")
		}

tempfit<-hansenfit
otree<-hansenfit[[1]] 
if(!is.null(alpha)){
	if(length(alpha)==1&n_traits>1)alpha<-rep(alpha,n_traits)
	for(i in 1:length(tempfit))tempfit[[i]]@sqrt.alpha<-sqrt(alpha[i])
	}
if(!is.null(sigma_squared)){
	if(length(sigma_squared)==1&n_traits>1)sigma_squared<-rep(sigma_squared,n_traits)
	for(i in 1:length(tempfit))tempfit[[i]]@sigma<-sqrt(sigma_squared[i])
	}
n_regimes<-length(unique(shifts))
if(sample_optima){
	if(is.null(optima_distrib)){
		optima_distrib<-rbind(sapply(tempfit,function(x)mean(x@theta$x)),sapply(tempfit,function(x)sd(x@theta$x)))
		}else{
		if(length(optima_distrib)==2)optima_distrib<-matrix(optima_distrib,nrow=2,ncol=n_traits)
	}
	optima<-matrix(NA,ncol=n_traits,nrow=n_regimes,dimnames=list(sort(unique(shifts)),NULL))
	for(i in 1:n_traits){
		if(optima_type=="rnorm")optima[,i]<-rnorm(n_regimes,mean=optima_distrib[1,i],sd=optima_distrib[2,i])
		if(optima_type=="runif")optima[,i]<-runif(n_regimes,min=optima_distrib[1]-optima_distrib[2]/2,max=optima_distrib[1]+optima_distrib[2]/2)
		if(optima_type=="even")optima[,i]<-sample(seq(optima_distrib[2],by=optima_distrib[2]/2,length.out=n_regimes))
		}
	}
if(!is.null(optima)){
	if(class(optima)=="numeric")optima<-matrix(optima,ncol=1)
	if(dim(optima)[1]!=n_regimes)stop("If setting optima must either sample or set `n_regime`")
	for(i in 1:length(tempfit))tempfit[[i]]@theta[[1]]<-optima[,i]
	}

}else if(type=="hansen-paint"){

Letters<-c(letters,paste("z",letters,sep=""),paste("zz",letters,sep=""),paste("zzz",letters,sep="")) 

if(!is.null(dat))n_traits<-dim(dat)[2]
if(is.null(n_traits))n_traits<-1

if(class(phy)%in%c("phylo","ouchtree")==FALSE)stop("`phy` must either be a `phylo` object or an `ouchtree` object")
if(class(phy)=="phylo"){
	if(is.null(dat))dat<-as.data.frame(matrix(rnorm(n_traits*ntaxa),ncol=n_traits,dimnames=list(phy$tip.label,NULL)))
	olist<-convertTreeData(phy,dat) 
	otree<-olist[[1]];odata<-olist[[2]]
	if(!is.null(shifts))warning("`phy` converted to `ouchtree` object: `shifts` provided must correspond to `@nodes` in this format",call.=FALSE)
	}else{
	otree<-phy
	odata<-as.data.frame(matrix(rnorm(n_traits*otree@nnodes),ncol=n_traits,dimnames=list(as(otree,"data.frame")$nodes,NULL)))
	odata[-otree@term,]<-NA
}

if(!is.null(shifts)){ #if shifts and branches specified, no sampling involved
	n_regimes<-length(unique(shifts))
}else{
if(!is.null(n_per_regime)){  #if provided, n_shifts & n_regimes known
	n_shifts<-sum(n_per_regime)   
	n_regimes<-length(n_per_regime)
	}
	if(is.null(n_shifts)|n_shifts==1){    	#if n_shifts=1 (default), placed ancestrally and no sampling involved
		n_shifts<-1
		shifts<-c("1"="a")
	}else{
		shifts<-character(n_shifts)
	if(is.null(n_conv_shifts)&is.null(n_regimes)&is.null(n_per_regime)){  #no convergence specified
		n_per_regime<-rep(1,n_shifts)
		shifts[]<-Letters[1:length(shifts)]
	}else{
		if(is.null(n_per_regime)){   #if don't need to sample shifts per regime, skip to building shifts
		if(!is.null(n_conv_shifts)&!is.null(n_regimes))stop("provide only one of `n_conv_shifts` or `n_regimes`")
		if(!is.null(n_regimes)){
			n_conv_regimes<-sample(min(1,n_shifts-n_regimes):min(floor(n_regimes/2),n_shifts-n_regimes),1)
			if(n_conv_regimes>0){
			n_per_conv_reg<-2+as.numeric(rmultinom(1,n_shifts-2*n_conv_regimes-(n_regimes-n_conv_regimes),prob=rep(1,n_conv_regimes)))
			}else{
				n_per_conv_reg<-NULL
				}
			n_per_regime<-c(rep(1,n_regimes-n_conv_regimes),n_per_conv_reg)
			n_per_regime<-n_per_regime[n_per_regime>0]
		}else{
			n_conv_regimes<-sample(1:floor(n_conv_shifts/2),1)
			n_per_regime<-c(rep(1,n_shifts-n_conv_shifts),2+as.numeric(rmultinom(1,n_conv_shifts-2*n_conv_regimes,prob=rep(1,n_conv_regimes))))
			n_regimes<-length(n_per_regime)
		}
	}
	shifts[]<-c("a",sample(rep(Letters[2:length(n_per_regime)],times=n_per_regime[-1])))
	}	#moved this close brace here - was after the shift sampling step...

	
if(is.null(pshift_timefactor))pshift_timefactor<-1
brtimes<-getBranchTimes(otree)
probshift<-1+brtimes*(pshift_timefactor-1)/max(as(otree,"data.frame")$times)

while(1){
	names(shifts)<-c("1",sample(otree@nodes[-1],n_shifts-1,prob=probshift[-1]))
	checknodes<-0
	for(i in 2:length(shifts)){
		checknodes[i]<-any(names(shifts[which(shifts==shifts[i])])%in%ouchDescendants(names(shifts)[i],otree))
	}
if(sum(checknodes)==0|no_nested==FALSE)break
		}	#one too many '}'s here - skipped sampling if is.null(n_regimes)&is.null(n_conv_shifts)???
}	}

regs<-repaint(otree,regshifts=shifts)
tempfit<-apply(odata,2,function(x)hansen(x,otree,regimes=regs,sqrt.alpha=0.2,sigma=1))



if(!is.null(alpha)){
	if(length(alpha)==1&n_traits>1)alpha<-rep(alpha,n_traits)
	for(i in 1:length(tempfit))tempfit[[i]]@sqrt.alpha<-sqrt(alpha[i])
	}
if(!is.null(sigma_squared)){
	if(length(sigma_squared)==1&n_traits>1)sigma_squared<-rep(sigma_squared,n_traits)
	for(i in 1:length(tempfit))tempfit[[i]]@sigma<-sqrt(sigma_squared[i])
	}
if(is.null(optima)){
	if(is.null(optima_distrib))optima_distrib<-c(0,1)
	if(optima_type=="rnorm")optima<-matrix(rnorm(n_traits*n_regimes,mean=optima_distrib[1],sd=optima_distrib[2]),ncol=n_traits,dimnames=list(sort(unique(shifts)),NULL))
	if(optima_type=="runif")optima<-matrix(runif(n_traits*n_regimes,min=optima_distrib[1]-optima_distrib[2]/2,max=optima_distrib[1]+optima_distrib[2]/2),ncol=n_traits,dimnames=list(sort(unique(shifts)),NULL))
	if(optima_type=="even"){
		optima<-matrix(NA,nrow=n_regimes,ncol=n_traits,dimnames=list(sort(unique(shifts)),NULL))
		for(i in 1:n_traits)optima[,i]<-sample(seq(optima_distrib[2],by=optima_distrib[2]/2,length.out=n_regimes))
#	for(i in 1:n_traits)optima[,i]<-sample(seq(optima_distrib[1]-optima_distrib[2]/2,optima_distrib[1]+optima_distrib[2]/2,length.out=n_regimes))
		}
	}else{
		if(class(optima)=="numeric")optima<-matrix(optima,ncol=1)
		if(any(dim(optima)!=c(n_regimes,n_traits)))stop("Optima must be provided as a matrix with dimensions [n_regimes, n_traits]")
		}

for(i in 1:length(tempfit))tempfit[[i]]@theta$x[]<-optima[,i]

}

if(type%in%c("hansen-paint","hansen-fit")){
	simdat<-as.data.frame(matrix(unlist(sapply(tempfit,simulate)),ncol=n_traits))
	#added this switch so ouchtree or phylo trees could both be used
#	if(class(phy)=="ouchtree"){
#		rownames(simdat)<-as(otree,"data.frame")$labels[
#	}else{
	if(class(phy)=="phylo"){
		simdat<-simdat[match(phy$tip.label,as(otree,"data.frame")$labels),,drop=F]
		rownames(simdat)<-phy$tip.label
	}
	names(simdat)<-paste("V",1:n_traits,sep="")
	shifttimes<-getBranchTimes(tempfit[[1]])[as.numeric(names(shifts))]
	tip_regimes<-as(tempfit[[1]],"data.frame")
	tips<-((dim(tip_regimes)[1]+1)/2):dim(tip_regimes)[1]
	tip_regimes<-data.frame(regs=tip_regimes[tips,5],row.names=tip_regimes$labels[tips])	
	hfit<-tempfit
	}

return(list(data=simdat,optima=optima,savedshifts=shifts,regimes=tip_regimes,shifttimes=shifttimes,fit=hfit))
}
