# Functions for compute spatial variance-covariance matrix and other spatial SEM steps
 
# Updated Fan 22, 2014

# Citation: Lamb, E.G., K. Mengersen, K.J. Stewart, U. Attanayake, S.D. Siciliano. 2014. Spatially explicit 
# structural equation modeling. Ecology xx:xxx-xxx.

# Functions here were tested using R 3.0.2
# libraries "lavaan", "gplots", and "mgcv" are required


# calc.dist
	#
	# calculation of a distance matrix
	# datafile is a dataset with the first two columns X,Y coordinates
#' @export
calc.dist<-function(datafile){
	aa.tmp <- na.omit(datafile)					# omit rows with NA's (can impute later)
	aa.x = length(aa.tmp[,1])					# no. rows
	aa.dist=0								# compute distance lags
		for (k in 1:aa.x){
		for (l in k:aa.x){
		aa.dist=c(aa.dist,sqrt((aa.tmp[k,1]-aa.tmp[l,1])^2+(aa.tmp[k,2]-aa.tmp[l,2])^2))
		}}
	aa.dist=aa.dist[2:length(aa.dist)]
}


#make.bin(dist.mat,type="n.bins",p.dist=50,n.bins=10,s.size=100)

	# generates cut off values for lag distance bins and corresponding bin names
	# The function has three default parameter values available, if user does not want to specify:
		#Inference distance as a percentage(p.dist) = 50%
		#Number of bins (n.bins) = 10
		#Sample size (s.size) = 100
	#dist.mat is a distance matrix produced by calc.dist
	#Can use type="ALL","n.bins" OR "s.size" to control parameter values.
	#The function produces a list object containing (1.)binsize and (2.)binname
	#These two vectors (binsize and binname) will be used by make.covar to calculate variance covariance matrices for each lag distance bin
	
	#Special note:
	#User specified number of lag distance bins OR sample size 
	#will be used to calculate initial cutoff value of each lag distance bin.
	#However, if the cutoff point is in between a lag distance bin, 
	#real cutoff will apply at the upper margin of the particular bin. 
	#Therefore, resulting number of bins are less than or equal AND 
	#resulting sample sizes are greater than or equal to the value specified by the user.

#' @export
make.bin<-function(dist.mat,type="n.bins",p.dist=50,n.bins=10,s.size=100){
	if((type=="ALL")||(type=="")){
		cat("Inference distance = ",p.dist,"%",sep="",fill = TRUE)
		cat("Notification:inference distance = specifi. by user AND number of bins = ALL",fill = TRUE)	
		infr.dist.frac<-p.dist/100
		dist.max<-max(dist.mat)
		dist.infr<-round((dist.max*infr.dist.frac),digits=2)
		dist.temp<-sort(unique(dist.mat))
		dist.logic<-ifelse(dist.temp>dist.infr,FALSE,TRUE)
		binsize<-dist.temp[dist.logic==TRUE]
		cat("Bin size: ",binsize,sep=" ",fill = TRUE)
		colnames<-paste("Bin",1:(length(binsize)),sep="")
		binname.m<-matrix(colnames)
		binname.c<-as.character(binname.m[-length(binname.m),])
		binname<-c(binname.c)
		cat("Bin names: ",binname,sep=" ",fill = TRUE)
	}
	if(type=="n.bins"){
		s.size<-0
		cat("Inference distance = ",p.dist,"%",sep="",fill = TRUE)
		cat("Specified number of bins = ",n.bins,sep="",fill = TRUE)
		cat("Notification:inference distance = specifi. by user AND number of bins = specifi. by user",fill = TRUE)
	}
	if(type=="s.size"){
		n.bins<-0
		cat("Inference distance = ",p.dist,"%",sep="",fill = TRUE)
		cat("Specified sample size = ",s.size,sep="",fill = TRUE)
		cat("Notification:inference distance = specifi. by user AND sample size = specifi. by user",fill = TRUE)
	}
	if((type=="n.bins")||(type=="s.size")){
	infr.dist.frac<-p.dist/100
	dist.max<-max(dist.mat)
	dist.infr<-round((dist.max*infr.dist.frac),digits=2)
	dist.temp<-sort(unique(dist.mat))
	dist.logic<-ifelse(dist.temp>dist.infr,FALSE,TRUE)
	binsize.freq.tot<-as.data.frame(table(rbind(dist.mat)))
	binsize.freq<-binsize.freq.tot[dist.logic==TRUE,]
	number.zeros<-binsize.freq.tot[1,2]
	sample.tot<-sum(binsize.freq$Freq)-number.zeros
	sample.points<-matrix(as.data.frame(lapply(binsize.freq,function(Var1)rep(Var1,binsize.freq$Freq)))[,-2])
	if(n.bins!=0){
		s.size<-round((sample.tot/n.bins),digits=0)
	}else
		n.bins<-round((sample.tot/s.size),digits=0)
	sample.vector.null<-matrix(0,1,n.bins)
	null.length<-length(sample.vector.null)
	k<-number.zeros
	for(i in 1:null.length){
		temp.sum<-k+s.size
		if(temp.sum<=sum(binsize.freq$Freq)){
			val.temp.sum.1<-sample.points[temp.sum,]
			if(temp.sum!=sum(binsize.freq$Freq)){
				val.temp.sum.2<-sample.points[(temp.sum+1),]
			}else
				val.temp.sum.2<-0
			if(val.temp.sum.1==val.temp.sum.2){
				repeat
					{
					val.temp.sum.1<-sample.points[temp.sum,]
						if(temp.sum!=sum(binsize.freq$Freq)){
							val.temp.sum.2<-sample.points[(temp.sum+1),]
						}else
							val.temp.sum.2<-0
					if(val.temp.sum.1==val.temp.sum.2){
						temp.sum<-temp.sum+1
					}else
						temp.sum<-temp.sum
					if(val.temp.sum.1!=val.temp.sum.2)break
					sample.vector.null[1,i]<-temp.sum
					k<-temp.sum
					}
			}else
				sample.vector.null[1,i]<-temp.sum
				k<-temp.sum
		}else
			sample.vector.null[1,i]<-0
			k<-temp.sum
	}
	temp.logic<-ifelse(sample.vector.null>0,TRUE,FALSE)
	sample.omit.zeros.temp<-sample.vector.null[1,1:(length(temp.logic[temp.logic==TRUE]))]
	sample.omit.zeros<-c(1,sample.omit.zeros.temp)
	binsize<-sort(dist.mat)[sample.omit.zeros]
	cat("Bin size: ",binsize,sep=" ",fill = TRUE)
	colnames<-paste("Bin",1:(length(binsize)),sep="")
	binname.m<-matrix(colnames)
	binname.c<-as.character(binname.m[-length(binname.m),])
	binname<-c(binname.c)
	cat("Bin names: ",binname,sep=" ",fill = TRUE)
	}
return(list(	
bin_size<-binsize,  #vector of bin size cut off distances
bin_names<-binname  #vector of bin names for each lag distance bin
))
}



# plotbin - used to visualize sample sizes within bins
	# 
	# dist.mat is a distance matrix produced by calc.dist
	#
	# binsize is a vector of lag distances starting at 0. The following values are the upper limits of each distance bin
		# binsize should have n+1 elements where n is the number of lag distance bins desired
		
#' @export
	
plotbin<-function(dist.mat,binsize) {
	par(mfrow=c(1,2),mar=c(7,7,2,2),mgp=c(5,1,0))
	hist(dist.mat[(dist.mat>0)==T],col="gray",las=2,las=2,xlab="distance",ylab="# sample pairs",cex.lab=1.5,main="Histogram of all Pairs")
	barplot(tapply(dist.mat[dist.mat<=max(binsize)],cut(dist.mat[dist.mat<=max(binsize)],breaks=binsize),length),las=2,xlab="lag distance bins",ylab="# sample pairs",cex.lab=1.5,main="Frequency in Selected Bins")
	abline(h=0)
	par(mfrow=c(1,1),mar=c(5,4,2,2),mgp=c(3,1,0))
}

#make.covar - Calculate covariance matrices
	
	# calculates variance covariance matrices for each lag distance bin and for a flat (non-spatial) bin
	# produces list object with [[1]]bin.summary, [[2]] variable names [[3]] flat covariance matrix, [[3]][,,i] covariance matrices for each bin i
	# summary of bins printed
	#
	# datafile is a dataset with the first two columns in data file X,Y coordinates
	# dist.mat is a distance matrix produced by calc.dist
	# binsize is a vector of lag distances starting at 0. The following values are the upper limits of each distance bin
		# binsize should have n+1 elements where n is the number of lag distance bins desired
	# binname is a vector of same length as aa.bin with names for each lag distance bin

#' @export
make.covar<-function(datafile,dist.mat,binsize,binname) {
	aa.tmp <- na.omit(datafile)					# omit rows with NA's (can impute later)
	aa.x = length(aa.tmp[,1])					# no. rows
	aa.y = length(aa.tmp[1,])					# no. covariates
	aa.y2 = aa.y-2							# no. covariates, excluding X,Y
	varnames<-as.factor(names(aa.tmp[,c(3:aa.y)]))
	names(aa.tmp)<-seq(1:aa.y)					# shorten names for convenience
	aa.z = length(binsize)
	aa.z1 = aa.z-1
	aa.t <- matrix(0,aa.x,aa.y2)
	for (i in 3:aa.y){			
		aa.t[,i-2] = (aa.tmp[,i]-mean(aa.tmp[,i]))}	# standardise (x-xbar) 
	aa.tgt = array(0,dim=c(aa.y2,aa.y2,aa.z1)) 		# set up for covar matrices 
	aa.lagcnt=rep(0,aa.z1)						# compute covariances for each bin size; no. entries in each lag group
	aa.c = 0								# internal counter
	for (k in 1:aa.x){    						# loop over observations
		for (l in k:aa.x){
			aa.c = aa.c + 1
 			aa.cut<-cut(dist.mat[aa.c],breaks=binsize, labels=1:aa.z1)
			 aa.lagcnt[aa.cut] <- aa.lagcnt[aa.cut]+1
	for (i in 1:aa.y2){  						# for each cell in cov matrix
		for (j in 1:aa.y2){	
		aa.tgt[i,j,aa.cut] <- aa.tgt[i,j,aa.cut] + (aa.t[k,i]-aa.t[l,i])*(aa.t[k,j]-aa.t[l,j])	# compute covar	
}}
}}  										# for k,l
	aa.cov <- aa.tgt[,,1:aa.z1]
		for (m in 1:aa.z1){					# divide by 2N 
		aa.cov[,,m] <- aa.tgt[,,m]/(2*aa.lagcnt[m])
}
	aa.cov.0 <- cov(aa.t[,1:aa.y2])				# standard (no-lag) covariance matrix
	
	binlow<-binsize[1:(length(binsize)-1)]			#produces summary of the bins
	binhigh<-binsize[2:(length(binsize))]
	samplesize<-aa.lagcnt
	bin.summary<-data.frame(binname,binlow,binhigh,samplesize)
	flatrow<-data.frame("binflat",0,0,aa.x)
	colnames(flatrow)<-colnames(bin.summary)
	print(bin.summary<-rbind(bin.summary,flatrow))
	return(list(bin.summary,varnames,aa.cov.0,aa.cov))
}


# runModels(spatial_model,covdata)
#		given an sem model (spatial_model) specified using lavaan syntax,
#			and a list object containing covariance matrices and other descriptors,
#			as produced by make.covar, fits the sem model to the non-spatial covariance matrix, 	
#			and for the covariance matrices corresponding to all lag distance bins found there.
#		Produces a list object containing:
#			[1] a table of model fit estimates for each model. 
#					Values available in this table include the following: 
#					"fmin","chisq","df","pvalue","baseline.chisq","baseline.df","baseline.pvalue",
#					"cfi","tli","nnfi", "rfi","nfi","pnfi","ifi","rni","logl","unrestricted.logl",
#					"npar","aic","bic","ntotal","bic2","rmsea","rmsea.ci.lower","rmsea.ci.upper",
#					"rmsea.pvalue","rmr","rmr_nomean","srmr","srmr_nomean","cn_05","cn_01","gfi","agfi",
#					"pgfi","mfi","ecvi". See the lavaan documentation for an explanation of each value.
#			[2] table containing a vector of parameter numbers and a character vector containing the names 
#					of the paths included in each model.
#			[3] unstandardized path coefficient estimates
#			[4] standard error of unstandardized path coefficient estimates
#			[5] p-values for each unstandardized path coefficient estimate
#			[6] standardized parameter estimates for each path in each model
#			[7] character vector containing list of names of dependent variables within the models
#			[8] r-square values for each dependent variable in each model
#			[9] names of each path for which there is a modification index value
#			[10] modification index values for each potential path addition for each model
#			[11] a copy of the bin.summary table in the input covdata object

#' @export
#' @import lavaan
runModels<-function(spatial_model,covdata){
	bin.summary<-covdata[[1]]
	model.fit<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(model.fit)<-bin.summary[,1]
	par.name<-matrix(numeric(1*length(bin.summary[,1])),ncol=1)
	unst.est<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(unst.est)<-bin.summary[,1]
	est.se<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(est.se)<-bin.summary[,1]
	p.val<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(p.val)<-bin.summary[,1]
	std.est<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(std.est)<-bin.summary[,1]
	varname.r.square<-matrix(numeric(1*length(bin.summary[,1])),ncol=1)
	r.square<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(r.square)<-bin.summary[,1]
	mod.parname<-matrix(numeric(1*length(bin.summary[,1])),ncol=1)
	mod.indices<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(mod.indices)<-bin.summary[,1]
	for (i in 2:length(bin.summary[,1])-1) {
		bin_name<-bin.summary[i,1]
		covmatrix<-data.frame(covdata[[4]][,,i])
		names(covmatrix)<-covdata[[2]]
		covmatrix<-as.matrix(round(covmatrix,8))#rounding needed to make matrix symmetrical 
												#difference in number ofdigits saved between 
												#upper and lower; as.matrix needed to convert 
												#dataframe to matrix for sem
		semmodel<-sem(spatial_model,sample.cov=covmatrix,sample.nobs=bin.summary[i,4])
		if(inspect(semmodel,"converged")==TRUE) {
			model.fit[i]<-data.frame(fitMeasures(semmodel))
			unst.est[i]<-data.frame(parameterEstimates(semmodel)[,4])
			est.se[i]<-data.frame(parameterEstimates(semmodel)[,5])
			p.val[i]<-data.frame(parameterEstimates(semmodel)[,7])
			std.est[i]<-data.frame(parameterEstimates(semmodel,standardized=T)[,11])
			if(i==1) {varname.r.square<-(rownames(data.frame((inspect(semmodel,"rsquare")))))}
			r.square[i]<-data.frame(inspect(semmodel,"rsquare"))
			if(i==1) {mod.parname<-data.frame(paste(modificationIndices(semmodel)[,1],modificationIndices(semmodel)[,2],modificationIndices(semmodel)[,3]))
			colnames(mod.parname)<-c("parameter name")}
			mod.indices[i]<-data.frame(modificationIndices(semmodel)[,4])
}}
		covmatrix<-data.frame(covdata[[3]])#nonspatial (flat) matrix
		names(covmatrix)<-covdata[[2]]
		covmatrix<-as.matrix(round(covmatrix,8)) 
		semmodel<-sem(spatial_model,sample.cov=covmatrix,sample.nobs=bin.summary[length(bin.summary[,1]),4])
		par.name<-data.frame(seq(1:length(parameterEstimates(semmodel)[,1])),paste(parameterEstimates(semmodel)[,1],parameterEstimates(semmodel)[,2],parameterEstimates(semmodel)[,3]))
				colnames(par.name)<-c("parameter.number","parameter.name")
		model.fit[length(bin.summary[,1])]<-data.frame(fitMeasures(semmodel))
		unst.est[length(bin.summary[,1])]<-data.frame(parameterEstimates(semmodel)[,4])
		est.se[length(bin.summary[,1])]<-data.frame(parameterEstimates(semmodel)[,5])
		p.val[length(bin.summary[,1])]<-data.frame(parameterEstimates(semmodel)[,7])
		std.est[length(bin.summary[,1])]<-data.frame(parameterEstimates(semmodel,standardized=T)[,11])
		r.square[length(bin.summary[,1])]<-data.frame(inspect(semmodel,"rsquare"))
		mod.indices[length(bin.summary[,1])]<-data.frame(modificationIndices(semmodel)[,4])
	return(list(
	model_fit<-data.frame(model.fit,
	row.names=names(fitMeasures(semmodel))), # dataframe containing model fit indices for all bins
	par_name<-data.frame(par.name), # parameter names in rows
	unst_est<-data.frame(unst.est), # unstandardized parameter estimates in rows, column for each lag distance bin
	est_se<-data.frame(est.se), # std error of unstandardized parameter estimates in rows, column for each lag distance bin
	p_val<-data.frame(p.val), # p-value for each unstandardized parameter estimates in rows, column for each lag distance bin
	std_est<-data.frame(std.est),# standardized parameter estimates in rows, column for each lag distance bin
	varname_r_square<-data.frame(varname.r.square),#names of dependent variables
	r_square<-data.frame(r.square),#returns r2 values for dependent observed variables
	mod_parname<-data.frame(mod.parname),# list of all parameters for which there is a modification index in rows, column for each lag distance bin
	mod_indices<-data.frame(mod.indices),	# modification indices in rows, column for each lag distance bin
	bin.summary
	))
}	

# modelsummary(spatial_model_results)
#	spatial_model_results = object output from run.Models 
#
#	extracts basic model summary information from the bin.summary
#	file and the object created by run.Models in a readable format
#' @export
modelsummary<-function(spatial_model_results){
	bin.summary<-data.frame(spatial_model_results[[11]])
	model.summary <- data.frame(bin.summary$binname,round(t(data.frame(spatial_model_results[1])),digits=3))
	model.summary<-merge(bin.summary,model.summary,by.x="binname",by.y="bin.summary.binname")
	print(format(model.summary[order(model.summary[,3]),c("binname","binlow","binhigh","chisq","df","pvalue","cfi",
	"npar","aic","bic","rmsea","rmsea.ci.lower","rmsea.ci.upper","rmsea.pvalue")],digits=3))
}

# bin.results(spatial_model_results,bin="binflat") 
#	spatial_model_results = object output by run.Models 
#	bin= name of the bin that results are desired for. Defaults to flat (nonspatial) model

#	extracts path coefficients, standard errors, and standardized coefficients for a particular bin
#	in a readable format
#' @export
bin.results<-function(spatial_model_results,bin="binflat") {
	bin.summary <-data.frame(spatial_model_results[[11]])
	bincol<-match(bin,bin.summary[,1],nomatch=0)
	if(bincol==0) {print("Error: Bin name does not match")}
	if(bincol>0) {
		results<-data.frame(data.frame(spatial_model_results[2])[,2],round(data.frame(spatial_model_results[3])[,bincol],digits=3),
			round(data.frame(spatial_model_results[4])[,bincol],digits=3),round(data.frame(spatial_model_results[5])[,bincol],digits=3),
			round(data.frame(spatial_model_results[6])[,bincol],digits=3))
		names(results)<-c("pathname","unstd_coeff","se_unstd_coeff","pvalue","std_coeff")
		return(list(paste("Summary results for",bin),results))}
}

# bin.rsquare(spatial_model_results,bin="binflat") 
#	spatial_model_results = object output by run.Models
#	bin= name of the bin that results are desired for. Defaults to nonspatial (flat) model
#	extracts rsquare values for dependent variables for a particular binn in a readable format
#' @export
bin.rsquare<-function(spatial_model_results,bin="binflat") {
	bin.summary <-data.frame(spatial_model_results[[11]])
	bincol<-match(bin,bin.summary[,1],nomatch=0)
	if(bincol==0) {print("Error: Bin name does not match")}
	if(bincol>0) {
	results<-data.frame(data.frame(spatial_model_results[7]),round(data.frame(spatial_model_results[8])[,bincol],digits=3))
	names(results)<-c("parname","rsquare")
	return(list(paste("r-square values for dependent variables in",bin),results))}
	}	


# Step 6 plot results
#	plotmodelfit(spatial_model_results,plots="all",add.line="none",rmsea_err=T,pch=16,lwd=2,lty=1,cex=1,cex.lab=1)
#   working version of model fit plotting with optional lines, and formatting
#	spatial_model_results = object output by summary.models.lavaan
#	plots - options for selecting fit indices to plot
#		"all" indicates all of chi square, cfi, rmsea, and srmr to be plotted 
#		"chi", "cfi", "rmsea", "srmr" select individual plots
#	add.lines. options for plotting a fit line 
#		"none" indicates no line
#		"step" plots straight line segments between points
#		"smooth" plots a smoothed curve fit using function lowess. Smoothed lines do not include the flat model
#	rmsea_err. should confidence limits for rmsea be plotted
#	pch, lwd, lty cex, cex.lab, cex.axis, cex.main options for formatting points and fit lines.
#' @export
#' @import gplots
plotmodelfit<-function(spatial_model_results,plots="all",add.line="none",rmsea_err=T,pch=16,lwd=2,lty=1,cex=1,cex.lab=1,cex.axis=1,cex.main=1.5) {
	plot.chi<-ifelse(plots=="chi",T,ifelse(plots=="all",T,F))
	plot.cfi<-ifelse(plots=="cfi",T,ifelse(plots=="all",T,F))
	plot.rmsea<-ifelse(plots=="rmsea",T,ifelse(plots=="all",T,F))
	plot.srmr<-ifelse(plots=="srmr",T,ifelse(plots=="all",T,F))
	plot.all<-ifelse(plots=="all",T,F)
	if(plot.all==T) {par(mfrow=c(2,2))}
	if(plot.all==F) {par(mfrow=c(1,1))}
	bin.summary<-data.frame(spatial_model_results[[11]])
	bin.summary_noflat<-bin.summary[-length(bin.summary[,1]),]
	fitindices<-spatial_model_results[[1]]
	fitindices_noflat<-fitindices[,-length(fitindices)]
	bin.sequence<-order(bin.summary[,3])
	bin.sequence_noflat<-bin.sequence[c(2:length(bin.sequence))]
	add.line.step<-ifelse(add.line=="step",T,F)
	add.line.smooth<-ifelse(add.line=="smooth",T,F)
	if(plot.chi==T) {
	plot(matrix(bin.summary[,3]),as.numeric(matrix(fitindices[2,])),main="Chi Sq. vs Lag dist",pch=pch,xlab="Lag Distances",ylab="Model Chi Square Values",cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,cex.main=cex.main)
		abline(h=qchisq(0.95,fitindices[3,1]),lwd=2,lty=3)
		if(add.line.step==T) {lines(bin.summary[,3][bin.sequence],as.numeric(matrix(fitindices[2,]))[bin.sequence],lwd=lwd,lty=lty)}
		if(add.line.smooth==T) {lines(lowess(bin.summary_noflat[,3][bin.sequence_noflat],as.numeric(matrix(fitindices_noflat[2,]))[bin.sequence_noflat]),lwd=lwd,lty=lty)}
	}
	if(plot.cfi==T) {
	plot(matrix(bin.summary[,3]),as.numeric(matrix(fitindices[8,])),main="CFI vs Lag dist",pch=pch,xlab="Lag Distances",ylab="Model CFI Values",cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,cex.main=cex.main)
		abline(h=0.9,lwd=2,lty=3)
		if(add.line.step==T) {lines(bin.summary[,3][bin.sequence],as.numeric(matrix(fitindices[8,]))[bin.sequence],lwd=lwd,lty=lty)}
		if(add.line.smooth==T) {lines(lowess(bin.summary_noflat[,3][bin.sequence_noflat],as.numeric(matrix(fitindices_noflat[8,]))[bin.sequence_noflat]),lwd=lwd,lty=lty)}
	}
	if(plot.rmsea==T) {
	if(rmsea_err==T) {plotCI(matrix(bin.summary[,3]),as.numeric(matrix(fitindices[23,])),uiw=as.numeric(matrix(fitindices[25,])),liw=as.numeric(matrix(fitindices[24,])),
		pch=pch,gap=0,main="RMSEA vs Lag Dist",xlab="Lag Distances",ylab="Model RMSEA Values",cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,cex.main=cex.main)}
	if(rmsea_err==F) {plot(matrix(bin.summary[,3]),as.numeric(matrix(fitindices[23,])),pch=pch,main="RMSEA vs Lag Dist",xlab="Lag Distances",ylab="Model RMSEA Values",cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,cex.main=cex.main)}
		abline(h=0.05,lwd=2,lty=3)
		if(add.line.step==T) {lines(bin.summary[,3][bin.sequence],as.numeric(matrix(fitindices[23,]))[bin.sequence],lwd=lwd,lty=lty)}
		if(add.line.smooth==T) {lines(lowess(bin.summary_noflat[,3][bin.sequence_noflat],as.numeric(matrix(fitindices_noflat[23,]))[bin.sequence_noflat]),lwd=lwd,lty=lty)}
	}
	if(plot.srmr==T) {
	plot(matrix(bin.summary[,3]),as.numeric(matrix(fitindices[29,])),main="SRMR vs Lag dist",pch=pch,xlab="Lag Distances",ylab="Model SRMR Values",cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,cex.main=cex.main)
		abline(h=0.08,lwd=2,lty=3)
		if(add.line.step==T) {lines(bin.summary[,3][bin.sequence],as.numeric(matrix(fitindices[29,]))[bin.sequence],lwd=lwd,lty=lty)}
		if(add.line.smooth==T) {lines(lowess(bin.summary_noflat[,3][bin.sequence_noflat],as.numeric(matrix(fitindices_noflat[29,]))[bin.sequence_noflat]),lwd=lwd,lty=lty)}
	}
	par(mfrow=c(1,1),mar=c(4,4,2,2))
}

#	plotpath(spatial_model_results,path.type="directed",add.line="none",rmsea_err=T,pch=16,lwd=2,lty=1)
#   
#	spatial_model_results = object output by run.Models
#	path.type= options for selecting which paths to plot
#		"directed" = only directed paths plotted
#		"undirected" = only undirected correlations plotted
#		"both" = all paths plotted
#		"user" = allows user to specify particular paths and a particular order for plotting
#			argument selectpath must also be provided with path.type="user"
#	selectpath = required when path.type="user"
#		usage is selectpath==c(5,18,16,23,29) where values refer to path numbers
#		path numbers can be obtained using spatial_model_results[2]
#	add.lines. options for plotting a fit line 
#		"none" indicates no line
#		"step" plots straight line segments between points
#		"smooth" plots a smoothed curve fit using function lowess
#	add.error. should standard error bars be added for each path coefficient
#	pcut - p-value cutoff above which points with non significant p-values are shaded grey
#		set pcut=1 to have all points black
#	pch, lwd, lty options for formatting points and fit lines.


#' @export
#' @import gplots
plotpath<-function(spatial_model_results,path.type="directed",selectpath="none selected",add.line="none",add.error=T,pcut=0.05,pch=16,lwd=2,lty=1,cex.main=1.2) {
	path.type.dir<-ifelse(path.type=="directed",T,F)
	path.type.undir<-ifelse(path.type=="undirected",T,F)
	path.type.both<-ifelse(path.type=="both",T,F)
	path.type.user<-ifelse(path.type=="user",T,F)
	if(path.type.dir==T) {selected.path<-grep(" ~ ",spatial_model_results[[2]][,2])}
	if(path.type.undir==T) {selected.path<-grep(" ~~ ",spatial_model_results[[2]][,2])}
	if(path.type.both==T) {selected.path<-spatial_model_results[[2]][,1]}
	if(path.type.user==T) {selected.path<-selectpath}
	bin.summary<-data.frame(spatial_model_results[[11]])
	bin.summary_noflat<-bin.summary[-length(bin.summary[,1]),]
	bin.sequence<-order(bin.summary[,3])
	bin.sequence_noflat<-bin.sequence[c(2:length(bin.sequence))]
	path.name<-spatial_model_results[[2]][selected.path,2]
	est<-spatial_model_results[[3]][selected.path,]
	se<-spatial_model_results[[4]][selected.path,]
	pval<-spatial_model_results[[5]][selected.path,]
	add.line.step<-ifelse(add.line=="step",T,F)
	add.line.smooth<-ifelse(add.line=="smooth",T,F)
	par(mfrow=n2mfrow(length(selected.path)),mar=c(2.5,2.5,1.5,1))
	for(i in 1:length(selected.path)){
		if(add.error==T & sum(se[i,])>0) {
			plotCI(matrix(bin.summary[,3]),as.numeric(matrix(est[i,])),
				uiw=as.numeric(matrix(se[i,])),liw=as.numeric(matrix(se[i,])),
				pch=pch,main=path.name[i],cex=ifelse(is.na(pval[i,])==T,1,ifelse(pval[i,]<pcut,1.2,1)),
				col=ifelse(is.na(pval[i,])==T,1,ifelse(pval[i,]<pcut,201,24)),lwd=2,cex.main=cex.main,gap=0)}
		if(add.error==T & sum(se[i,])==0) {
			plot(matrix(bin.summary[,3]),as.numeric(matrix(est[i,])),
				pch=pch,main=path.name[i],cex=ifelse(is.na(pval[i,])==T,1,ifelse(pval[i,]<pcut,1.2,1)),
				col=ifelse(is.na(pval[i,])==T,1,ifelse(pval[i,]<pcut,201,24)),lwd=2,cex.main=cex.main)}
		if(add.error==F) {
			plot(matrix(bin.summary[,3]),as.numeric(matrix(est[i,])),
				pch=pch,main=path.name[i],cex=ifelse(is.na(pval[i,])==T,1,ifelse(pval[i,]<pcut,1.2,1)),
				col=ifelse(is.na(pval[i,])==T,1,ifelse(pval[i,]<pcut,201,24)),lwd=2,cex.main=cex.main)}
		abline(h=0,lwd=1,lty=1)
		est.sequence<-as.numeric(matrix(est[i,]))[bin.sequence]
		est.sequence_noflat<-as.numeric(matrix(est[i,-length(bin.summary[,1])]))[bin.sequence_noflat]
		if(add.line.step==T) {lines(bin.summary[,3][bin.sequence],est.sequence,lwd=lwd,lty=lty)}
		if(add.line.smooth==T) {lines(lowess(bin.summary_noflat[,3][bin.sequence_noflat],est.sequence_noflat),lwd=lwd,lty=lty)}
}
	par(mfrow=c(1,1))
}


#	gam.path(spatial_model_results,path.type="directed",selectpath="none selected",
#		plot.points=T,se.plot=T,lwd.pred=2,lty.pred=1,lwd.se=2,lty.se=3,cex=1,cex.axis=1,cex.lab=1
#		,xlab="Lag Distance",ylab="Unst. Path Coeff.",yaxt="s",xaxt="s")
#
#	library "mgcv" is required
   
#	This function fits and prints some of the gam results for each path in the model as well as producing figures
   
#	spatial_model_results = object output by summary.models
#	path.type= options for selecting which paths to plot
#		"directed" = only directed paths plotted
#		"undirected" = only undirected correlations plotted
#		"both" = all paths plotted
#		"user" = allows user to specify particular paths and a particular order for plotting
#			argument selectpath must also be provided with path.type="user"
#		selectpath = required when path.type="user"
#		usage is selectpath==c(5,18,16,23,29) where values refer to path numbers
#		path numbers can be obtained using spatial_model_results[2]
#	plot.points = should points be plotted on figure
#	se.plot = should standard errors for the prediction lines be printed
#	additional arguments formatting figure
#' @export
#' @import mgcv gplots
gam.path<-function(spatial_model_results,path.type="directed",selectpath="none selected",plot.points=T,se.plot=T,lwd.pred=2,lty.pred=1,lwd.se=2,lty.se=3,cex=1,cex.axis=1,cex.lab=1,xlab="Lag Distance",ylab="Unst. Path Coeff.",yaxt="s",xaxt="s") {
	path.type.dir<-ifelse(path.type=="directed",T,F)
	path.type.undir<-ifelse(path.type=="undirected",T,F)
	path.type.both<-ifelse(path.type=="both",T,F)
	path.type.user<-ifelse(path.type=="user",T,F)
	if(path.type.dir==T) {selected.path<-grep(" ~ ",spatial_model_results[[2]][,2])}
	if(path.type.undir==T) {selected.path<-grep(" ~~ ",spatial_model_results[[2]][,2])}
	if(path.type.both==T) {selected.path<-spatial_model_results[[2]][,1]}
	if(path.type.user==T) {selected.path<-selectpath}
	bin.summary<-data.frame(spatial_model_results[[11]])
	bin.summary_noflat<-bin.summary[-length(bin.summary[,1]),]
	bin.sequence<-order(bin.summary[,3])
	bin.sequence_noflat<-bin.sequence[c(2:length(bin.sequence))]
	path.name<-spatial_model_results[[2]][selected.path,2]
	est_noflat<-spatial_model_results[[3]][selected.path,c(1:length(bin.summary[,1])-1)]
	spatial_predictor<-matrix(bin.summary_noflat[,3])
	gam.dev.exp<-matrix(numeric(1*length(path.name)),ncol=1)
	gam.s.edf<-matrix(numeric(1*length(path.name)),ncol=1)
	gam.s.Ref.df<-matrix(numeric(1*length(path.name)),ncol=1)
	gam.s.F<-matrix(numeric(1*length(path.name)),ncol=1)
	gam.s.pvalue<-matrix(numeric(1*length(path.name)),ncol=1)
	max.path.value<-matrix(numeric(1*length(path.name)),ncol=1)
	par(mfrow=n2mfrow(length(selected.path)))
	pd <- data.frame(spatial_predictor=seq(min(spatial_predictor),max(spatial_predictor),length=100))
	for(i in 1:length(selected.path)){
		gammodel<-gam(as.numeric(matrix(est_noflat[i,]))~s(spatial_predictor))
		sum_model<-summary(gammodel)
		gam.dev.exp[i]<-sum_model$dev.exp
		gam.s.edf[i]<-sum_model[24]$s.table[1]
		gam.s.Ref.df[i]<-sum_model[24]$s.table[2]
		gam.s.F[i]<-sum_model[24]$s.table[3]
		gam.s.pvalue[i]<-sum_model[24]$s.table[4]
		pv <- predict(gammodel,newdata=pd,type="link",se=TRUE)
		max.path.value[i]<-pd[pv$fit==max(abs(pv$fit))*sign(pv$fit),1]#spatial predictor value gam with largest absolute value
		plot(spatial_predictor,as.numeric(matrix(est_noflat[i,])),type="n",main=path.name[i],xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,cex.axis=cex.axis,cex.lab=cex.lab)
		if(plot.points==T) {points(spatial_predictor,as.numeric(matrix(est_noflat[i,])),pch=16,cex=cex)}
		se_upp<-pv$fit+1.96*pv$se
		se_low<-pv$fit-1.96*pv$se
		abline(h=0,lwd=1,lty=1)
		lines(pd$spatial_predictor,pv$fit,lwd=lwd.pred,lty=lty.pred)
		if(se.plot==T) {
			lines(pd$spatial_predictor,se_upp,lwd=lwd.se,lty=lty.se)
			lines(pd$spatial_predictor,se_low,lwd=lwd.se,lty=lty.se)
		}}
	par(mfrow=c(1,1))
	gamresults<-data.frame(path.name,round(gam.dev.exp,digits=3),round(gam.s.edf,digits=3),round(gam.s.Ref.df,digits=3),round(gam.s.F,digits=3),round(gam.s.pvalue,digits=3),round(max.path.value,digits=3))
	colnames(gamresults)<-c("path.name","gam.dev.exp","gam.s.edf","gam.s.Ref.df","gam.s.F","gam.s.pvalue","max.path.value")
	return(gamresults)
}


# avg.modindices() extracts modification indices for all models in the object produced by runModels() 
# 		and summarizes the mod indices by taking the mean for each possible additional path
#		across all lag distance bins. The flat model is not included in these calculations
#	 modcut eliminates printing of average MI values below the cutoff. The default is 4
# 
#' @export
avg.modindices<-function(spatial_model_results,modcut=4){
	mod.indices<-spatial_model_results[[10]][,names(spatial_model_results[[10]])!="binflat"]
	rownames(mod.indices)<-spatial_model_results[[9]][,1]
	mod.indices<-na.omit(mod.indices)
	mean.modindices<-rowSums(mod.indices)/length(mod.indices)
	return(round(mean.modindices[mean.modindices>=modcut],4))
}
