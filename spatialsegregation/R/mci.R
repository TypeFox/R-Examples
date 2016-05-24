# mci.R
# 
# Mean compositional information for multitype spatial point pattern and graphs
# 
# Author: Tuomas Rajala <tuomas.rajala@jyu.fi>
###############################################################################

# mciF

mciF<-function(X, r=NULL, target=NULL, ...)
{
	# check that X is ppp-object
	verifyclass(X, "ppp")
	if(length(levels(X$marks))<2) stop("Use only on a multitype point pattern (data.frame-marks not yet supported). ")
	# if no target given, calculate for all types
	funtype<-"Mean compositional information"
	if(is.null(target))
		targeti <- 0
	# else convert to an integer
	else
	{
		if(!is.factor(X$marks))warning("Marks of X are not in factor form. Transforming.")
		X$marks<-as.factor(X$marks)
		targeti<- which( levels(X$marks)  == target)
#		targeti<-which( union(X$marks, NULL) == target)
		if(length(targeti)!=1) stop("Target type not one of pattern types.")
		funtype<-paste(funtype, "of type", target)
	}
	
	# use the main calc function
	res<-segregationFun(X=X, fun="mci", r, funpars=targeti, ...)
	
	
	
	
	# theoretical values in CSR
	theo<-rep(0, length(res$parvec))
			
	# create the fv-object
	mci.final<-fv(data.frame(theo=theo,par=res$parvec), 
			argu="par",
			alim=range(res$parvec),
			ylab=substitute(MCI, NULL),
			desc=c("CSR values","Parameter values"),
			valu="theo",
			fmla=".~par",
			unitname=res$unitname,
			fname=funtype
	)
#	return(res)
	omarks<-order(union(X$marks[res$included],NULL)) # right order of marks
	# center on CSR
	sum0<-summary(X[res$included])
	reallythere<-sum0$marks[,3]>0
	lvec<-sum0$marks[reallythere,3][omarks]
	wmean<-function(x, w) sum(w*x)/sum(w)
	MCIMean<-apply(res$v, 1, wmean, w=lvec)
	
	if(targeti==0)
	{
		# the values from calculation
		tw<-res$v
		
		# set the names right, and don't forget to check border correction inclusion (might drop some types off)
		colnames(tw)<-union(X$marks[res$included],NULL)
		mci.final<-bind.fv(x=mci.final,
				y=tw,
				desc=paste("Typewise",funtype,"for type",colnames(tw)),
				labl=colnames(tw)
		)
		
		mci.final<-bind.fv(x=mci.final,
				y=data.frame("MCIMean"=MCIMean),
				desc=paste("Mean MCI over types"),
				labl="MCIMean",
				preferred="MCIMean"
		)		 
	}
	
	# if target type given add the values for the target type
	else
	{
		mci.final<-bind.fv(x=mci.final,
				y=data.frame("MCI"=res$v[,1]),
				desc=funtype,
				labl="MCI",
				preferred="MCI"
		)
	}
	
	# attach the frequencies too
	attr(mci.final,"frequencies")<-freqs(X[res$included])
	
	# and some notes
	attr(mci.final,"neighbourhoodType")<-res$ntype
	attr(mci.final,"note")<-res$note
	
	# return 
	mci.final
}

###############################################################################
###############################################################################
## R-version: Geometric neighbourhood only, compute per species and then weighted mean
#mciFR<-function(pp,  r=seq(0,0.3, length=50), ...)
#{
#	require(spatgraphs)
#	sum0<-summary(pp)
#	truelythere<-sum0$marks[,3]>0
#	mvec<-levels(pp$marks)[truelythere]#union(pp$marks, NULL)
#	
#	CSR<-pvec<-NULL
#	lvec<-sum0$marks[truelythere,3]
#	# estimate from data
#	Ivec<-NULL
#	for(r1 in r)
#	{
#		# CSR value
#		pvec<-exp(-pi*lvec*r1^2)
#		CSR<-c(CSR,-sum((1-pvec)*log(1-pvec)+pvec*log(pvec)))
#		
#		# data
#		g<-spatgraph(pp, type="geo", par=r1, ...)
#		Ivecr<-NULL
#		
#		for(mark in mvec) # mean over one mark
#		{
#			I1<-n<-0
#			for(i in (1:pp$n)[pp$marks==mark]) # go through points of one mark
#			{
#				n<-n+1
#				I1i<-0
#				for(z in 1:length(mvec)) # check how many types present
#				{
#					m<-mvec[z]
#					Fi<-1*(m%in%pp$marks[g$edges[[i]]])
#					I1i<-I1i+Fi*log(1-pvec[z]) + (1-Fi)*log(pvec[z])
#				}
#				I1<-I1-I1i#c(I1,I1i)
##				cat(i,":",paste(Fi),"\n",sep="")
#			}
##			cat("\n")
#			Ivecm<-I1/n
#			Ivecr<-rbind(Ivecr, Ivecm)
#		}
#		Ivec<-cbind(Ivec, Ivecr)
#	}
#	
#	rownames(Ivec)<-mvec
#	
#	# weighted mean
#	wmean<-function(x, w) sum(w*x)/sum(w)
#	
#	Ivec<-t(apply(Ivec, 1, function(x)x-CSR))
#	v<-apply(Ivec, 2, wmean, w=lvec)
#	list(v=v, CSR=CSR, Ivec=Ivec, r=r, pvec=pvec, lvec=lvec)
#}
#
##eof
