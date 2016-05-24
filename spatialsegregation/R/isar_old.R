## isar.R
## 
## ISAR function for multitype spatial point pattern and graphs
## see: T. Wiegand et al.: How individual species structure diversity in tropical forests. PNAS, November 16 2007
#
## Author: Tuomas Rajala <tarajala@maths.jyu.fi>
## last changed: 060909
################################################################################
## 
#
#
#isarF.old<-function(pp, parvec=1:20, r=parvec, graph_type="knn", type=NULL, v2=FALSE, ...)
##ISAR for various graphs
##If target type is not given (type=NULL) take over all types mean
##type as integer 1,...,S=number of types
#
#{
#	
#	note2<-paste("ISAR for type",type)
#	if(is.null(type))
#	{
#		type<-0
#		note2<-"ISAR over all types."
#	}
#	note3<-NULL
#	if(v2)note3<-"Degree weighted values, E_x [I(x)/deg(x)]."
#
##   main function call
#	res<-segregationFun(pp, fpar=c(type,as.integer(v2)), graph_type=graph_type, graph_parvec=parvec, funtype=4, ...)
#	
##   the poisson values
#	if(graph_type=="geometric")mdeg<-function(l,k)pi*l*k^2
#	if(graph_type=="knn")mdeg<-function(l,k)k
#	if(graph_type=="delauney")mdeg<-function(l,k) 6
#	if(graph_type=="gabriel")mdeg<-function(l,k) 4
#	poisson<-NULL
#	sum0<-summary(pp)
#	l<-sum0$marks[,3]
#	for(i in 1:length(parvec)) poisson<-c(poisson, sum(1-(1-l/sum(l))^mdeg(sum(l),parvec[i]) )/ifelse(v2,mdeg(sum(l),parvec[i]),1) )
#	
#	#other
#	if(pp$markformat=="none")pp$markformat<-"vector" # TODO: fix this bug
#	f<-freqs(pp[res$included])
#	w<-f/sum(f)
#	if(type!=0)w<-w[which(rownames(sum0$marks)==type)]
#	
#	# old class, swiched to fv-class from spatstat for compatibility
#	#segfcl(list(I=apply(res$v,2,mean,na.rm=TRUE),sd=apply(res$v,2,sd,na.rm=TRUE),typewise=res$v,
#	#				Iw=apply(res$v,2,weighted.mean,w=w,na.rm=TRUE), frequencies=f, 
#	#				gtype=graph_type,par=res$parvec,note=res$note,note2=note2,note3=note3,poisson=poisson))
#	
#	tw<-t(res$v)
#	colnames(tw)<-paste("type",1:length(f),sep="")
#	isar<-data.frame(theo=poisson,		 
#					typewise=tw,
#					I=apply(res$v,2,mean,na.rm=TRUE), #Iw=apply(res$v,2,weighted.mean,w=w,na.rm=TRUE), 
#					 par=res$parvec )
#	attr(isar,"frequencies")<-f
#	attr(isar,"graphType")<-graph_type
#	attr(isar,"note1")<-res$note
#	attr(isar,"note2")<-note2
#	attr(isar,"note3")<-note3
#	fv(x=isar, argu="par", ylab=substitute(ISAR, NULL), valu="I", fmla=".~par", alim=range(parvec),
#			labl=names(isar), desc=NULL, unitname=NULL, fname="ISAR")
#}
#
#
#
################################################################################
#
################################################################################
#isar_index.old<-function(pp, graph_type="knn", graph_par=4, type=NULL, ...)
##if type=NULL compute the arithmetic mean over all types
##type as integer 1,...,S=number of types
#{
#	if(is.null(type)){type<-0;nvec<-paste("Mean ISAR over all types.")}
#	else nvec<-paste("ISAR for type",type)
#	I0<-isarF(pp=pp, graph_type=graph_type, parvec=graph_par[1], type=type, ...)
#	I<-I0$I
#	names(I)<-nvec
#	I
#}