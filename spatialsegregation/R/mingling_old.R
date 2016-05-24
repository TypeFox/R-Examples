## mingling.R
## 
## Spatial mingling index for multitype spatial point pattern and graphs
## 
## Author: Tuomas Rajala <tarajala@maths.jyu.fi>
################################################################################
#
## minglingF
#
#
#minglingF.old<-function(pp, parvec=1:20, graph_type="knn", type=NULL, ratio=FALSE, ...)
##Mingling index for geometric graph with various range-parameters
##If target type is not given (type=NULL) take over all types mean
##type as integer 1,...,S=number of types
##If ratio=T  use M = (1 - M)/(lambda_t/lambda) instead.
#{
#	
#	note2<-paste("Mingling index for type",type,".")
#	if(is.null(type))
#	{
#		type<-0
#		note2<-"Mean mingling index over all types."
#	}
#	res<-segregationFun(pp, fpar=c(as.numeric(type),as.numeric(ratio)), graph_type=graph_type, graph_parvec=parvec, funtype=1, ...)
#	poisson<-1
#	sum0<-summary(pp)
#	if(!ratio) poisson<-mean(1-sum0$marks[,3]/sum0$int)
#
#	if(pp$markformat=="none")pp$markformat<-"vector" # TODO: fix this bug
#	f<-freqs(pp[res$included])
#	w<-f/sum(f)
#	
#	segfcl(list(M=apply(unname(res$v),2,mean),sd=apply((res$v),2,sd), 
#			 typewise=(res$v), Mw=apply(unname(res$v),2,weighted.mean,w=w),gtype=graph_type, 
#			 frequencies=f,
#			 par=res$parvec,
#			 note=res$note,note2=note2, poisson=poisson))
#}
#
################################################################################
#
#
################################################################################
#mingling_index.old<-function(pp, graph_type="knn", graph_par=4, type=NULL, ...)
##if type=NULL compute the arithmetic mean over all types
##type as integer 1,...,S=number of types
#{
#	if(is.null(type)){type<-0;nvec<-paste("Mean Mingling index over all types.")}
#	else nvec<-paste("Mingling index for type",type)
#	M0<-minglingF(pp=pp, graph_type=graph_type, parvec=graph_par[1], type=type, ...)
#	M<-M0$M	
#	names(M)<-nvec
#	M
#}
#
##eof
