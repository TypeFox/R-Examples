## Spatial Shannon index functional
##
##
## Author: Tuomas Rajala <tarajala@maths.jyu.fi>
################################################################################
#
#shannonF.old<-function(pp, parvec=1:20, graph_type="knn", v2=FALSE, ...)
##Shannon index for graphs, with possibly various range-parameters
#{
#	sum0<-summary(pp)
#	note<-""
#	parvec1<-parvec
#		
#	res<-segregationFun(pp=pp, fpar=ifelse(v2,1,0), graph_type=graph_type, graph_parvec=parvec, funtype=2, ...) # calculates only the pii_tau-vector
#	
#	#lambs<-summary(pp[as.logical(res$included)])$marks[,3]
#	#lambs<-lambs[lambs>0]
#	#lamb0<-sum(lambs)
#	#eglobal<-sum(lambs/lamb0*log(lambs/lamb0,length(lambs)))
#	eglobal<--shannon_index(pp,spatial=FALSE)
#	f<-function(pvec) # to calculate the -pi_tau*log(pi_tau)-sum
#	{
#		a<-0
#		S<-length(pvec)
#		for(i in 1:S)
#			a<-a + ifelse(pvec[i]>0, pvec[i]*log(pvec[i],base=S), 0)
#		
#		1-a/eglobal
#	}
#	if(v2)
#	{
#		H<-apply(-res$v,2,mean)
#		typewise<--res$v
#		note2<-"Spatial Shannon index, v2"
#	}
#	else
#	{
#		H<-unname(apply(res$v,2,f))
#		typewise<-res$v
#		note2<-"Spatial Shannon index"
#	}	
#	
#	segfcl(list(H=H,typewise=typewise, global=-eglobal, par=res$parvec,note=res$note, note2=note2, gtype=graph_type, poisson=0))
#}
#
#
##Shannon index
#
#shannon_index.old<-function(pp, graph_type="knn", graph_par=4, spatial=TRUE, ...)
#{
#	if(!spatial)
#	{			 #the traditional aspatial Shannon index
#		sum0<-summary(pp)
#		ints<-sum0$marks[,3]
#		m<-union(pp$marks,NULL)
#		pii<-ints/sum(ints)
#		H<- (-1)*sum( pii * log(pii, length(m)),na.rm=TRUE)
#		names(H)<-"Non-spatial Shannon index"
#	}
#	if(spatial)
#	{			   #spatial Simpson index for a set of edgelists	
#		H<-shannonF(pp=pp, parvec=graph_par, graph_type=graph_type, ...)$H
#		names(H)<-"Spatial Shannon index"
#	}
#	H
#}
