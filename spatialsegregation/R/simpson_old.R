#
################################################################################
## Spatial Simpson index functional
##
##
## Author: Tuomas Rajala <tarajala@maths.jyu.fi>
################################################################################
#
#simpsonF.old<-function(pp,parvec=1:20, graph_type="knn", ...)
##Simpson index for graphs, with possibly various range-parameters
##notice that if 'relative'=T we take the normalised range which eases comparisons
#{
#	res<-segregationFun(pp, fpar=NULL, graph_type=graph_type, graph_parvec=parvec, funtype=3, ...)
#	note2<-"Spatial Simpson index"	
#	segfcl(list(S=1-unname(colSums(res$v)), typewise=res$v, par=res$parvec,gtype=graph_type,note=res$note, note2=note2, poisson=simpson_index(pp,spatial=FALSE)))
#}
#
##Simpson index
#
#simpson_index.old<-function(pp, graph_type="knn", graph_par=4, spatial=TRUE, ...)
#{
#	
#	if(!spatial)
#	{			#the traditional aspatial Simpson index
#		sum0<-summary(pp)
#		ints<-sum0$marks[,3]
#		m<-union(pp$marks,NULL)
#		pii<-ints/sum(ints)
#		n<-pii*pp$n
#		N<-pp$n
#		S<- 1 - sum( n*(n-1)  )/(N*(N-1))#Simpson index of diversity
#		names(S)<-"Non-spatial Simpson index"
#	}
#	if(spatial)
#	{			   #spatial Simpson index for a set of edgelists	
#		S<-simpsonF(pp, parvec=graph_par, graph_type=graph_type, ...)$S
#		#simpsonF2(pp,parvec=seq(0,3,length=30), graph_type="geometric", relative=T, toroidal=FALSE, doDists=T, dbg=F)
#		names(S)<-"Spatial Simpson index"
#	}
#	S
#}
