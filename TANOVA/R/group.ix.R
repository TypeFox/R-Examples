###############################################################################
# TANOVA: test.R
# 
# TODO: Add comment
#
# Author: Weihong
# Mar 20, 2009, 2009
###############################################################################
group.ix<-function(f1,f2){
	n1<-nlevels(as.factor(f1))
	n2<-nlevels(as.factor(f2))
	ix<-list()
	k<-1
	for (i in 1:n1){
		for (j in 1:n2){
			ix[[k]]<-which(f1==i & f2==j)
			k<-k+1
		}
	}
	return (ix)
}