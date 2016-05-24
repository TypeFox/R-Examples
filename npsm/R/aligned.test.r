aligned.test<-function(x,y,g,scores=Rfit::wscores) {

############################################################
#
#  input: 
#    x : N by p design matrix of covariates
#    y : N by 1 response vector
#    g : N by 1 vector representing treatment assignment
#    scores (default = Wilcoxon) : which scores to use
#
#  output:
#    statistic : the aligned rank test statistic
#    p.value : the p-value based on the asymptotic distribution
#
############################################################

	k<-length(unique(g))
	N<-length(y)

	w<-matrix(model.matrix(~as.factor(g))[,2:k],ncol=(k-1))

	fitr<-rfit(y~x)
	A<-qr.qty(fitr$qrx1,w)[(fitr$qrx1$rank+1):N,]
	apai<-chol2inv(chol(crossprod(A)))
	S2<-crossprod(w,getScores(scores,resid(fitr)))

	tstat<-t(S2)%*%apai%*%S2
	pval<-pchisq(tstat,k-1,lower.tail=FALSE)

	res<-list(statistic=tstat, p.value=pval)

	class(res)<-'rank.test'
	res


}
