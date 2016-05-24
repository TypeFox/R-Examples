# 
# Author: Ian
###############################################################################

summarylm<-function(object,correlation=FALSE,symbolic.cor = FALSE,white.adjust=FALSE,...){
	if(is.logical(white.adjust) && !white.adjust)
		return(stats::summary.lm(object,correlation,symbolic.cor,...))
	if(is.logical(white.adjust))
		white.adjust='hc3'
	ans<-list()
	ans$Estimate<-coef(object)
	ans[['Std..Error']]<-sqrt(diag(hccm(object,type=white.adjust)))
	ans[['t.value']]<-ans$Estimate/ans[['Std..Error']]
	ans[['p-value']]<-2*pt(abs(ans[['t.value']]),object$df.residual,lower.tail=FALSE)
	ans<-as.data.frame(ans)
	result<-list('Summary Table'=ans)
	if(correlation){
		result[['correlation']]<-hccm(object,type=white.adjust)/outer(ans[['Std..Error']],ans[['Std..Error']])
	}
	result
}
