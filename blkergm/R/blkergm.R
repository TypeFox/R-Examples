blkergm = function(formula,offset.coef=NULL,target.stats=NULL,eval.loglik=TRUE,estimate=c("MLE","MPLE"),control=control.ergm(),verbose=FALSE,...){
	model = ergm(formula=formula,
                 constraints=~.,
                 offset.coef=offset.coef,
                 target.stats=target.stats,
                 eval.loglik=eval.loglik,
                 estimate=estimate,
                 control=control,
                 verbose=verbose,...)
    model$call = match.call()
    class(model)=c("blkergm","ergm")
    model
}