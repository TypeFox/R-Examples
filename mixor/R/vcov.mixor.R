vcov.mixor <-
function(object, ...) {
    dimnames(object$varcov)[[1]]<-dimnames(object$varcov)[[2]]<-dimnames(object$Model)[[1]]
	object$varcov
}
