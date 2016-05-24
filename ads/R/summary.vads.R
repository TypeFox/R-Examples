summary.vads<-function(object,...) {
	UseMethod("summary.vads")
}

summary.vads.dval<-function (object,...) {
    res<-list(spp=summary(eval(object$call$p)))
	res$gsize<-c(diff(res$spp$window$xrange)/object$call$nx,diff(res$spp$window$yrange)/object$call$ny)
	res$ldens<-data.frame(apply(object$dval,2,summary))
	names(res$ldens)<-as.character(object$r)
	class(res)<-"summary.dval"
	return(res)
}

print.summary.dval<-function(x,...) {
	cat(paste("Multiscale first-order local density:\n"))
	print(x$spp)
	cat(paste("grid size:",x$gsize[1],"X",x$gsize[2],"\n"))
	cat("local density:\n")
	print(t(x$ldens))
}

summary.vads.kval<-function (object,...) {
    res<-list(spp=summary(eval(object$call$p)))
	res$lneig<-data.frame(apply(object$nval,2,summary))
	names(res$lneig)<-as.character(object$r)
	class(res)<-"summary.kval"
	return(res)
}

print.summary.kval<-function(x,...) {
	cat(paste("Multiscale second-order local neighbour density:\n"))
	print(x$spp)
	cat("local neighbour density:\n")
	print(t(x$lneig))
}

summary.vads.k12val<-function (object,...) {
    res<-list(spp=summary(eval(object$call$p)))
	res$marks<-object$marks
	res$lneig<-data.frame(apply(object$n12val,2,summary))
	names(res$lneig)<-as.character(object$r)
	class(res)<-"summary.k12val"
	return(res)
}

print.summary.k12val<-function(x,...) {
	cat(paste("Multiscale bivariate second-order local neighbour density:\n"))
	print(x$spp)
	cat("mark 1:",x$marks[1],"\n")
	cat("mark 2:",x$marks[2],"\n")
	cat("bivariate local neighbour density:\n")
	print(t(x$lneig))
}
