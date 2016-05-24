nonzero.coef <-
function(object, s) {
	estimates<-coef(object,s=s)
	keep<-estimates[estimates!=0]
	names(keep)<-dimnames(estimates)[[2]][estimates!=0]
	keep
	}

