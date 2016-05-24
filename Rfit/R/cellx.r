cellx<-function(X) {
	X<-pasteColsRfit(t(X))
	W<-suppressWarnings(model.matrix(~0+X))
	matrix(W,ncol=ncol(W))
}
