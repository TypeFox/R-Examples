ica.resolution <-
function(D, k)
{
	ICA.res <- JADE(D,k, center=F)
	w.sign <- diag(sign(as.matrix(ICA.res$S)[apply(as.matrix(abs(ICA.res$S)),2, which.max),]))
	w.sign[w.sign==0] <- -1
	ICA.res$S <- sweep(as.matrix(ICA.res$S),2,w.sign,"*")
	ICA.res
}
