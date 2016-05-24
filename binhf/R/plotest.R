`plotest` <-
function(l, plot.it=FALSE,verbose=FALSE){

x <- l$x
truep <- l$truep

if (plot.it==TRUE)	{
	plot(x, truep, type="l", ylim=c(0,1))

	lines(x, l$fhat, col=2)
	lines(x, l$fhata, col=3)
	lines(x, l$fhatf, col=4)
	}

hfn <- norm(truep, l$fhat)
an <- norm(truep, l$fhata)
fn <- norm(truep, l$fhatf)

if (verbose==TRUE)	{
	cat("HF norm is ", hfn, "\n")
	cat("Ans norm is ", an, "\n")
	cat("Free norm is ", fn, "\n")
	}

c(hfn, an, fn)

}

