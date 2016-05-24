func1 <- function(b,rl1,rl,delta,ni,maxiter)
{
	cb <- abs(rl1-rl)
	ca <- 0
	ca <- sum(delta*delta)
	b <- b+delta
	ni <- ni+1
	result <- list(ca=ca,cb=cb,b=b,ni=ni)
}