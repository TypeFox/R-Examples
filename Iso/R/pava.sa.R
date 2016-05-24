pava.sa <- function(y,w=NULL,decreasing=FALSE,long.out=FALSE,stepfun=FALSE)
{
#
# Function 'pava.sa' (stand-alone pava). To perform isotonic
# regression for a simple (increasing) linear ordering using the ``pool
# adjacent violators algorithm''.  This version is programmed in raw
# R; i.e.  it does not invoke dynamically loaded fortran.  If
# long.out is TRUE then the result returned consists of a list
# containing the fitted values, the final weights, and a set of
# indices `tr', made up of the smallest index in each level set, which
# thus keeps track of the level sets.  If in addition stepfun is TRUE,
# then the step function represention of the isotonic regression is
# added to the forgoing list. If stepfun is TRUE and long.out is FALSE
# then only the stepfunction representation is returned.  If stepfun
# and long.out are both FALSE then only the fitted values are
# returned.
# 

if(decreasing) y <- rev(y)
n <- length(y)
if(is.null(w))
	w <- rep(1,n)
else if(decreasing) w <- rev(w)
r <- rep(1,n)
repeat {
	stble <- TRUE
	i <- 1
	while(i < n) {
		if(y[i] > y[i+1]) {
			stble <- FALSE
			www <- w[i] + w[i+1]
			ttt <- (w[i]*y[i] + w[i+1]*y[i+1])/www
			y[i+1] <- ttt
			w[i+1] <- www
			y <- y[-i]
			w <- w[-i]
			r[i+1] <- r[i] + r[i+1]
			r <- r[-i]
			n <- n-1
		}
		i <- i+1
	}
if(stble) break
}
y  <- rep(y,r)
if(decreasing) y <- rev(y)
if(long.out | stepfun) {
	if(decreasing) r <- rev(r)
	tr <- rep(tapply(1:length(y),rep(1:length(r),r),min),r)
}
if(long.out) {
	if(decreasing) w <- rev(w)
	w <- rep(w,r)
	lout <- list(y=y,w=w,tr=tr)
}
if(stepfun) {
	knots <- 1+which(diff(tr)!=0)
	y0    <- c(y[1],y[knots])
	h     <- stepfun(knots,y0)
}
ntype <- 1+sum(c(long.out,stepfun)*(1:2))
switch(ntype,y,lout,h,c(lout,list(h=h)))
}
