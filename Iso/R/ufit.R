ufit <- function(y,lmode=NULL,x=NULL,w=NULL,lc=TRUE, rc=TRUE,
                 type=c("raw","stepfun","both"))
{
#
# Function `ufit'.  Calculates the isotonic unimodal fit to a data
# sequence y, with mode at ``lmode''.  If lmode==NULL, then it determines
# the optimal (least squares) location for the mode, and the fit
# for this optimum value.  The optimum mode may, by virtue of the
# nature of the procedure, be taken to be one of the points x_i, i =
# 1, ..., n.  NOTE that the optimum will occur at one of the
# midpoints (x_i + x_{i+1})/2, i = 1, ..., n-1 AND at one of the two
# adjacent points, i.e. either at x_i or at x_{i+1}.  If x is null, x
# is taken to be an equispaced sequence on [0,1].
#
type <- match.arg(type)

n <- length(y)
if(is.null(w)) w <- rep(1,n)
if(is.null(x)) x <- seq(0,1,length=n)
if(is.null(lmode)) lmode <- -1
else {
	mode.save <- lmode
	k1 <- sum(x <= lmode)
	k2 <- n + 1 - sum(x >= lmode)
	lmode <- (k1+k2)/2
}

rslt <- .Fortran(
		"ufit",
		yk=as.double(y),
		wk=as.double(w),
		xmode=as.double(lmode),
		y=double(n),
		w=double(n),
		mse=double(1),
		y1=double(n),
		w1=double(n),
		y2=double(n),
		w2=double(n),
		ind=integer(n),
		kt=integer(n),
		n=as.integer(n),
		goof=logical(1),
		PACKAGE="Iso"
		)
if(rslt$goof) stop('Goof in unimode subroutine called by ufit subroutine.\n')
imode <- if(lmode < 0) rslt$xmode else lmode
lmode <- if(lmode < 0) x[imode] else mode.save
ys <- rslt$y
if(type%in%c("stepfun","both")) {
	kind <- 1+which(diff(ys)!=0)
	if(!(n%in%kind)) kind <- c(kind,n)
	y0    <- c(ys[1],ys[kind])
	h     <- stepfun(x[kind],y0)
}
i     <- floor(imode)
if(!lc) ys[i]   <- NA
if( (!rc) & (i < n) ) ys[i+1] <- NA
switch(type,raw=list(x=x,y=ys,mode=lmode,mse=rslt$mse),
            stepfun=h,
            both=list(x=x,y=ys,mode=lmode,mse=rslt$mse,h=h))
}
