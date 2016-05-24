#############################################################
#	All functions in this file are copyrighted
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

"rimtrf" <-
function(x,n=nrow(x),intch=1,tau=1.e-6) {
np <- ncol(x)
mdx <- nrow(x)
if (missing(x)) x <- matrix(single(1),mdx,np)
k <- integer(1)
sf <- single(np)
sg <- single(np)
sh <- single(np)
ip <- integer(np)
f.res <- .Fortran("rimtrf",
x=to.single(x),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
intch=to.integer(intch),
tau=to.single(tau),
k=to.integer(k),
sf=to.single(sf),
sg=to.single(sg),
sh=to.single(sh),
ip=to.integer(ip))
list(x=f.res$x,k=f.res$k,sf=f.res$sf,sg=f.res$sg,sh=f.res$sh,ip=f.res$ip)
}
