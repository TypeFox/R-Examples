#'@exportPattern "^[a-zA-Z][[:alpha:]]*"
#'@useDynLib wrswoR

#'@title Faster weighted sampling without replacement
#'@description \R's default sampling without replacement using
#'  \code{\link[base]{sample.int}} seems to require quadratic run time,
#'  e.g., when using weights drawn from a uniform distribution. For large
#'  sample sizes, this is too slow.  This package contains several
#'  alternative implementations.
#'@details Implementations are adapted from
#'  \url{http://stackoverflow.com/q/15113650/946850}.
#'
#'@name wrswoR-package
#'@aliases wrswoR-package wrswoR
#'@docType package
#'@author Kirill Müller
#'@references Efraimidis, Pavlos S., and Paul G. Spirakis. "Weighted
#'random sampling with a reservoir." \emph{Information Processing Letters} 97,
#'no. 5 (2006): 181-185.
#'
#'Wong, Chak-Kuen, and Malcolm C. Easton. "An efficient method for
#'weighted sampling without replacement." \emph{SIAM Journal on Computing} 9,
#'no. 1 (1980): 111-113.
#'
#'
#'@keywords package
#'@examples
#'sample_int_rej(100, 50, 1:100)
NULL



