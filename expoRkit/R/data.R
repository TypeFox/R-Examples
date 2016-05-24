##' Australia Economic Model data, 1968-68.
##'
##' This sparse matrix of order 2,529 and with 90,158 non-zero
##' elements was used in Sidje (1998) to illustrate the application of
##' Expokit.
##' 
##' @name orani
##' @docType data
##' @keywords data
##' @references Sidje, R. B. (1998) Expokit. Software Package for Computing Matrix
##' Exponentials. ACM Trans. Math. Softw. 24(1), 130-156.
##' @seealso \code{\link{expv}}
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @examples
##' data(orani)  ## Load the data as a 'dgCMatrix' (CCS format)
##' v <- rep(1, 2529)
##' ### Solving a system of 2529 coupled linear differential equations
##' system.time(wCCS <- expv(orani, v = v, t = 10))
##' oraniCOO <- as(orani, "TsparseMatrix")  ## Coerce to COO format
##' ### In this case, the COO format gives a slight increase in
##' ### computational time as reported in Sidje (1998).
##' system.time(wCOO <- expv(oraniCOO, v = v, t = 10))
##' 
##' print(cbind(wCCS[1:5], wCOO[1:5]), digits = 14)
NULL
