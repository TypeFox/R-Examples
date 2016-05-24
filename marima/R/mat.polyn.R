##' @title pol.inv
##' 
##' @description Calculation of left inverse of matrix polynomial. The
##' leading term
##' is expected to be the (k by k) identity matrix. This is checked
##' and the proper leading unity terms are taken into account.
##'
##' phi = matrix polynomial coefficients = I,phi1,phi2,...phi(p)
##' dim(phi) = c(k,k,p+1)
##'    k   = dimension of coefficient matrices (k by k)
##'    L  = order of polynomial (length=1+L , including leading unity matrix)
##'
##' @param phi polynomium to invert
##' @param L   order of inverse polynomium
##'
##' @return  left inverse of phi of order L (L+1 terms including leading
##' unity matrix)
##'
##' @examples
##' set.seed(4711)
##' p2<-check.one(array(rnorm(32),dim=c(4,4,2)))
##' pi2<-pol.inv(p2,L=12)
##' short.form(pi2)
##'
##' @export

pol.inv <- function(phi, L) {
    aphi <- check.one(phi)
    k <- dim(aphi)[1]
    # if (dim(phi)[3]==dim(aphi)[3]) {up<-L+1}
    up <- L + 1
    ip <- dim(aphi)[3]
    polinv <- array(0, dim = c(k, k, up))
    polinv[, , 1] <- diag(k)
    if (up > 1) {
        if (ip > 1) {
            for (LL in 2:up) {
                for (J in 2:LL) {
                  if (J <= ip) {
                    help <- polinv[, , LL + 1 - J] %*% aphi[, , J]
                    polinv[, , LL] <- polinv[, , LL] - help
                  }
                }
            }
        }
    }
    return(polinv)
}


##' @title pol.mul
##' 
##' @description Calculation of product of two matrix
##' polynomials (arrays).
##'
##' If one or both leading unity matrices (of eta and theta) are
##' missing, they are (it is)
##' generated (and taken into account).
##'
##' @param eta first matrix polynomial
##' @param theta  second matrix olynomial
##' @param L order of output polynomial (length = L+1)
##'
##' @return  matrix polynomial product af eta and theta
##'
##' @examples
##' set.seed(4711)
##' p1<-check.one(matrix(rnorm(16),nrow=4))
##' p2<-check.one(array(rnorm(32),dim=c(4,4,2)))
##' p12<-pol.mul(p1,p2,L=(2+3))
##' short.form(p12)
##'
##' @export

pol.mul <- function(eta, theta, L) {
    eta <- check.one(eta)
    theta <- check.one(theta)
    if (L < 1) {
        L <- dim(eta)[3] + dim(theta)[3] - 1
    }
    L <- L + 1
    k <- dim(eta)[1]
    iet <- dim(eta)[3]
    ith <- dim(theta)[3]
    polmul <- array(0, dim = c(k, k, L))
    for (LL in 1:L) {
        polmul[, , LL] <- 0
        for (J in 1:LL) {
            if (J <= iet) {
                if (LL + 1 - J <= ith) {
                  help <- eta[, , J] %*% theta[, , LL + 1 - J]
                  polmul[, , LL] <- polmul[, , LL] + help
                }
            }
        }
    }
    return(polmul)
}

##' @title rand.shock
##' 
##' @description Calculation of random shock form for arma model
##'
##' @param ar.poly autoregressive matrix part of model
##' @param ma.poly moving average matrix part of model
##' @param L order of return polynomial  (length=L+1 including
##' leading unity matrix)
##'
##' @return random shock form of arma model up to order L (array(k,k,L+1))
##'
##' @examples
##' set.seed(4711)
##' p1  <- check.one(matrix(rnorm(16),nrow=4))
##' p2  <- check.one(array(rnorm(32),dim=c(4,4,2)))
##' randshock <- rand.shock(ar.poly=p1,ma.poly=p2,L=6)
##' short.form(randshock)
##'
##' @export

rand.shock <- function(ar.poly, ma.poly, L) {
    ar <- check.one(ar.poly)
    ma <- check.one(ma.poly)
    invar <- pol.inv(ar, (L + 1))
    rand.shock <- pol.mul(invar, ma, L)
    return(rand.shock)
}

##' @title inverse.form
##'
##' @description Calculation of inverse form for arma model
##'
##' @param ar.poly =autoregressive matrix part of model (array(k,k,ar-order)).
##' @param ma.poly =moving average matrix part of model (array(k,k,ma-order)).
##' @param L  =order of return polynomial (length=L+1 including leading
##' unity matrix).
##'
##' @return inverse form for arma model up to order L (array(k,k,L+1)).
##'
##' @examples
##' set.seed(4711)
##' p1  <- check.one(matrix(rnorm(16),nrow=4))
##' p2  <- check.one(array(rnorm(32),dim=c(4,4,2)))
##' inverse <- inverse.form(ar.poly=p1,ma.poly=p2,L=6)
##' short.form(inverse)
##'
##' @export

inverse.form <- function(ar.poly, ma.poly, L) {
    ar <- check.one(ar.poly)
    ma <- check.one(ma.poly)
    inma <- pol.inv(ma, (L + 1))
    inverse.form <- pol.mul(inma, ar, L)
    return(inverse.form)
}

##' @title check.one
##' 
##' @description Function to check and insert leading unity matrix
##' if NOT present.
##'
##' @param polyn (k,k,...) matrix polynomium with or without leading
##' unity matrix.
##'
##' @return polyn (array) with a leading unity matrix being
##' inserted if not present.
##'
##' @examples
##' set.seed(4711)
##' X<-array(rnorm(32),dim=c(4,4,2))
##' X<-check.one(X)
##' short.form(X)
##'
##' @export

check.one <- function(polyn = NULL) {

    if (is.null(polyn)) {
        stop("No input polynomial to 'check.one'. \n")
    }
    d <- dim(polyn)[1]
    dd <- dim(polyn)
    if (is.null(d)) {
        cat("Empty polynomial in call to check.one \n")
    }
    if (length(dim(polyn)) < 2) {
        cat("Argument not polynomial in call to check.one")
        cat(dd, " ", c(polyn), "\n")
    }
    pol2 <- polyn
    if (length(dim(polyn)) < 3) {
        pol2 <- array(polyn, dim = c(d, d, 1))
    }
    if (max(abs(pol2[, , 1] - diag(1, d))) > 1e-16) {
        pol2 <- lead.one(pol2, +1)
    }
    pol2[, , 1] <- diag(1, d)  # set leading polyn = unity exactly.
    return(pol2)
}

lead.one <- function(polyn, add = 0) {

    #
    # title: lead.one
    # 
    # description:  Function to add or remove leading identity matrix from
    # matrix polynomial. Is used by function 'check.one'.
    #
    # param: polyn matrix polynomial = array of ar- or ma-coefficients
    # param: add +1 or -1 for addition or removal of leading unity matrix
    #
    # return: polyn (an array) with or without a leading unity matrix
    #
    
    d <- dim(polyn)[1]
    dd <- dim(polyn)
    if (is.null(d)) {
        cat("Empty polynomial in call to lead.one,\n")
    }
    if (length(dim(polyn)) < 2) {
        cat("Argument not polynomial in call to lead.one")
        cat(dd, c(polyn), "\n")
    }
    l <- dim(polyn)[3]
    d <- dim(polyn)[1]
    if (add == 1) {
        out.pol <- array(0, dim = c(d, d, (l + 1)))
        out.pol[, , 1] <- diag(1, d)
        out.pol[, , 2:(l + 1)] <- polyn
    }
    if (add == -1) {
        out.pol <- array(0, dim = c(d, d, (l - 1)))
        out.pol[, , 1:(l - 1)] <- polyn[, , 2:l]
    }
    return(out.pol)
}

##' @title pol.order
##' 
##' @description Function to evaluate (significant)
##' order of matrix polynomium.
##'
##' @param polyn the polynomium the order of which is determined.
##' 
##' @param digits number of significant digits to be considered (values
##' smaller than 10^(-digits) are taken to be 0 (zero).
##'
##' @return pol.order order of polynomium polyn.
##' (exclusive the leading unity matrix if present.
##' pol.order=0 corresponds to the k by k unity matrix)
##'
##' @examples
##' pol        <- array(1e-8*rnorm(96),dim=c(4,4,6))
##' pol[,,1:3] <- array(rnorm(48),dim=c(4,4,3))
##' pol.order(polyn=pol,digits=12)
##' pol.order(polyn=pol,digits=4)
##'
##' @export

pol.order <- function(polyn = NULL, digits = 12) {
    pol.order <- NULL
    if (is.array(polyn)) {
        if (dim(polyn)[1] == dim(polyn)[2]) {
            polyn <- check.one(polyn)
            d <- dim(polyn)
            for (i in 1:d[3]) {
                S <- round(sum(abs(polyn[, , i])), digits = digits)
                if (S > 0) {
                  pol.order <- (i - 1)
                }
            }
        }
    }
    return(pol.order)
}

