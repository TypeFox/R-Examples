checkargs <- function(A, k, alpha, beta, maxit, tol, gamma, ini.W, ini.zeta, my.seed, variation, quiet, nreg, p, n){
    if(!is.matrix(A)) {stop("The input A is not a matrix. Consider as.matrix(A)?")}
    if(!is.numeric(A)) {stop("The input A is not a numeric matrix.")}
    if(!all(A >= 0)) {stop("Not all entries are non-negative.")}
    if(length(k) != 1 | !is.wholenumber(k) | k <= 0) {stop("Something is wrong about k. It should be a positive integer.")}
    if(length(alpha) != 1 | !is.numeric(alpha) | alpha < 0) {stop("Argument 'alpha' must be non-negative.")}
    if(length(beta) != 1 | !is.numeric(beta) | beta < 0) {stop("Argument 'beta' must be non-negative.")}
    if(length(maxit) != 1 | !is.wholenumber(maxit) | maxit <= 0) {stop("Something is wrong about maxit. It should be a positive integer.")}
    if(length(tol) != 1 | !is.numeric(tol) | tol < 0 | tol > 1) {stop("tol must be a number in [0,1)")}
    if(gamma != FALSE){
        if(length(gamma) != 1 | !is.numeric(gamma) | gamma < 0 | gamma > 1) {stop("Argument 'gamma' must be a number in [0,1), or 'FALSE'.")}
    }
    if(!is.null(ini.W)){
        if(!is.matrix(ini.W) | !is.numeric(ini.W) | !all(ini.W >= 0) | nrow(ini.W) != nrow(A) | ncol(ini.W) != k){stop("ini.W must be a p by k non-negative numeric matrix.")}
    }
    if(!is.null(ini.zeta)){
        if(!is.matrix(ini.zeta) | !is.logical(ini.zeta) | nrow(ini.zeta) != nrow(A) | ncol(ini.zeta) != ncol(A)) {stop("ini.zeta must be a logical matrix of the same size as A.")}
        if(sum(c(!ini.zeta)) > round(gamma * p * n)) {stop("ini.zeta contains too many FALSES (outliers); Increase trimming percentage?")}
    }
    if(!is.null(my.seed)){
        if(length(my.seed) != 1 | !is.numeric(my.seed) | !is.wholenumber(my.seed) | my.seed <= 0) {stop("Something is wrong about my.seed. It should be a positive integer.")}
    }
    if(length(variation) != 1 | !is.character(variation) | !(variation %in% c('col', 'row', 'cell', 'smooth', 'row', 'rowsmooth'))) {stop("Argument variation must be one of the following strings 'col', 'row', 'cell', 'smooth', 'row' or 'rowsmooth'")}
    if(length(quiet) != 1 | !is.logical(quiet)) {stop("The argument 'quiet' must be logical.")}
    if(length(nreg) != 1 | !is.numeric(nreg) | !is.wholenumber(nreg) | nreg <= 0) {stop("The argument 'nreg' must be a positive ineger.")}
}
