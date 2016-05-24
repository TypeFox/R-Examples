#' Difference in Deviances
#'
#' @param X a data matrix.
#' @param LF_alt Observed logistic factors.
#' @param LF_null Null logistic factors.
#' @author Wei Hao
#'
#' @importFrom lfa af
#'
#' @keywords internal
devdiff = function(X, LF_alt, LF_null = NULL){
    if(is.null(LF_null)){
        LF_null = matrix(1, ncol(X), 1)
    }

    m = nrow(X)

    F_alt  = af(X, LF_alt)
    F_null = af(X, LF_null)

    sapply(1:m, function(i){devdiff_snp(X[i,], F_alt[i,], F_null[i,])})
}

#' @keywords internal
devdiff_snp = function(snp, p1, p0){
    devalt  = sum(snp*log(p1) + (2-snp)*log(1-p1))
    devnull = sum(snp*log(p0) + (2-snp)*log(1-p0))

    -2*(devnull-devalt)
}

#' Difference in Deviances, in Parallel
#'
#' @param X a data matrix.
#' @param LF_alt Observed logistic factors.
#' @param LF_null Null logistic factors.
#' @author Wei Hao
#'
#' @importFrom lfa af
#' @importFrom lfa af_snp
#' @importFrom parallel mclapply
#'
#' @keywords internal
devdiff_parallel = function(X, LF_alt, LF_null=NULL, numcores=1){
    if(is.null(LF_null)){
        LF_null = matrix(1, ncol(X), 1)
    }

    m = nrow(X)

    devdiff_parallel_snp = function(snp, LF_alt, LF_null){
        p0 = af_snp(snp, LF_null)
        p1 = af_snp(snp, LF_alt)

        devalt  = sum(snp*log(p1) + (2-snp)*log(1-p1))
        devnull = sum(snp*log(p0) + (2-snp)*log(1-p0))

        -2*(devnull-devalt)
    }


    simplify2array(mclapply(1:m, function(i){devdiff_parallel_snp(X[i,], LF_alt, LF_null)}, mc.cores=numcores))
}

#' Mcfadden's Pseudo R-sqaured
#'
#' @param X a data matrix.
#' @param LF_alt Observed logistic factors.
#' @param LF_null Null logistic factors.
#' @author Wei Hao
#'
#' @importFrom lfa af
#'
#' @keywords internal
pseudo_Rsq = function(X, LF_alt, LF_null = NULL){
    if(is.null(LF_null)){
        LF_null = matrix(1, ncol(X), 1)
    }

    m=nrow(X)

    F_alt  = af(X, LF_alt)
    F_null = af(X, LF_null)

    sapply(1:m, function(i){mcfadden_Rsq_snp(X[i,], F_alt[i,], F_null[i,])})
}

#' @keywords internal
mcfadden_Rsq_snp = function(snp, p1, p0){
    #check for p's = 0 or 1
    IND = (p0 != 0) & (p0 != 1) & (p1 != 0) & (p1 != 1)
    p1 = p1[IND]
    p0 = p0[IND]
    snp = snp[IND]

    llalt  = sum(snp*log(p1) + (2-snp)*log(1-p1))
    llnull = sum(snp*log(p0) + (2-snp)*log(1-p0))

    1 - (llalt/llnull)
}

#' Efron's Pseudo R-sqaured
#'
#' @param X a data matrix.
#' @param LF Observed logistic factors.
#' @author Wei Hao
#'
#' @importFrom lfa af
#'
#' @keywords internal
efron_Rsq = function(X, LF){
    m=nrow(X)

    F = af(X, LF)

    sapply(1:m, function(i){efron_Rsq_snp(X[i,], F[i,])})
}

#' @keywords internal
efron_Rsq_snp = function(snp, p1){
    IND = (p1 != 0) & (p1 != 1)
    p1 = p1[IND]
    snp = snp[IND]

    y = as.numeric(c(snp > 0, snp == 2))
    p = c(p1, p1)
    ybar = mean(y)

    1 - sum( (y-p)^2 )/sum( (y-ybar)^2 )
}
