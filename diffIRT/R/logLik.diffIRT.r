logLik.diffIRT = function (object, ...) {
    LL = object$totLL/-2
    attr(LL, "df") = object$npars
    attr(LL, "nobs") = object$N
    class(LL) = "logLik"
    LL    
}