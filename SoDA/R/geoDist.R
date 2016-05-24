"geoDist" <-
function(lat1, lon1, lat2, lon2, NAOK = TRUE, DUP = TRUE) {
    n <- unique(c(length(lat1), length(lon1), length(lat2),
                  length(lon2)))
    nok <- n
    if(length(n)>1)
        stop("Need all arguments of the same length:  got ",
             paste(n, collapse=", "))
    if(n < 1)
      return(numeric())
    nas <- is.na(lat1) | is.na(lat2) | is.na(lon1) | is.na(lon2)
    if(NAOK) {
        if(any(nas)) {
            ok <- !nas
            lat1 <- lat1[ok]; lon1 <- lon1[ok];
            lat2 <- lat2[ok]; lon2 <- lon2[ok];
            nok <- sum(ok)
        }
    }
    else if(any(nas))
      stop("NA values found but not allowed")
    res <- .Fortran("GEODISTV",
                    as.double(lat1),as.double(lon1),
                    as.double(lat2),as.double(lon2),
                     dist = double(nok), as.integer(nok),
                    DUP = DUP, PACKAGE = "SoDA")$dist
    if(NAOK && any(nas)) {
        value <- rep(NA, n)
        value[ok] <- res
        value
    }
    else
        res
}

