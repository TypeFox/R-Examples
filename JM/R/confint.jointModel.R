confint.jointModel <-
function (object, parm = c("all", "Longitudinal", "Event"), level = 0.95, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    parm <- match.arg(parm)    
    cf <- switch(parm, 
        "Longitudinal" = fixef(object),
        "Event" = fixef(object, "Event"),
        "all" = {
            cy <- fixef(object)
            names(cy) <- paste("Y.", names(cy), sep = "")
            ct <- fixef(object, "Event")
            names(ct) <- paste("T.", names(ct), sep = "")
            c(cy, ct)
        })
    pnames <- names(cf)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- format.perc2(a, 3)
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(cf), 3L), dimnames = list(names(cf), 
            c(pct[1], "est.", pct[2])))
    ses <- sqrt(diag(vcov(object))) 
    ii <- switch(parm,
       "Longitudinal" = grep("Y.", names(ses), fixed = TRUE)[seq_along(cf)],
       "Event" = grep("T.", names(ses), fixed = TRUE)[seq_along(cf)],
       "all" = {
            iy <- grep("Y.", names(ses), fixed = TRUE)
            it <- grep("T.", names(ses), fixed = TRUE)
            c(iy[-length(iy)], it[seq(1, length(cf) - length(iy) + 1)])
       } 
    )
    ses <- ses[ii]
    ci[, c(1,3)] <- cf + ses %o% fac
    ci[, 2] <- cf
    ci
}
