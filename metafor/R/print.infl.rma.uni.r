print.infl.rma.uni <-
function (x, digits, ...) 
{
    if (class(x) != "infl.rma.uni") 
        stop("Argument 'x' must be an object of class \"infl.rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(digits)) 
        digits <- x$digits
    inf <- round(x$inf, digits)
    dfb <- round(x$dfb, digits)
    inf$inf <- ifelse(!is.na(x$is.infl) & x$is.infl, "*", "")
    any.na <- apply(is.na(inf), 1, any) | apply(is.na(dfb), 1, 
        any)
    if (any(any.na)) {
        if (na.act == "na.omit") {
            inf <- inf[x$not.na, ]
            dfb <- dfb[x$not.na, , drop = FALSE]
            out <- list(inf = inf, dfb = dfb)
        }
        if (na.act == "na.exclude" || na.act == "na.pass") {
            out <- list(inf = inf, dfb = dfb)
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    else {
        out <- list(inf = inf, dfb = dfb)
    }
    if (x$p == 1) {
        out <- cbind(inf[-ncol(inf)], dfb, inf[ncol(inf)])
        colnames(out)[ncol(out) - 1] <- "dfb"
    }
    print(out)
}
