polyreg <-
function (x, y, w, degree = 1, monomial = FALSE, ...) 
{
    x <- polybasis(x, degree, monomial)
    y <- as.matrix(y)                   # just making sure ...
    if (iswt <- !missing(w)) {
        if (any(w <= 0)) 
            stop("only positive weights")
        w <- sqrt(w)
        y <- y * w
        x <- x * w
    }
    qrx <- qr(x)
    coef <- as.matrix(qr.coef(qrx, y))
    fitted <- qr.fitted(qrx, y)
    if ((df <- qrx$rank) < ncol(x)) 
        coef[qrx$pivot, ] <- coef
    if (iswt) 
        fitted <- fitted/w
    structure(list(fitted.values = fitted, coefficients = coef, 
        degree = degree, monomial = monomial, df = df), class = "polyreg")
}

