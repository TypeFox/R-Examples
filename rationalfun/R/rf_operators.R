##' Simplify a rational function
##'
##' Simplify a rational function by dropping terms whose coefficients
##' are close to zero, and then reducing it to an irreducible form.
##' @export
##' @param x an object of class "rationalfun"
##' @param \dots currently not used in this function
##' @return A new object of class "rationalfun"representing the simplified
##' rational function.
##' @examples # (x + 1) / (x^2 + 2 * x + 1) ==> 1 / (x + 1)
##' r <- rationalfun(c(1, 1), c(1, 2, 1))
##' simplify(r)
simplify <- function(x, ...)
{
    numer.coef <- coef(x$numerator);
    denom.coef <- coef(x$denominator);
    numer.coef <- round(numer.coef, 12);
    denom.coef <- round(denom.coef, 12);
    numer <- polynomial(numer.coef);
    denom <- polynomial(denom.coef);
    
    gcd <- .GCD(numer, denom);
    numer <- numer / gcd;
    denom <- denom / gcd;
    ration <- rationalfun(round(coef(numer), 12), round(coef(denom), 12));
    return(ration);
}

.add <- function(e1, e2, ...)
{
    if(missing(e2)) return(e1);
    e1 <- simplify(e1);
    e2 <- simplify(e2);
    p1 <- e1$numerator;
    q1 <- e1$denominator;
    p2 <- e2$numerator;
    q2 <- e2$denominator;
    gcd <- .GCD(q1, q2);
    l1 <- q1 / gcd;
    l2 <- q2 / gcd;
    denom <- q1 * q2 / gcd;
    numer <- p1 * l2 + p2 * l1;
    return(rationalfun.poly(numer, denom));
}

.subtract <- function(e1, e2, ...)
{
    if(missing(e2)) return(rationalfun.poly(-e1$numerator, e1$denominator));
    return(.add(e1, .subtract(e2)));
}

.multiply <- function(e1, e2, ...)
{
    e1 <- simplify(e1);
    e2 <- simplify(e2);
    p1 <- e1$numerator;
    q1 <- e1$denominator;
    p2 <- e2$numerator;
    q2 <- e2$denominator;
    e1 <- rationalfun.poly(p1, q2);
    e2 <- rationalfun.poly(p2, q1);
    e1 <- simplify(e1);
    e2 <- simplify(e2);
    numer <- e1$numerator * e2$numerator;
    denom <- e1$denominator * e2$denominator;
    return(rationalfun.poly(numer, denom));
}

.divide <- function(e1, e2, ...)
{
    e2 <- rationalfun.poly(e2$denominator, e2$numerator);
    return(.multiply(e1, e2));
}

.power <- function(e1, e2, ...)
{
    pp <- e1$numerator;
    qq <- e1$denominator;
    numer <- pp^e2;
    denom <- qq^e2;
    return(rationalfun.poly(numer, denom));
}

##' Operators for rational functions
##'
##' Basic arithmetic operators for rational functions.
##' @S3method Ops rationalfun
##' @method Ops rationalfun
##' @param e1 an object of class "rationalfun"
##' @param e2 for \code{"^"}, a positive integer; in other cases,
##' an object of class "rationalfun"
##' @return A new object of "rationalfun" class.
##' @export
##' @examples r1 <- rationalfun(c(1, 2), c(1, 2, 1))
##' r2 <- rationalfun(c(1, 1), c(1, -2, 1))
##' r1 + r2
##' r1 * r2
##' r1^2
Ops.rationalfun <- function(e1, e2)
{
	if(missing(e2))
		return(switch(.Generic,
					  "+" = e1,
					  "-" = .subtract(e1),
					  stop("unsupported unary operation")));
	e1.op.e2 <- switch(.Generic,
					   "+" = .add(e1, e2),
					   "-" = .subtract(e1, e2),
					   "*" = .multiply(e1, e2),
					   "/" = .divide(e1, e2),
					   "^" = .power(e1, e2),
					   stop("unsupported operation on rational functions"));
	return(e1.op.e2);
}
