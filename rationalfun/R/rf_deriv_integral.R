##' Differentiate a rational function
##'
##' Calculate the derivative of a rational function. The returned value
##' result is still an object of class "rationalfun".
##' @S3method deriv rationalfun
##' @method deriv rationalfun
##' @export
##' @param expr an object of class "rationalfun"
##' @param \dots not used in this function
##' @return An object of class "rationalfun" representing
##' the derivative of the original rational function.
##' @seealso \code{\link[polynom]{deriv.polynomial}},
##' \code{\link[stats]{deriv}}
##' @examples # (x + 1) / (x^2 + x + 1)
##' r <- rationalfun(c(1, 1), c(1, 1, 1))
##' deriv(r)
deriv.rationalfun <- function(expr, ...)
{
    expr <- simplify(expr);
    pp <- expr$numerator;
    qq <- expr$denominator;
    numer <- deriv(pp) * qq - pp * deriv(qq);
    denom <- qq * qq;
    return(rationalfun.poly(numer, denom));
}

integraterf.transcendent <- function(x, ...)
{
    numer <- x$numerator;
    denom <- x$denominator;
    numercoef <- coef(numer);
    if(numercoef %*% numercoef < 1e-20) return(substitute(0 * x));
    
    ddenom <- deriv(denom);
    coefrf <- rationalfun(coef(numer), coef(ddenom));
    
    fac <- solve(denom);
    fac.re <- Re(fac);
    fac.im <- Im(fac);
    real.index <- abs(fac.im) < 1e-10;
    real.roots <- fac.re[real.index];
    fac.re <- fac.re[!real.index];
    fac.im <- fac.im[!real.index];
    conju.index <- duplicated(fac.re) & duplicated(abs(fac.im));
    conju.roots.re <- fac.re[conju.index];
    conju.roots.im <- fac.im[conju.index];
    conju.roots.inner <- conju.roots.re^2 + conju.roots.im^2;
    
    real.coefs <- predict(coefrf, real.roots);
    conju.coefs <- predict(coefrf, complex(real = conju.roots.re,
                                           imaginary = conju.roots.im));
    conju.coefs.re <- Re(conju.coefs);
    conju.coefs.im <- Im(conju.coefs);
    
    expr <- NULL;
    if(length(real.coefs))
    {
        term0 <- substitute(log(abs(ex)),
                           list(ex = .linear(real.roots[1])));
        term <- call("*", real.coefs[1], term0);
        expr <- term;
        for(i in seq_along(real.coefs)[-1])
        {
            term0 <- substitute(log(abs(ex)),
                                list(ex = .linear(real.roots[i])));
            term <- call("*", abs(real.coefs[i]), term0);
            op <- if(real.coefs[i] > 0) "+" else "-";
            expr <- call(op, expr, term);
        }
    }
    if(length(conju.coefs))
    {
        term01 <- substitute(log(ex),
                             list(ex = .quadratic(-2 * conju.roots.re[1],
                                                  conju.roots.inner[1])));
        if(is.null(expr))
        {
            expr <- call("*", conju.coefs.re[1], term01);
        } else {
            term1 <- call("*", abs(conju.coefs.re[1]), term01);
            op <- if(conju.coefs.re[1] > 0) "+" else "-";
            expr <- call(op, expr, term1);
        }
        term02 <- substitute(atan(ex),
                             list(ex = .frac(conju.roots.re[1],
                                             conju.roots.im[1])));
        term2 <- call("*", abs(2 * conju.coefs.im[1]), term02);
        op <- if(conju.coefs.im[1] > 0) "-" else "+";
        expr <- call(op, expr, term2);
        for(i in seq_along(conju.coefs)[-1])
        {
            term01 <- substitute(log(ex),
                                 list(ex = .quadratic(-2 * conju.roots.re[i],
                                                 conju.roots.inner[i])));
            term1 <- call("*", abs(conju.coefs.re[i]), term01);
            op <- if(conju.coefs.re[i] > 0) "+" else "-";
            expr <- call(op, expr, term1);
            term02 <- substitute(atan(ex),
                                 list(ex = .frac(conju.roots.re[i],
                                                conju.roots.im[i])));
            term2 <- call("*", abs(2 * conju.coefs.im[i]), term02);
            op <- if(conju.coefs.im[i] > 0) "-" else "+";
            expr <- call(op, expr, term2);
        }
    }
    return(expr);
}

##' Integrate a rational function
##'
##' Calculate the integral of a rational function. See "Details".
##'
##' The returned value is a function call with argument named "x".
##' That is, the integral is an expression in R with an explicit form,
##' which could be evaluated directly by calling \code{\link{eval}()},
##' or indirectly using the \code{\link{int2fun}()} function.
##'
##' The algorithm is based on the Hermite-Ostrogradski formula which is
##' discussed in the reference. See the article for more details.
##' @S3method integral rationalfun
##' @method integral rationalfun
##' @export
##' @param expr an object of class "rationalfun"
##' @param \dots not used in this function
##' @return A function call representing the explicit form of the
##' integral.
##' @seealso \code{\link[polynom]{integral.polynomial}}
##' @references T. N. Subramaniam, and Donald E. G. Malm,
##' How to Integrate Rational Functions,
##' \emph{The American Mathematical Monthly},
##' Vol. 99, No.8 (1992), 762-772.
##' @examples # (x + 1) / (x^2 + x + 1)
##' r <- rationalfun(c(1, 1), c(1, 1, 1))
##' expr <- integral(r)
##' # Evaluate the call directly
##' eval(expr, list(x = 2))
##' # Use int2fun()
##' f <- int2fun(expr)
##' f(2)
integral.rationalfun <- function(expr, ...)
{
    x <- expr;
    pp <- x$numerator;
    qq <- x$denominator;
    poly.part <- pp / qq;
    poly.expr <- parse(text = as.character(poly.part))[[1]];
    pp <- pp %% qq;
    qq1 <- .GCD(qq, deriv(qq));
    if(.degree(qq1) == 0)
    {
        x <- rationalfun.poly(pp, qq);
        x <- simplify(x);
        transcendent.expr <- integraterf.transcendent(x);
        expr <- call("+", transcendent.expr, poly.expr);
        return(expr);
    }
    qq2 <- qq / qq1;
    ss <- deriv(qq1) * qq2 / qq1;
    # P = P1' * Q2 - P1 * S + P2 * Q1
    np <- .degree(qq2);
    nq <- .degree(qq1);
    xx <- seq(0.295782039048272, 5.277857726515673, length.out = np + nq);
    p1p <- outer(xx, 0:(nq - 1), function(x, k) k * x^(k - 1));
    term1 <- sweep(p1p, 1, predict(qq2, xx), "*");
    p1 <- outer(xx, 0:(nq - 1), function(x, k) x^k);
    term2 <- sweep(p1, 1, predict(ss, xx), "*");
    p2 <- outer(xx, 0:(np - 1), function(x, k) x^k);
    term3 <- sweep(p2, 1, predict(qq1, xx), "*");
    A <- cbind(term1 - term2, term3);
    coefs <- solve(A, predict(pp, xx));
    p1.coef <- coefs[1:nq];
    p2.coef <- coefs[-(1:nq)];
    
    rational.part <- rationalfun(p1.coef, coef(qq1));
    rational.part <- simplify(rational.part);
    rational.expr <- parse(text = as.character(rational.part))[[1]];
    transcendent.part <- rationalfun(p2.coef, coef(qq2));
    transcendent.part <- simplify(transcendent.part);
    transcendent.expr <- integraterf.transcendent(transcendent.part);
    expr <- call("+", transcendent.expr, rational.expr);
    expr <- call("+", expr, poly.expr);
    return(expr);
}

##' Convert a call to a function
##'
##' Convert a function call to a function in R. In this package, the
##' function is typically used to convert the result of
##' \code{\link{integral.rationalfun}()} to a function with one
##' argument.
##' @export
##' @param expr a function call, typically returned by
##' \code{\link{integral.rationalfun}()}.
##' @return A function with one argument which could be a real
##' or complex vector.
##' @seealso \code{\link[polynom]{integral.polynomial}}
##' @examples x <- rationalfun(c(-6, -1, -8, 15, -1, 8, -9, 2),
##'                            c(8, 12, 16, 4, 4))
##' int <- integral(x)
##' fun <- int2fun(int)
##' fun(c(0, 1))
int2fun <- function(expr)
{
    f <- function(x) NULL;
    body(f) <- expr;
    return(f);
}

