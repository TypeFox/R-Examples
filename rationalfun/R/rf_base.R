##' Construction of rational functions
##'
##' Construction of rational functions.
##'
##' A rational function object could be constructed either by calling
##' \code{rationalfun()} or by calling \code{rationalfun.poly()}.
##'
##' \code{rationalfun()} constructs a rational function from the
##' coefficient vectors of the numerator and the denominator.
##' For example, consider a rational function
##' \eqn{R(x) = P(x) / Q(x)} where
##' \deqn{P(x) = p_1 + p_2 x + p_3 x^2 + \dots + p_k x^{k-1}}{P(x) = p[1] + p[2] * x + p[3] * x^2 + ... + p[k]* x^(k-1)} and
##' \deqn{Q(x) = q_1 + q_2 x + q_3 x^2 + \dots + q_m x^{m-1}}{Q(x) = q[1] + q[2] * x + q[3] * x^2 + ... + q[m]* x^(m-1)},
##' you may call \code{rationalfun(p[1:k], q[1:m])} to build the object.
##'
##' For \code{rationalfun.poly()}, it receives two objects of class
##' "polynomial" from the \pkg{polynom} package, representing the
##' polynomials of the numerator and the denominator respectively.
##' Use this function if you already have objects of "polynomial"
##' class, typically by calling \code{\link[polynom]{polynomial}()},
##' \code{\link[polynom]{poly.calc}()} or
##' \code{\link[polynom]{poly.orth}()}.
##'
##' \code{rfun()} and \code{rfun.poly()} are aliases of
##' \code{rationalfun()} and \code{rationalfun.poly()} in order to
##' type fewer letters.
##' 
##' The value returned by \code{rationalfun()} and
##' \code{rationalfun.poly()}
##' is an object of class "rationalfun". You can coerce the object to
##' a function, by calling \code{\link{as.function.rationalfun}()}, or to a
##' character string, by calling \code{\link{as.character.rationalfun}()}.
##'
##' Objects of "ratioanlfun" class support basic operators including
##' \code{"+"}, \code{"-"}, \code{"*"}, \code{"/"} and \code{"^"}.
##' To evaluate a rational function at a given vector, use
##' \code{\link{predict.rationalfun}()}. To compute the derivative and
##' integral in \strong{explicit} form, call
##' \code{\link{deriv.rationalfun}()} and
##' \code{\link{integral.rationalfun}()} respectively.
##' @param numer in \code{rationalfun()}, the coefficient vector of
##' the numerator; in \code{rationalfun.poly()}, an object of class
##' "polynom" in \pkg{polynom} package representing the numerator
##' @param denom similar to \code{numer}, but for the denominator
##' @return An object of class "rationalfun".
##' @seealso \code{\link[polynom]{polynomial}},
##' \code{\link[polynom]{poly.calc}}, \code{\link[polynom]{poly.orth}}
##' @importFrom stats deriv predict
##' @importFrom polynom integral
##' @export
##' @rdname rationalfun
##' @examples # (x + 1) / (x^2 + 2 * x + 3)
##' r1 <- rationalfun(c(1, 1), c(3, 2, 1))
##' print(r1)
##' # Construct from objects of "polynomial" class
##' if(require(polynom)) {
##'     p1 <- poly.calc(c(1, 2))
##'     p2 <- polynomial(rep(1, 5))
##'     r2 <- rfun.poly(p1, p2)
##'     print(r2)
##' }
# Functions in "polynom" package will handle most of the exceptions
rationalfun <- function(numer = c(0, 1), denom = c(1, 1, 1))
{
    fun <- list(numerator = polynomial(numer),
                denominator = polynomial(denom));
    if(.is_zero_polynomial(fun$denominator))
        stop("denominator should not be 0");
    return(structure(fun, class = "rationalfun"));
}

##' @rdname rationalfun
##' @export
# Alias of "rationalfun". Allow you to type fewer letters.
rfun <- rationalfun;

##' @rdname rationalfun
##' @export
rationalfun.poly <- function(numer = polynomial(c(0, 1)),
                             denom = polynomial(c(1, 1, 1)))
{
    if(!is.polynomial(numer) | !is.polynomial(denom))
        stop("numer and denom should objects of class 'polynomial'");
    fun <- list(numerator = numer,
                denominator = denom);
    if(.is_zero_polynomial(fun$denominator))
        stop("denominator should not be 0");
    return(structure(fun, class = "rationalfun"));
}

##' @rdname rationalfun
##' @export
# Alias of "rationalfun.poly"
rfun.poly <- rationalfun.poly;

##' Convert object to function
##'
##' This function converts an object of class "rationalfun" to a function.
##' @S3method as.function rationalfun
##' @method as.function rationalfun
##' @param x an object of class "rationalfun"
##' @param \dots not used in this function
##' @return A function with one argument which could be a real
##' or complex vector.
##' @export
##' @seealso \code{\link[polynom]{as.function.polynomial}}
##' @examples r <- rationalfun(c(1, 1), c(3, 2, 1))
##' r
##' f <- as.function(r)
##' f
##' f(1:10)
##' f(1:10 + 2i)
as.function.rationalfun <- function(x, ...)
{
    numer.expr <- .poly2expr(x$numerator, "numer");
    denom.expr <- .poly2expr(x$denominator, "denom");
    ex <- call("{");
    len <- 1;
    for(i in 2:length(numer.expr))
    {
        ex[[i + len - 1]] <- numer.expr[[i]];
    }
    len <- length(ex);
    for(i in 2:length(denom.expr))
    {
        ex[[i + len - 1]] <- denom.expr[[i]];
    }
    ex[[length(ex) + 1]] <- call("/", as.name("numer"), as.name("denom"));
    f <- function(x) NULL;
    body(f) <- ex;
    return(f);
}

##' Convert object to character
##'
##' This function converts an object of class "rationalfun" to a
##' character string.
##' @S3method as.character rationalfun
##' @method as.character rationalfun
##' @param x an object of class "rationalfun"
##' @param \dots not used in this function
##' @return A character string representing the rational function.
##' @export
##' @seealso \code{\link[polynom]{as.character.polynomial}}
##' @examples r <- rationalfun(c(1, 1), c(3, 2, 1))
##' as.character(r)
as.character.rationalfun <- function(x, ...)
{
    numer <- as.character(x$numerator);
    denom <- as.character(x$denominator);
    ration <- sprintf("(%s) / (%s)", numer, denom);
    return(ration);
}

##' Evaluate a rational function
##'
##' Evaluate a rational function at a real or complex vector.
##' @S3method predict rationalfun
##' @method predict rationalfun
##' @param object an object of class "rationalfun"
##' @param newdata a vector at which evaluation is requested.
##' @param \dots not used in this function
##' Both real and complex vectors are accepted.
##' @return A vector of evaluated results.
##' @export
##' @seealso \code{\link[polynom]{predict.polynomial}}
##' @examples r <- rationalfun(c(1, 1), c(3, 2, 1))
##' predict(r, 1:10)
predict.rationalfun <- function(object, newdata, ...)
{
    numer.eval <- predict(object$numerator, newdata);
    denom.eval <- predict(object$denominator, newdata);
    return(numer.eval / denom.eval);
}

##' Print a rational function
##'
##' Print a rational function in a fraction form.
##' @S3method print rationalfun
##' @method print rationalfun
##' @param x an object of class "rationalfun"
##' @param \dots not used in this function
##' @return Invisible, the object itself.
##' @export
##' @seealso \code{\link[polynom]{print.polynomial}}
##' @examples r <- rationalfun(c(1, 1), c(3, 2, 1))
##' print(r)
print.rationalfun <- function(x, ...)
{
    numer <- as.character(x$numerator);
    denom <- as.character(x$denominator);
    numer.nch <- nchar(numer);
    denom.nch <- nchar(denom);
    nspace <- floor(abs(denom.nch - numer.nch) / 2);
    if(numer.nch < denom.nch)
    {
        cat(rep(" ", nspace), numer, "\n", sep = "");
        cat(rep("-", denom.nch), "\n", denom, "\n", sep = "");
    } else {
        cat(numer, "\n", sep = "");
        cat(rep("-", numer.nch), "\n", sep = "");
        cat(rep(" ", nspace), denom, "\n", sep = "");
    }
}

