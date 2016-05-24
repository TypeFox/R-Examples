# Code from "polynom" package with some modifications
.poly2expr <- function(x, var.name)
{
    a <- rev(coef(x));
    w <- as.name(var.name);
    v <- as.name("x");
    ex <- call("{", call("<-", w, 0));
    for(i in seq_along(a))
    {
        ex[[i + 2]] <- call("<-", w, call("+", a[1], call("*", v, w)));
        a <- a[-1];
    }
    return(ex);
}

.is_zero_polynomial <- function(x)
{
    cf <- coef(x);
	return(cf %*% cf < 1e-16);
}

.degree <- function(x) length(unclass(x)) - 1;

.GCD <- function(x, y)
{
    if(.is_zero_polynomial(y)) x
    else if(.degree(y) == 0) as.polynomial(1)
    else Recall(y, x %% y)
}

.LCM <- function(x, y)
{
    if(.is_zero_polynomial(x) || .is_zero_polynomial(y))
        return(as.polynomial(0))
    (x / .GCD(x, y)) * y
}
#####################################################

# Functions used to generate calls
# Expression of x - a
.linear <- function(a)
{
    expr = if(a > 0) substitute(x - a, list(a = a))
           else if(a < 0) substitute(x + a, list(a = -a))
           else substitute(x)
    return(expr);
}
# Expression of x^2 + b * x + c
.quadratic <- function(bb, cc)
{
    expr <- substitute(x^2);
    if(bb != 0)
    {
        op <- if(bb > 0) "+" else "-";
        expr <- call(op, expr, substitute(b * x, list(b = abs(bb))));
    }
    if(cc != 0)
    {
        op <- if(cc > 0) "+" else "-";
        expr <- call(op, expr, abs(cc));
    }
    return(expr);
}
# Expression of (x - a) / b
.frac <- function(a, b)
{
    expr <- call("/", .linear(a), abs(b));
    expr <- if(b >= 0) expr else call("-", expr);
    return(expr);
}

