#### Abstract out the .C() calls into these auxiliaries;
#### Makes it easier later to replace by calls to  (*, PACKAGE = "splines")

###--> ../tests/spline-ex.R  for investigating these
###    ~~~~~~~~~~~~~~~~~~~~

.splBasis <- function(ord, knots, ncoef, xo, derivs = rep(0, n))
{
    ## Purpose: encapsulate .C("spline_basis", ..)
    ## ----------------------------------------------------------------------
    ## Arguments: from result of B-spline fit, see ?cobs , etc
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 20 Feb 2002, 14:46

    ord <- as.integer(ord)
    new.knots <- c(rep.int(knots[1], ord-1:1),
                   knots,
                   rep.int(knots[length(knots)], ord-1:1))

    if(ord + length(knots) != ncoef + 2)
        warning(".splBasis(): (ord,length(knots),ncoef)=",
                paste(ord,length(knots),ncoef, sep=", "),
                " -- not ``matching'' ?\n")

    n <- length(xo <- as.double(xo))

    .C(spline_basis,
       as.double(new.knots),
       as.integer(ncoef),
       ord, # "order"
       xo, # "xvals"
       derivs = as.integer(derivs),
       n,
       design = array(0, c(ord, n)),# "basis"
       offsets = integer(n))[c("design","offsets")]
}

.splValue <- function(degree, knots, coef, xo)
{
    ## Purpose: encapsulate .C("spline_value", ..)
    ## ----------------------------------------------------------------------
    ## Arguments: from result of B-spline fit, see ?cobs , etc
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 20 Feb 2002, 13:48

    degree <- as.integer(degree)
    ord <- as.integer(degree + 1)
    new.knots <- c(rep.int(knots[1], degree),
                   knots,
                   rep.int(knots[length(knots)], degree))
    derivs <- as.integer(0)
    n <- length(xo)
    .C(spline_value,
       as.double(new.knots),
       as.double(coef),
       length(coef),
       ord,
       as.double(xo),
       as.integer(n),
       derivs,
       y = double(n))$y
} ## .splValue
