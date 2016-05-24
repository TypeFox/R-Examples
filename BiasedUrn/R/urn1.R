# Package BiasedUrn, file urn1.R 
# R interface to univariate noncentral hypergeometric distributions

# *****************************************************************************
#    dFNCHypergeo
#    Mass function, Fisher's NonCentral Hypergeometric distribution
# *****************************************************************************
dFNCHypergeo <-
function(x, m1, m2, n, odds, precision=1E-7)  {
   stopifnot(is.numeric(x), is.numeric(m1), is.numeric(m2),
   is.numeric(n), is.numeric(odds), is.numeric(precision));
   .Call("dFNCHypergeo", 
   as.integer(x),         # Number of red balls drawn, scalar or vector
   as.integer(m1),        # Number of red balls in urn
   as.integer(m2),        # Number of white balls in urn
   as.integer(n),         # Number of balls drawn from urn
   as.double(odds),       # Odds of getting a red ball among one red and one white
   as.double(precision),  # Precision of calculation
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    dWNCHypergeo
#    Mass function, Wallenius' NonCentral Hypergeometric distribution
# *****************************************************************************
dWNCHypergeo <-
function(x, m1, m2, n, odds, precision=1E-7 ) {
   stopifnot(is.numeric(x), is.numeric(m1), is.numeric(m2),
   is.numeric(n), is.numeric(odds), is.numeric(precision));
   .Call("dWNCHypergeo", 
   as.integer(x),         # Number of red balls drawn, scalar or vector
   as.integer(m1),        # Number of red balls in urn
   as.integer(m2),        # Number of white balls in urn
   as.integer(n),         # Number of balls drawn from urn
   as.double(odds),       # Odds of getting a red ball among one red and one white
   as.double(precision),  # Precision of calculation
   PACKAGE = "BiasedUrn");
}   


# *****************************************************************************
#    pFNCHypergeo
#    Cumulative distribution function for
#    Fisher's NonCentral Hypergeometric distribution
# *****************************************************************************
pFNCHypergeo <-
function(x, m1, m2, n, odds, precision=1E-7, lower.tail=TRUE) {
   stopifnot(is.numeric(x), is.numeric(m1), is.numeric(m2), is.numeric(n),
   is.numeric(odds), is.numeric(precision), is.vector(lower.tail));
   .Call("pFNCHypergeo", 
   as.integer(x),          # Number of red balls drawn, scalar or vector
   as.integer(m1),         # Number of red balls in urn
   as.integer(m2),         # Number of white balls in urn
   as.integer(n),          # Number of balls drawn from urn
   as.double(odds),        # Odds of getting a red ball among one red and one white
   as.double(precision),   # Precision of calculation
   as.logical(lower.tail), # TRUE: P(X <= x), FALSE: P(X > x)
   PACKAGE = "BiasedUrn");
}

# *****************************************************************************
#    pWNCHypergeo
#    Cumulative distribution function for
#    Wallenius' NonCentral Hypergeometric distribution
# *****************************************************************************
pWNCHypergeo <-
function(x, m1, m2, n, odds, precision=1E-7, lower.tail=TRUE) {
   stopifnot(is.numeric(x), is.numeric(m1), is.numeric(m2), is.numeric(n),
   is.numeric(odds), is.numeric(precision), is.vector(lower.tail));
   .Call("pWNCHypergeo", 
   as.integer(x),          # Number of red balls drawn, scalar or vector
   as.integer(m1),         # Number of red balls in urn
   as.integer(m2),         # Number of white balls in urn
   as.integer(n),          # Number of balls drawn from urn
   as.double(odds),        # Odds of getting a red ball among one red and one white
   as.double(precision),   # Precision of calculation
   as.logical(lower.tail), # TRUE: P(X <= x), FALSE: P(X > x)
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    qFNCHypergeo
#    Quantile function for
#    Fisher's NonCentral Hypergeometric distribution.
#    Returns the lowest x for which P(X<=x) >= p when lower.tail = TRUE
#    Returns the lowest x for which P(X >x) <= p when lower.tail = FALSE
# *****************************************************************************
# Note: qWNCHypergeo if more accurate than qFNCHypergeo when odds = 1
qFNCHypergeo <-
function(p, m1, m2, n, odds, precision=1E-7, lower.tail=TRUE) {
   stopifnot(is.numeric(p), is.numeric(m1), is.numeric(m2), is.numeric(n),
   is.numeric(odds), is.numeric(precision), is.vector(lower.tail));
   .Call("qFNCHypergeo", 
   as.double(p),           # Cumulative probability
   as.integer(m1),         # Number of red balls in urn
   as.integer(m2),         # Number of white balls in urn
   as.integer(n),          # Number of balls drawn from urn
   as.double(odds),        # Odds of getting a red ball among one red and one white
   as.double(precision),   # Precision of calculation
   as.logical(lower.tail), # TRUE: P(X <= x), FALSE: P(X > x)
   PACKAGE = "BiasedUrn");
}   


# *****************************************************************************
#    qWNCHypergeo
#    Quantile function for
#    Wallenius' NonCentral Hypergeometric distribution.
#    Returns the lowest x for which P(X<=x) >= p when lower.tail = TRUE
#    Returns the lowest x for which P(X >x) <= p when lower.tail = FALSE
# *****************************************************************************
qWNCHypergeo <-
function(p, m1, m2, n, odds, precision=1E-7, lower.tail=TRUE) {
   stopifnot(is.numeric(p), is.numeric(m1), is.numeric(m2), is.numeric(n),
   is.numeric(odds), is.numeric(precision), is.vector(lower.tail));
   .Call("qWNCHypergeo", 
   as.double(p),           # Cumulative probability
   as.integer(m1),         # Number of red balls in urn
   as.integer(m2),         # Number of white balls in urn
   as.integer(n),          # Number of balls drawn from urn
   as.double(odds),        # Odds of getting a red ball among one red and one white
   as.double(precision),   # Precision of calculation
   as.logical(lower.tail), # TRUE: P(X <= x), FALSE: P(X > x)
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    rFNCHypergeo
#    Random variate generation function for
#    Fisher's NonCentral Hypergeometric distribution.
# *****************************************************************************
rFNCHypergeo <-
function(nran, m1, m2, n, odds, precision=1E-7) {
   stopifnot(is.numeric(nran), is.numeric(m1), is.numeric(m2),
   is.numeric(n), is.numeric(odds), is.numeric(precision));
   .Call("rFNCHypergeo", 
   as.integer(nran),       # Number of random variates desired
   as.integer(m1),         # Number of red balls in urn
   as.integer(m2),         # Number of white balls in urn
   as.integer(n),          # Number of balls drawn from urn
   as.double(odds),        # Odds of getting a red ball among one red and one white
   as.double(precision),   # Precision of calculation
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    rWNCHypergeo
#    Random variate generation function for
#    Wallenius' NonCentral Hypergeometric distribution.
# *****************************************************************************
rWNCHypergeo <-
function(nran, m1, m2, n, odds, precision=1E-7) {
   stopifnot(is.numeric(nran), is.numeric(m1), is.numeric(m2),
   is.numeric(n), is.numeric(odds), is.numeric(precision));
   .Call("rWNCHypergeo", 
   as.integer(nran),       # Number of random variates desired
   as.integer(m1),         # Number of red balls in urn
   as.integer(m2),         # Number of white balls in urn
   as.integer(n),          # Number of balls drawn from urn
   as.double(odds),        # Odds of getting a red ball among one red and one white
   as.double(precision),   # Precision of calculation
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    meanFNCHypergeo
#    Calculates the mean of
#    Fisher's NonCentral Hypergeometric distribution.
# *****************************************************************************
meanFNCHypergeo <- function(
   m1,                  # Number of red balls in urn
   m2,                  # Number of white balls in urn
   n,                   # Number of balls drawn from urn
   odds,                # Odds of getting a red ball among one red and one white
   precision=1E-7) {    # Precision of calculation
   stopifnot(is.numeric(m1), is.numeric(m2), is.numeric(n), 
   is.numeric(odds), is.numeric(precision));
   .Call("momentsFNCHypergeo", as.integer(m1), as.integer(m2),         
   as.integer(n), as.double(odds), as.double(precision),
   as.integer(1),       # 1 for mean, 2 for variance
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    meanWNCHypergeo
#    Calculates the mean of
#    Wallenius' NonCentral Hypergeometric distribution.
# *****************************************************************************
meanWNCHypergeo <- function(
   m1,                  # Number of red balls in urn
   m2,                  # Number of white balls in urn
   n,                   # Number of balls drawn from urn
   odds,                # Odds of getting a red ball among one red and one white
   precision=1E-7) {    # Precision of calculation
   stopifnot(is.numeric(m1), is.numeric(m2), is.numeric(n), 
   is.numeric(odds), is.numeric(precision));
   .Call("momentsWNCHypergeo", as.integer(m1), as.integer(m2),         
   as.integer(n), as.double(odds), as.double(precision),
   as.integer(1),       # 1 for mean, 2 for variance
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    varFNCHypergeo
#    Calculates the variance of
#    Fisher's NonCentral Hypergeometric distribution.
# *****************************************************************************
varFNCHypergeo <- function(
   m1,                  # Number of red balls in urn
   m2,                  # Number of white balls in urn
   n,                   # Number of balls drawn from urn
   odds,                # Odds of getting a red ball among one red and one white
   precision=1E-7) {    # Precision of calculation
   stopifnot(is.numeric(m1), is.numeric(m2), is.numeric(n), 
   is.numeric(odds), is.numeric(precision));
   .Call("momentsFNCHypergeo", as.integer(m1), as.integer(m2),         
   as.integer(n), as.double(odds), as.double(precision),
   as.integer(2),       # 1 for mean, 2 for variance
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    varWNCHypergeo
#    Calculates the variance of
#    Wallenius' NonCentral Hypergeometric distribution.
# *****************************************************************************
varWNCHypergeo <- function(
   m1,                  # Number of red balls in urn
   m2,                  # Number of white balls in urn
   n,                   # Number of balls drawn from urn
   odds,                # Odds of getting a red ball among one red and one white
   precision=1E-7) {    # Precision of calculation
   stopifnot(is.numeric(m1), is.numeric(m2), is.numeric(n), 
   is.numeric(odds), is.numeric(precision));
   .Call("momentsWNCHypergeo", as.integer(m1), as.integer(m2),         
   as.integer(n), as.double(odds), as.double(precision),
   as.integer(2),       # 1 for mean, 2 for variance
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    modeFNCHypergeo
#    Calculates the mode of
#    Fisher's NonCentral Hypergeometric distribution.
# *****************************************************************************
# Note: The result is exact regardless of the precision parameter.
# The precision parameter is included only for analogy with modeWNCHypergeo.
modeFNCHypergeo <- function(
   m1,                  # Number of red balls in urn
   m2,                  # Number of white balls in urn
   n,                   # Number of balls drawn from urn
   odds,                # Odds of getting a red ball among one red and one white
   precision=0) {       # Precision of calculation
   stopifnot(is.numeric(m1), is.numeric(m2), is.numeric(n), 
   is.numeric(odds));
   .Call("modeFNCHypergeo", as.integer(m1), as.integer(m2),         
   as.integer(n), as.double(odds), 
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    modeWNCHypergeo
#    Calculates the mode of
#    Fisher's NonCentral Hypergeometric distribution.
# *****************************************************************************
modeWNCHypergeo <- function(
   m1,                  # Number of red balls in urn
   m2,                  # Number of white balls in urn
   n,                   # Number of balls drawn from urn
   odds,                # Odds of getting a red ball among one red and one white
   precision=1E-7) {    # Precision of calculation
   stopifnot(is.numeric(m1), is.numeric(m2), is.numeric(n), 
   is.numeric(odds), is.numeric(precision));
   .Call("modeWNCHypergeo", as.integer(m1), as.integer(m2),         
   as.integer(n), as.double(odds), as.double(precision),
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    oddsFNCHypergeo
#    Estimate odds ratio from mean for
#    Fisher's NonCentral Hypergeometric distribution
# *****************************************************************************
# Uses Cornfield's approximation. Specified precision is ignored.
oddsFNCHypergeo <-
function(mu, m1, m2, n, precision=0.1)  {
   stopifnot(is.numeric(mu), is.numeric(m1), is.numeric(m2),
   is.numeric(n), is.numeric(precision));
   .Call("oddsFNCHypergeo", 
   as.double(mu),         # Observed mean of x1
   as.integer(m1),        # Number of red balls in urn
   as.integer(m2),        # Number of white balls in urn
   as.integer(n),         # Number of balls drawn from urn
   as.double(precision),  # Precision of calculation
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    oddsWNCHypergeo
#    Estimate odds ratio from mean for
#    Wallenius' NonCentral Hypergeometric distribution
# *****************************************************************************
oddsWNCHypergeo <-
function(mu, m1, m2, n, precision=0.1)  {
   stopifnot(is.numeric(mu), is.numeric(m1), is.numeric(m2),
   is.numeric(n), is.numeric(precision));
   .Call("oddsWNCHypergeo", 
   as.double(mu),         # Observed mean of x1
   as.integer(m1),        # Number of red balls in urn
   as.integer(m2),        # Number of white balls in urn
   as.integer(n),         # Number of balls drawn from urn
   as.double(precision),  # Precision of calculation
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    numFNCHypergeo
#    Estimate number of balls of each color from experimental mean for
#    Fisher's NonCentral Hypergeometric distribution
# *****************************************************************************
# Uses Cornfield's approximation. Specified precision is ignored.
numFNCHypergeo <-
function(mu, n, N, odds, precision=0.1)  {
   stopifnot(is.numeric(mu), is.numeric(n), is.numeric(N),
   is.numeric(odds), is.numeric(precision));
   .Call("numFNCHypergeo", 
   as.double(mu),         # Observed mean of x1
   as.integer(n),         # Number of balls sampled
   as.integer(N),         # Number of balls in urn before sampling
   as.double(odds),       # Odds of getting a red ball among one red and one white
   as.double(precision),  # Precision of calculation (ignored)
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    numWNCHypergeo
#    Estimate number of balls of each color from experimental mean for
#    Wallenius' NonCentral Hypergeometric distribution
# *****************************************************************************
# Uses approximation. Specified precision is ignored.
numWNCHypergeo <-
function(mu, n, N, odds, precision=0.1)  {
   stopifnot(is.numeric(mu), is.numeric(n), is.numeric(N),
   is.numeric(odds), is.numeric(precision));
   .Call("numWNCHypergeo", 
   as.double(mu),         # Observed mean of x1
   as.integer(n),         # Number of balls sampled
   as.integer(N),         # Number of balls in urn before sampling
   as.double(odds),       # Odds of getting a red ball among one red and one white
   as.double(precision),  # Precision of calculation (ignored)
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    minHypergeo
#    Minimum of x for central and noncentral Hypergeometric distributions
# *****************************************************************************
minHypergeo <- function(m1, m2, n) {
   stopifnot(m1>=0, m2>=0, n>=0, n<=m1+m2);
   max(n-m2, 0);
}


# *****************************************************************************
#    maxHypergeo
#    Maximum of x for central and noncentral Hypergeometric distributions
# *****************************************************************************
maxHypergeo <- function(m1, m2, n) {
   stopifnot(m1>=0, m2>=0, n>=0, n<=m1+m2);
   min(m1, n);
}   
