# Package BiasedUrn, file urn2.R 
# R interface to multivariate noncentral hypergeometric distributions

# *****************************************************************************
#    dMFNCHypergeo
#    Mass function for
#    Multivariate Fisher's NonCentral Hypergeometric distribution
# *****************************************************************************
dMFNCHypergeo <-
function(
   x,                   # Number of balls drawn of each color, vector or matrix
   m,                   # Number of balls of each color in urn, vector
   n,                   # Number of balls drawn from urn, scalar
   odds,                # Odds for each color, vector
   precision=1E-7) {    # Precision of calculation, scalar
   stopifnot(is.numeric(x), is.numeric(m), is.numeric(n), is.numeric(odds), is.numeric(precision));
   
   # Convert x to integer vector or matrix without loosing dimensions:
   if (is.matrix(x)) {   
      xx <- matrix(as.integer(x), nrow=dim(x)[1], ncol=dim(x)[2]);
   }
   else {
      xx <- as.integer(x);
   }
   .Call("dMFNCHypergeo", xx, as.integer(m), as.integer(n),         
   as.double(odds), as.double(precision), PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    dMWNCHypergeo
#    Mass function for
#    Multivariate Wallenius' NonCentral Hypergeometric distribution
# *****************************************************************************
dMWNCHypergeo <-
function(
   x,                   # Number of balls drawn of each color, vector or matrix
   m,                   # Number of balls of each color in urn, vector
   n,                   # Number of balls drawn from urn, scalar
   odds,                # Odds for each color, vector
   precision=1E-7) {    # Precision of calculation, scalar
   stopifnot(is.numeric(x), is.numeric(m), is.numeric(n), is.numeric(odds), is.numeric(precision));
   
   # Convert x to integer vector or matrix without loosing dimensions:
   if (is.matrix(x)) {   
      xx <- matrix(as.integer(x), nrow=dim(x)[1], ncol=dim(x)[2]);
   }
   else {
      xx <- as.integer(x);
   }
   .Call("dMWNCHypergeo", xx, as.integer(m), as.integer(n),         
   as.double(odds), as.double(precision), PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    rMFNCHypergeo
#    Random variate generation function for
#    Multivariate Fisher's NonCentral Hypergeometric distribution.
# *****************************************************************************
rMFNCHypergeo <-
function(nran, m, n, odds, precision=1E-7) {
   stopifnot(is.numeric(nran), is.numeric(m),
   is.numeric(n), is.numeric(odds), is.numeric(precision));
   .Call("rMFNCHypergeo", 
   as.integer(nran),       # Number of random variates desired, scalar
   as.integer(m),          # Number of balls of each color in urn, vector
   as.integer(n),          # Number of balls drawn from urn, scalar
   as.double(odds),        # Odds for each color, vector
   as.double(precision),   # Precision of calculation, scalar
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    rMWNCHypergeo
#    Random variate generation function for
#    Multivariate Wallenius' NonCentral Hypergeometric distribution.
# *****************************************************************************
rMWNCHypergeo <-
function(nran, m, n, odds, precision=1E-7) {
   stopifnot(is.numeric(nran), is.numeric(m),
   is.numeric(n), is.numeric(odds), is.numeric(precision));
   .Call("rMWNCHypergeo", 
   as.integer(nran),       # Number of random variates desired, scalar
   as.integer(m),          # Number of balls of each color in urn, vector
   as.integer(n),          # Number of balls drawn from urn, scalar
   as.double(odds),        # Odds for each color, vector
   as.double(precision),   # Precision of calculation, scalar
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    momentsMFNCHypergeo
#    Calculates the mean and variance of the
#    Multivariate Fisher's NonCentral Hypergeometric distribution.
#    Results are returned as a data frame.
# *****************************************************************************
momentsMFNCHypergeo <- function(
   m,                   # Number of balls of each color in urn, vector
   n,                   # Number of balls drawn from urn, scalar
   odds,                # Odds for each color, vector
   precision = 0.1) {   # Precision of calculation, scalar
   stopifnot(is.numeric(m), is.numeric(n), 
   is.numeric(odds), is.numeric(precision));
   res <- .Call("momentsMFNCHypergeo", as.integer(m), 
   as.integer(n), as.double(odds), as.double(precision),
   PACKAGE = "BiasedUrn");
   # Convert result to data frame
   colnames(res) <- list("xMean","xVariance")
   as.data.frame(res);   
}


# *****************************************************************************
#    momentsMWNCHypergeo
#    Calculates the mean and variance of the
#    Multivariate Wallenius' NonCentral Hypergeometric distribution.
#    Results are returned as a data frame.
# *****************************************************************************
momentsMWNCHypergeo <- function(
   m,                   # Number of balls of each color in urn, vector
   n,                   # Number of balls drawn from urn, scalar
   odds,                # Odds for each color, vector
   precision = 0.1) {   # Precision of calculation, scalar
   stopifnot(is.numeric(m), is.numeric(n), 
   is.numeric(odds), is.numeric(precision));
   res <- .Call("momentsMWNCHypergeo", as.integer(m), 
   as.integer(n), as.double(odds), as.double(precision),
   PACKAGE = "BiasedUrn");
   # Convert result to data frame
   colnames(res) <- list("xMean","xVariance")
   as.data.frame(res);   
}


# *****************************************************************************
#    meanMFNCHypergeo
#    Calculates the mean of the
#    Multivariate Fisher's NonCentral Hypergeometric distribution.
# *****************************************************************************
meanMFNCHypergeo <- function(
   m,                   # Number of balls of each color in urn, vector
   n,                   # Number of balls drawn from urn, scalar
   odds,                # Odds for each color, vector
   precision = 0.1) {   # Precision of calculation, scalar
   momentsMFNCHypergeo(m, n, odds, precision)$xMean
}


# *****************************************************************************
#    meanMWNCHypergeo
#    Calculates the mean of the
#    Multivariate Wallenius' NonCentral Hypergeometric distribution.
# *****************************************************************************
meanMWNCHypergeo <- function(
   m,                   # Number of balls of each color in urn, vector
   n,                   # Number of balls drawn from urn, scalar
   odds,                # Odds for each color, vector
   precision = 0.1) {   # Precision of calculation, scalar
   momentsMWNCHypergeo(m, n, odds, precision)$xMean
}


# *****************************************************************************
#    varMFNCHypergeo
#    Calculates the variance of the
#    Multivariate Fisher's NonCentral Hypergeometric distribution.
# *****************************************************************************
varMFNCHypergeo <- function(
   m,                   # Number of balls of each color in urn, vector
   n,                   # Number of balls drawn from urn, scalar
   odds,                # Odds for each color, vector
   precision = 0.1) {   # Precision of calculation, scalar
   momentsMFNCHypergeo(m, n, odds, precision)$xVariance
}


# *****************************************************************************
#    varMWNCHypergeo
#    Calculates the variance of the
#    Multivariate Wallenius' NonCentral Hypergeometric distribution.
# *****************************************************************************
varMWNCHypergeo <- function(
   m,                   # Number of balls of each color in urn, vector
   n,                   # Number of balls drawn from urn, scalar
   odds,                # Odds for each color, vector
   precision = 0.1) {   # Precision of calculation, scalar
   momentsMWNCHypergeo(m, n, odds, precision)$xVariance
}


# *****************************************************************************
#    oddsMFNCHypergeo
#    Estimate odds ratio from mean for the
#    Multivariate Fisher's NonCentral Hypergeometric distribution
# *****************************************************************************
# Uses Cornfield's approximation. Specified precision is ignored.
oddsMFNCHypergeo <-
function(mu, m, n, precision=0.1)  {
   stopifnot(is.numeric(mu), is.numeric(m), is.numeric(n), is.numeric(precision));
   # Convert mu to double vector or matrix without loosing dimensions:
   if (is.matrix(mu)) {   
      mux <- matrix(as.double(mu), nrow=dim(mu)[1], ncol=dim(mu)[2]);
   }
   else {
      mux <- as.double(mu);
   }
   .Call("oddsMFNCHypergeo", 
   mux,                   # Observed mean of each x, vector
   as.integer(m),         # Number of balls of each color in urn, vector
   as.integer(n),         # Number of balls drawn from urn, scalar
   as.double(precision),  # Precision of calculation, scalar
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    oddsMWNCHypergeo
#    Estimate odds ratio from mean for the
#    Multivariate Wallenius' NonCentral Hypergeometric distribution
# *****************************************************************************
# Uses approximation. Specified precision is ignored.
oddsMWNCHypergeo <-
function(mu, m, n, precision=0.1)  {
   stopifnot(is.numeric(mu), is.numeric(m), is.numeric(n), is.numeric(precision));
   # Convert mu to double vector or matrix without loosing dimensions:
   if (is.matrix(mu)) {   
      mux <- matrix(as.double(mu), nrow=dim(mu)[1], ncol=dim(mu)[2]);
   }
   else {
      mux <- as.double(mu);
   }
   .Call("oddsMWNCHypergeo", 
   mux,                   # Observed mean of each x, vector
   as.integer(m),         # Number of balls of each color in urn, vector
   as.integer(n),         # Number of balls drawn from urn, scalar
   as.double(precision),  # Precision of calculation, scalar
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    numMFNCHypergeo
#    Estimate number of balls of each color from experimental mean for
#    Multivariate Fisher's NonCentral Hypergeometric distribution
# *****************************************************************************
# Uses Cornfield's approximation. Specified precision is ignored.
numMFNCHypergeo <-
function(mu, n, N, odds, precision=0.1)  {
   stopifnot(is.numeric(mu), is.numeric(n), is.numeric(N), is.numeric(odds), is.numeric(precision));
   # Convert mu to double vector or matrix without loosing dimensions:
   if (is.matrix(mu)) {   
      mux <- matrix(as.double(mu), nrow=dim(mu)[1], ncol=dim(mu)[2]);
   }
   else {
      mux <- as.double(mu);
   }
   .Call("numMFNCHypergeo", 
   mux,                   # Observed mean of each x, vector
   as.integer(n),         # Number of balls drawn from urn, scalar
   as.integer(N),         # Number of balls in urn before sampling, scalar
   as.double(odds),       # Odds for each color, vector
   as.double(precision),  # Precision of calculation, scalar (ignored)
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    numMWNCHypergeo
#    Estimate number of balls of each color from experimental mean for
#    Multivariate Wallenius' NonCentral Hypergeometric distribution
# *****************************************************************************
# Uses approximation. Specified precision is ignored.
numMWNCHypergeo <-
function(mu, n, N, odds, precision=0.1)  {
   stopifnot(is.numeric(mu), is.numeric(n), is.numeric(N), is.numeric(odds), is.numeric(precision));
   # Convert mu to double vector or matrix without loosing dimensions:
   if (is.matrix(mu)) {   
      mux <- matrix(as.double(mu), nrow=dim(mu)[1], ncol=dim(mu)[2]);
   }
   else {
      mux <- as.double(mu);
   }
   .Call("numMWNCHypergeo", 
   mux,                   # Observed mean of each x, vector
   as.integer(n),         # Number of balls drawn from urn, scalar
   as.integer(N),         # Number of balls in urn before sampling, scalar
   as.double(odds),       # Odds for each color, vector
   as.double(precision),  # Precision of calculation, scalar (ignored)
   PACKAGE = "BiasedUrn");
}


# *****************************************************************************
#    minMHypergeo
#    Minimum of x for central and noncentral 
#    Multivariate Hypergeometric distributions
# *****************************************************************************
#    m = Number of balls of each color in urn, vector
#    n = Number of balls drawn from urn, scalar
minMHypergeo <- function(m, n)  {
   stopifnot(m>=0, n>=0, n<=sum(m));
   pmax(n - sum(m) + m, 0);
}


# *****************************************************************************
#    maxMHypergeo
#    Maximum of x for central and noncentral 
#    Multivariate Hypergeometric distributions
# *****************************************************************************
#    m = Number of balls of each color in urn, vector
#    n = Number of balls drawn from urn, scalar
maxMHypergeo <- function(m, n)  {
   stopifnot(m>=0, n>=0, n<=sum(m));
   pmin(m, n);
}   
