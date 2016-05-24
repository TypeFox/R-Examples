#' @encoding UTF-8
#' @title Stukel's test of the logistic link
#' @description The Stukel's test is an alternative to the goodness-of-fit test for logistic regression.
#'  It tests if significant change occurs in the model with the addition of new coefficients.
#'
#' @param object An object of class \code{glm}.
#'
#' @param alternative add both \code{z1} and \code{z2} to model or just one of them.
#'
#' @details Two new covariates, z1 and z2 are generated such that \deqn{z1 = 0.5 \* logit^{2} * I(pi >= 0.5)}, \deqn{z2 = - 0.5 \* logit^{2} \* I(pi <= 0.5)}, where \deqn{I(arg) = 1} if arg is \code{TRUE} and \deqn{I(arg) = 1} if \code{FALSE}.
#'
#' @note Adapted from program published by Brett Presnell's code available at the Florida University.
#'
#' @keywords Tests
#'
#'@references
#' Stukel, T.A. (1988) Generalized logistic models. \emph{Journal of the American Statistical Association} 83: 426-431.
#'@references
#' Hosmer, David W., et al (1997) A comparison of goodness-of-fit tests for the logistic regression model. \emph{Statistics in medicine} 16.9, 965-980.
#' @references
#'Allison, Paul (2014) \emph{Another Goodness-of-Fit Test for Logistic Regression}.
#' @export
#'
#' @importFrom stats predict vcov residuals pchisq
#'
`stukel` <- function(object, alternative = c("both", "alpha1", "alpha2")) {
  DNAME <- deparse(substitute(object))
  METHOD <- "Stukel's test of the logistic link"
  alternative <- match.arg(alternative)
  eta <- predict(object, type = "link")
  etasq <- 0.5 * eta * eta
  etapos <- eta > 0
  dv <- matrix(0, nrow = length(eta), ncol = 2)
  dv[etapos,1] <- etasq[etapos]
  dv[!etapos,2] <- - etasq[!etapos]
  colnames(dv) <- c("z1","z2")
  oinfo <- vcov(object)
  oX <- qr.X(object$qr)
  ImH <- - oX %*% oinfo %*% t(oX)
  diag(ImH) <- 1 + diag(ImH)
  wdv <- sqrt(object$weights) * dv
  qmat <- t(wdv) %*% ImH %*% wdv
  sc <- apply(dv * (object$weights * residuals(object, "working")), 2, sum)
  allstat <- c(sc * sc / diag(qmat), sc %*% solve(qmat) %*% sc)
  names(allstat) <- c("alpha1", "alpha2", "both")
  allpar <- c(1,1,2)
  names(allpar) <- names(allstat)
  allpval <- pchisq(allstat, allpar, lower.tail=FALSE)
  STATISTIC <- allstat[alternative]
  PARAMETER <- allpar[alternative]
  names(PARAMETER) <- "df"
  PVAL <- allpval[alternative]
  names(allpar) <- rep("df", 3)
  structure(list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = PVAL,
                 alternative = alternative,
                 method = METHOD, data.name = DNAME,
                 allstat = allstat, allpar = allpar, allpval = allpval
  ),
  class = "htest")
}### end -- stukel function
NULL






#' @encoding UTF-8
#' @title Variance Inflation Factor
#' @description Extracts Variance Inflation Factor from a model of class \dQuote{lm}
#' @param model a model object
#' @param \dots further arguments passed to or used by other methods.
#' @return A numeric value indicating the variance inflation in the model
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'#' @keywords Tests
#'
#' @importFrom stats coef vcov terms cov2cor model.matrix coefficients
#'
#' @keywords Models, Tests
#' @examples
#' data(mtcars)
#'
#' m1 <- lm(mpg ~ qsec + hp, data=mtcars)
#'
#' vif(m1)
#'
#' @export
`vif` <-
  function(model, ...) {
    if (any(is.na(coef(model))))
      stop ("there are aliased coefficients in the model")
    v <- vcov(model)
    assign <- attributes(model.matrix(model))$assign
    if (names(coefficients(model)[1]) == "(Intercept)") {
      v <- v[-1, -1]
      assign <- assign[-1]
    }
    else warning("Vifs may not be sensible without intercept.")
    parameters <- labels(terms(model))
    n.parameters <- length(parameters)
    if (n.parameters < 2) stop("model contains fewer than 2 parameters")
    R <- cov2cor(v)
    detR <- det(R)
    result <- matrix(0, n.parameters, 3)
    rownames(result) <- parameters
    colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
    for (parameters in 1:n.parameters) {
      subs <- which(assign == parameters)
      result[parameters, 1] <- det(as.matrix(R[subs, subs])) *
        det(as.matrix(R[-subs, -subs])) / detR
      result[parameters, 2] <- length(subs)
    }
    if (all(result[, 2] == 1)) result <- result[, 1]
    else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
    result
  }### end -- vif function
NULL





#' @encoding UTF-8
#' @title Geary's test for normality
#' @description This function computes an estimator of Geary's measure of kurtosis.
#' @param x the numeric vector.
#' @param na.rm A logical for NA values.
#' @details Null hypothesis is that the data obeys to normal distribution and that data should have kurtosis equal to 3.
#' @return statistic The Geary's test of statistic G.
#' @return p.value The significant probability of the null-hypothesis testing.
#' @author Daniel Marcelino \email{dmarcelino@@live.com}
#' @keywords Tests
#' @examples
#' set.seed(1234)
#' x = rnorm(1000)
#' geary(x)
#'
#' geary(20:50)
#'
#' y = c(0.269, 0.357, 0.2, 0.221, 0.275, 0.277, 0.253, 0.127, 0.246)
#'
#' stats::qqnorm(y)
#' @export
`geary` <- function(x, na.rm=TRUE) {
  if (any(i.na <- is.na(x))) {
    if(na.rm)
      x <- x[!i.na]
    else return(NA)
  }
  DNAME <- deparse(substitute(x))
  mu <- mean(x)
  n <- length(x)
  kurt <- n*sum( (x-mean(x))^4 )/(sum( (x-mean(x))^2 )^2);
  G <- sum(abs(x-mu))/sqrt(n*sum((x-mu)^2))
  pval <- (1-stats::pnorm((G-sqrt(2/pi))/sqrt(1-3/pi)*sqrt(n)))*2
  RVAL <- list(statistic = c(kurt = kurt, G = G), p.value = pval,
               method = "Geary's test for normality",
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
} ### end -- geary function
NULL




#' @encoding UTF-8
#' @title James-Stein shrunken estimates
#'
#' @description Computes James-Stein shrunken estimates of cell means given
#' a response variable (which may be binary) and a grouping indicator.
#' @references
#' Efron, Bradley and Morris, Carl (1977) ``Stein's Paradox in Statistics.'' \emph{Scientific American} Vol. 236 (5): 119-127.
#'
#' James, Willard and Stein, Charles (1961) ``Estimation with
#' Quadratic Loss.'' \emph{Proceedings of the Fourth Berkeley
#' Symposium on Mathematical Statistics and Probability}, Vol. 1: 361-379.
#'
#' @param y The response variable.
#' @param k The grouping factor.
#' @keywords Tests
#' @export
james.stein <- function(y, k)
{
  s <- !(is.na(y)|is.na(k))
  y <- y[s];
  k <- as.character(k[s])
  ## as.char -> unused levels OK
  k <- length(unique(k))
  if(k<3)
    stop("must have >=3 groups")

  .stats <- function(w) {
    bar <- mean(w)
    ss  <- sum((w-bar)^2)
    n <- length(w)
    ##if(n<2)
    ##  stop("a group has n<2")

    c(n=length(w), mean=bar, ss=ss, var=ss/n/(n-1))
  }

  Z <- .stats(y)
  st <- tapply(y, k, FUN=.stats)
  nams <- names(st)
  z <- matrix(unlist(st),ncol=4,byrow=TRUE)
  ssb <- .stats(z[,2])["ss"]
  shrink <- 1 - (k-3)*z[,4]/ssb
  shrink[z[,1]==1] <- 0
  shrink <- pmin(pmax(shrink,0),1)
  list(n=z[,1], mean=z[,2],
       shrunk.mean=structure(Z["mean"]*(1-shrink)+shrink*z[,2], names=nams),
       shrink=shrink)
}### end -- james.stein function
NULL




#' @encoding UTF-8
#' @title Jarque-Bera test for normality
#'
#' @description This function performs the Jarque-Bera test on the given data sample to determine if the data are sample drawn from a normal population.
#' @param x  A numeric vector of data.
#' @details The Jarque-Bera statistic is chi-square distributed with two degrees of freedom. Under the hypothesis of normality, data should be symmetrical (i.e. skewness should be equal to zero) and have skewness chose to three.
#' @keywords Tests
#' @references
#' Jarque, C. M., Bera, A. K. (1980) Efficient test for normality, homoscedasticity and serial independence of residuals,
#' Economic Letters, Vol. 6 Issue 3, 255-259.
#'
#' @examples
#' set.seed(1234)
#'
#' x <- rnorm(1000)
#'
#' jarque.bera(x)
#' @export
jarque.bera <- function(x)
{
  if ( !is.vector( x ) )
    stop( "argument x is not a vector" )
  if ( !is.numeric( x ) )
    stop( "argument x is not numeric" )
  DNAME <- deparse( substitute( x ) )
  n <- length(x)
  ALTERNATIVE <- "greater"
  METHOD <- "Jarque-Bera Normality Test"
  K <- kurtosis( x )
  S <- skewness( x )
  JB  <- ( n / 6 ) * ( S^2 + 0.25 * ( ( K - 3 )^2 ) )
  pval <- 1 - pchisq( JB, df=2 )
  JBVAL <- list( statistic=c(JB=JB), p.value=pval, alternative=ALTERNATIVE, method=METHOD,
                 data.name=DNAME )
  class( JBVAL ) <- "htest"
  return( JBVAL )
} ### end -- jarque.bera function
NULL




#' @encoding UTF-8
#' @title Jensen-Shannon Distance
#'
#' @description The Jensen-Shannon divergence or distance matrix stores the \eqn{n*(n-1)/2} pairwise distances/similarities between observations in an \eqn{n x p} matrix where n correspond to the independent observational units and p represent the covariates measured on each individual.
#' @param mat An n x p matrix.
#' @keywords Tests
#' @examples
#'
#' # create a matrix
#' n  = 10
#' m = matrix(runif(n*10), ncol = 10)
#' m = m/rowSums(m)
#'
#' jensen.shannon(m)
#'
#' @export
jensen.shannon <- function(mat) {
  kld = function(p,q) sum(ifelse(p == 0 | q == 0, 0, log(p/q)*p))
  res = matrix(0, nrow(mat), nrow(mat))
  for (i in 1:(nrow(mat) - 1)) {
    for (j in (i+1):nrow(mat)) {
      m = (mat[i,] + mat[j,])/2
      d1 = kld(mat[i,], m)
      d2 = kld(mat[j,], m)
      res[j,i] = sqrt(.5*(d1 + d2))
    }
  }
  res
}### end -- jensen.shannon function
NULL




#' @encoding UTF-8
#' @title Johnson-Neyman Regression
#' @description Probing Regression Interactions
#' @param y the dependent variable.
#' @param x the independent variable.
#' @param z the moderator variable.
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' @keywords Tests
#'@export
`johnson.neyman` = function(y,x,z){
  mod = stats::lm(y~x + z + x:z)
  gam1.v = stats::vcov(mod)["x","x"]
  gam3.v = stats::vcov(mod)["x:z","x:z"]
  gam1gam3.cv = stats::vcov(mod)["x","z"]
  df = length(y)-length(stats::coef(mod))
  t = stats::qt(.975,df)
  zz = seq(min(z),max(z),by=.01)
  se.w1 = sqrt(gam1.v + 2*zz*gam1gam3.cv + zz^2*gam3.v)
  w1.hat = stats::coef(mod)["x"]+stats::coef(mod)["x:z"]*zz
  z.tab = cbind(zz,t<abs(w1.hat/se.w1))
  ci.low = w1.hat - t*se.w1
  ci.upp = w1.hat + t*se.w1
  w1.tab = data.frame(w1.hat,z=z.tab[,1],z.tab[,2],ci.low,ci.upp)
  colnames(w1.tab) = c("Est","Z","Significant","LB95", "UB95")
  w1.tab[,3] = ifelse(w1.tab[,3]=="1","Yes","No")
  w1.tab
}### end -- johnson.neyman function
NULL





#' @encoding UTF-8
#' @title Lilliefors (Kolmogorov-Smirnov) test for normality
#'
#' @description Performs the Lilliefors (Kolmogorov-Smirnov) test for the composite hypothesis of normality. The Lilliefors (Kolomorov-Smirnov) test is the most famous EDF omnibus test for normality; compared to the Anderson-Darling test and the Cramer-von Mises test it is known to perform worse.
#'
#' @param x A numeric vector of data values, the number of observations must be greater than 4.
#'  Missing values are allowed.
#' @keywords Tests
#' @references
#' Thode Jr., H.C. (2002): Testing for Normality. Marcel Dekker, New York.
#'
#' @importFrom stats complete.cases pnorm
#' @export
#' @examples
#'
#' set.seed(1234)
#' x = rnorm(1000)
#' lilliefors(x)
#'
`lilliefors` <- function(x){
  DNAME <- deparse(substitute(x))
  x <- sort(x[stats::complete.cases(x)])
  n <- length(x)
  if (n < 5)
    stop("sample size must be greater than 4")
  p <- stats::pnorm((x - mean(x))/sd(x))
  Dplus <- max(seq(1:n)/n - p)
  Dminus <- max(p - (seq(1:n) - 1)/n)
  K <- max(Dplus, Dminus)
  if (n <= 100) {
    Kd <- K
    nd <- n
  }
  else {
    Kd <- K * ((n/100)^0.49)
    nd <- 100
  }
  pvalue <- exp(-7.01256 * Kd^2 * (nd + 2.78019) + 2.99587 *
                  Kd * sqrt(nd + 2.78019) - 0.122119 + 0.974598/sqrt(nd) +
                  1.67997/nd)
  if (pvalue > 0.1) {
    KK <- (sqrt(n) - 0.01 + 0.85/sqrt(n)) * K
    if (KK <= 0.302) {
      pvalue <- 1
    }
    else if (KK <= 0.5) {
      pvalue <- 2.76773 - 19.828315 * KK + 80.709644 *
        KK^2 - 138.55152 * KK^3 + 81.218052 * KK^4
    }
    else if (KK <= 0.9) {
      pvalue <- -4.901232 + 40.662806 * KK - 97.490286 *
        KK^2 + 94.029866 * KK^3 - 32.355711 * KK^4
    }
    else if (KK <= 1.31) {
      pvalue <- 6.198765 - 19.558097 * KK + 23.186922 *
        KK^2 - 12.234627 * KK^3 + 2.423045 * KK^4
    }
    else {
      pvalue <- 0
    }
  }
  RVAL <- list(statistic = c(D = K), p.value = pvalue, method = "Lilliefors (Kolmogorov-Smirnov) normality test",
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
} ### end -- lilliefors function
NULL





#' @encoding UTF-8
#' @title Bonett-Seier test of Geary's kurtosis
#'
#' @description Performs the Bonett-Seier test of Geary's measure of kurtosis for normally distributed data.
#' @param x A numeric vector of data values.
#' @param alternative A character string specifying the alternative hypothesis,
#'  must be one of '"two.sided"' (default), '"greater"' or '"less"'.
#'   You can specify just the initial letter
#'
#'   @details  Under the hypothesis of normality, data should have Geary's
#'    kurtosis equal to \code{sqrt(2/pi)} (0.7979). This test has such null
#'     hypothesis and is useful to detect a significant difference of Geary's
#'      kurtosis in normally distributed data.
#' @keywords Tests
#'  @references
#'  Bonett, D.G., Seier, E. (2002) A test of normality with high uniform power. Computational Statistics and Data Analysis, 40, 435-445.
#' @importFrom stats complete.cases  pnorm
#' @export
#' @examples
#' set.seed(1234)
#' x = rnorm(1000)
#' geary(x)
#' bonett.seier(x)
`bonett.seier` <-
  function (x, alternative=c("two.sided","less","greater"))
  {
    DNAME <- deparse(substitute(x))
    x <- sort(x[stats::complete.cases(x)])
    n <- length(x)
    s <- match.arg(alternative)
    alter <- switch(s, two.sided=0, less=1, greater=2)
    rho <- sqrt(sum((x-mean(x))^2)/n);
    tau <- sum(abs(x-mean(x)))/n;
    omega <- 13.29*(log(rho)-log(tau));
    z <- sqrt(n+2)*(omega-3)/3.54;
    pval <- stats::pnorm(z, lower.tail = FALSE)
    if (alter == 0) {
      pval <- 2*pval
      if (pval > 1) pval<-2-pval
      alt <- "kurtosis is not equal to sqrt(2/pi)"
    }
    else if (alter == 1)
    {
      alt <- "kurtosis is greater than sqrt(2/pi)"
    }
    else
    {
      pval <- 1-pval
      alt <- "kurtosis is lower than sqrt(2/pi)"
    }
    RVAL <- list(statistic = c(tau = tau, z = z), alternative = alt,
                 p.value = pval, method = "Bonett-Seier test for Geary kurtosis",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  } ### end -- bonett.seier function
NULL




#' @encoding UTF-8
#' @title Anscombe-Glynn test of kurtosis
#'
#' @description Performs the Anscombe-Glynn test of kurtosis for normal samples.
#'
#' @param x A numeric vector of data values.
#' @param alternative A character string specifying the alternative hypothesis,
#'  must be one of '"two.sided"' (default), '"greater"' or '"less"'. You can
#'   specify just the initial letter.
#' @details Under the hypothesis of normality, data should have kurtosis equal
#'  to 3.This test has such null hypothesis and is useful to detect a
#'  significant difference of kurtosis in normally distributed data.
#' @keywords Tests
#' @references Anscombe, F.J., Glynn, W.J. (1983) Distribution of kurtosis
#' statistic for normal statistics. Biometrika, 70, 1, 227-234
#' @importFrom stats complete.cases  pnorm
#' @export
#' @examples
#' set.seed(1234)
#' x = rnorm(1000)
#' kurtosis(x)
#' anscombe.glynn(x)
#'
`anscombe.glynn` <-
  function (x, alternative=c("two.sided","less","greater"))
  {
    DNAME <- deparse(substitute(x))
    x <- sort(x[stats::complete.cases(x)])
    n <- length(x)
    s <- match.arg(alternative)
    alter <- switch(s, two.sided=0, less=1, greater=2)
    b <- n*sum( (x-mean(x))^4 )/(sum( (x-mean(x))^2 )^2);
    eb2 <- 3*(n-1)/(n+1);
    vb2 <- 24*n*(n-2)*(n-3)/ ((n+1)^2*(n+3)*(n+5));
    m3 <- (6*(n^2-5*n+2)/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/(n*(n-2)*(n-3)));
    a <- 6 + (8/m3) * (2/m3 + sqrt(1 + 4/m3^2));
    xx <- (b-eb2)/sqrt(vb2);
    z <- ( 1-2/(9*a)-( (1-2/a) / (1+xx*sqrt(2/(a-4))) )^(1/3))/ sqrt(2/(9*a));
    pval <- stats::pnorm(z, lower.tail = FALSE)
    if (alter == 0) {
      pval <- 2*pval
      if (pval > 1) pval<-2-pval
      alt <- "kurtosis is not equal to 3"
    }
    else if (alter == 1)
    {
      alt <- "kurtosis is greater than 3"
    }
    else
    {
      pval <- 1-pval
      alt <- "kurtosis is lower than 3"
    }
    RVAL <- list(statistic = c(kurt = b, z = z), p.value = pval,
                 alternative = alt, method = "Anscombe-Glynn kurtosis test",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  }### end -- anscombe.glynn function
NULL




#' @encoding UTF-8
#' @title Anderson-Darling test for normality
#'
#' @description Performs the Anderson-Darling test for the composite hypothesis of normality. See details.
#'
#' @param x A numeric vector of data values, the number of observations must be greater than 7. Missing values are allowed.
#'  @details The Anderson-Darling test is the recommended EDF test by Stephens (1986) followed by the Cramer-von Mises test. Compared to the later, the Anderson-Darling gives more weight to the tails of the distribution.
#' @keywords Tests
#' @importFrom stats complete.cases  pnorm
#' @export
#' @examples
#' set.seed(1234)
#' x = rnorm(1000)
#' anderson.darling(x)
#'
anderson.darling <- function (x)
{
  DNAME <- deparse(substitute(x))
  x <- sort(x[stats::complete.cases(x)])
  n <- length(x)
  if (n < 8)
    stop("sample size must be greater than 7")
  logp1 <- stats::pnorm((x - mean(x))/sd(x), log.p = TRUE)
  logp2 <- stats::pnorm(-(x - mean(x))/sd(x), log.p = TRUE)
  h <- (2 * seq(1:n) - 1) * (logp1 + rev(logp2))
  A <- -n - mean(h)
  AA <- (1 + 0.75/n + 2.25/n^2) * A
  if (AA < 0.2) {
    pval <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
  }
  else if (AA < 0.34) {
    pval <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
  }
  else if (AA < 0.6) {
    pval <- exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
  }
  else if (AA < 10) {
    pval <- exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
  }
  else pval <- 3.7e-24
  RVAL <- list(statistic = c(A = A), p.value = pval, method = "Anderson-Darling normality test",
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
} ### end -- anderson.darling function
NULL






#' @title D'Agostino test of skewness
#' @description Performs the D'Agostino test for skewness in normally distributed data.
#' @param x A numeric vector of data values.
#' @param alternative A character string specifying the alternative hypothesis,
#' must be one of '"two.sided"' (default), '"greater"' or '"less"'.
#'  You can specify just the initial letter.
#'  @details Under the hypothesis of normality, data should be symmetrical (i.e. skewness should be equal to zero).
#'  This test has such null hypothesis and is useful to detect a significant skewness in normally distributed data.
#'
#' @keywords Tests
#' @references
#' D'Agostino, R.B. (1970). Transformation to Normality of the Null Distribution of G1. Biometrika, 57, 3, 679-681.
#' @importFrom stats complete.cases  pnorm
#' @export
#' @examples
#'
#' set.seed(1234)
#' x = rnorm(1000)
#'
#' skewness(x) # is data normal?
#'
#' agostino(x)
#'
`agostino` <-
  function (x, alternative=c("two.sided","less","greater"))
  {
    DNAME <- deparse(substitute(x))
    x <- sort(x[stats::complete.cases(x)])
    n <- length(x)
    s <- match.arg(alternative)
    alter <- switch(s, two.sided=0, less=1, greater=2)
    if ((n < 8 || n > 46340))
      stop("sample size must be between 8 and 46340")
    s3 <- (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
    y <- s3*sqrt((n+1)*(n+3)/(6*(n-2)))
    b2 <- 3*(n*n+27*n-70)*(n+1)*(n+3)/((n-2)*(n+5)*(n+7)*(n+9))
    w <- sqrt(-1+sqrt(2*(b2-1)));
    d <- 1/sqrt(log(w));
    a <- sqrt(2/(w*w-1));
    z <- d*log(y/a+sqrt((y/a)^2+1));
    pval <- stats::pnorm(z, lower.tail = FALSE)
    if (alter == 0) {
      pval <- 2*pval
      if (pval > 1) pval<-2-pval
      alt <- "data have a skewness"
    }
    else if (alter == 1)
    {
      alt <- "data have positive skewness"
    }
    else
    {
      pval <- 1-pval
      alt <- "data have negative skewness"
    }
    RVAL <- list(statistic = c(skew = s3, z = z), p.value = pval,
                 alternative = alt, method = "D'Agostino skewness test",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
  } ### end -- agostino function
NULL





#' @title Likelihood Ratio Test (G test) for Tables
#'
#' @description Computes the likelihood ratio test (G test) for contingency tables. Currently does not do Williams' and Yates' correction.
#'
#' @param x A vector or a matrix.
#' @param y A vector that is ignored if x is a matrix and required if x is a vector.
#' @param \dots Extra parameters pass to the \code{\link{table}} function.
#' @keywords Tests
#'
#' @references
#' Agresti, Alan (1996) \emph{Introduction to categorical data
#'  analysis}. NY: John Wiley and Sons
#'
#' Smithson, M.J. (2003) \emph{Confidence Intervals, Quantitative
#'  Applications in the Social Sciences Series, No. 140}.
#'  Thousand Oaks, CA: Sage. pp. 39-41.
#'
#' @examples
#' # 2000 General Social Survey-Sex and Party affiliation
#' # Agresti (1996) p. 38:
#'
#' gss <- data.frame(
#'    expand.grid(sex=c("female", "male"),
#'    party=c("dem", "indep", "rep")),
#'    count=c(762,484,327,239, 468,477))
#'
#' # expand it:
#' # GSS <- gss[rep(1:nrow(gss), gss[["count"]]),]
#' GSS = untable(gss, freq = "count")
#'
#' calc.LR(GSS$party, GSS$sex)
#'
#' @rdname calc.LR
#' @export
`calc.LR` <- function(x, y = NULL, ...) UseMethod("calc.LR")

#' @rdname calc.LR
#' @export
`calc.LR.default` <- function(x, y = NULL, ...) {
  DNAME <- deparse(substitute(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) { if (min(dim(x)) == 1) x <- as.vector(x)	}
  if (!is.matrix(x)) {
    if (is.null(y)) { stop("y must be a non-null vector") } else {
      if (length(x) != length(y)) { stop("x and y must have the same length") }
    }
    YDNAME <- deparse(substitute(y))
    ok <- complete.cases(x,y)
    x <- factor(x[ok])
    y <- factor(y[ok])
    if ((nlevels(x) < 2L) || (nlevels(y) < 2L)) { stop("'x' and 'y' must have at least 2 levels") }
    x <- table(x,y, ...)
    names(dimnames(x)) <- c(DNAME, YDNAME)
  }
  if (all(dim(x) == 2)) {
    result <- stats::chisq.test(x, correct = TRUE)
  } else {
      result <- suppressWarnings(stats::chisq.test(x, correct = FALSE))
      }

  G <- 2 * sum(result$obs * log(result$obs/result$exp), na.rm = TRUE)
  pvalue <- 1 - pchisq(G, result[[2]][[1]])
  return(list(statistics = G, df = result[[2]][[1]], p.value = pvalue))
} ### end -- LR function
NULL





#' @title Cramer's V Coefficient for Tables
#'
#' @description Computes the Cramer's V coefficient of association
#'  for tables. The Cramer's V is a measure of effect size for a chi-square goodness of fit test.
#'
#' @param x A vector or a matrix.
#' @param y A vector that is ignored if x is a matrix and required if x is a vector.
#' @param \dots Extra parameters pass to the \code{\link{table}} function.
#'
#' @note
#' # Bootstrap confidence intervals for Cramer's V
#' \url{http://support.sas.com/documentation/cdl/en/statugfreq/63124/PDF/default/statugfreq.pdf}, p. 1821
#' @keywords Association, Nominal, Ordinal Tests
#'
#' @references
#' Agresti, Alan (1996) \emph{Introduction to categorical data
#' analysis}. NY: John Wiley and Sons.
#'
#' @examples
#' # Consider an experiment with two conditions, each with 100 participants.
#' # Each participant chooses between one of following three parties.
#'
#' cond1 <- c(40, 25, 35)
#' cond2 <- c(25, 35, 45)
#' mat <- cbind(cond1, cond2)
#' rownames(mat) <- c( 'party1', 'party2', 'party3')
#'
#' # To test the null hypothesis that the distribution of preferences
#' # is identical in the two conditions, we run a chi-square test:
#' stats::chisq.test(mat) # still significant
#'
#' # However, if we want to estimate the effect size, we then use Cramer's V:
#' calc.CV(mat)
#'
#' # Agresti (2002), table 3.10, p. 104
#' # 1991 General Social Survey: The effect size of race on party identification.
#' gss <- data.frame(
#'    expand.grid(race=c("black", "white"),
#'    party=c("dem", "indep", "rep")),
#'    count=c(103,341,15,105,11,405))
#'
#' GSS = untable(gss, freq = "count")
#'
#' calc.CV(GSS$race, GSS$party)
#'
#' @export
#' @rdname calc.CV
`calc.CV` <- function(x, y = NULL, ...) UseMethod("calc.CV")

#' @rdname calc.CV
#' @export
`calc.CV.default` <- function(x, y = NULL, ...) {
  DNAME <- deparse(substitute(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) { if (min(dim(x)) == 1) x <- as.vector(x)	}
  if (!is.matrix(x)) {
    if (is.null(y)) { stop("y must be a non-null vector") } else {
      if (length(x) != length(y)) { stop("x and y must have the same length") }
    }
    YDNAME <- deparse(substitute(y))
    ok <- complete.cases(x,y)
    x <- factor(x[ok])
    y <- factor(y[ok])
    if ((nlevels(x) < 2L) || (nlevels(y) < 2L)) { stop("'x' and 'y' must have at least 2 levels") }
    x <- table(x,y, ...)
    names(dimnames(x)) <- c(DNAME, YDNAME)
  }
  if (all(dim(x) == 2)) {
    result <- suppressWarnings(stats::chisq.test(x, correct = FALSE))
    cramersV <- (prod(diag(result$obs)) - (result$obs[2,1]*result$obs[1,2]))/sqrt(prod(result$obs))
  } else {
    result <- suppressWarnings(stats::chisq.test(x, correct = FALSE))
    cramersV <- sqrt(result[[1]][[1]]/(sum(x)*(min(dim(x))-1)))
  }
  return(cramersV)
} ### end -- Cramers'V function
NULL





#' @title The Phi Coefficient for \code{2 x 2} Tables
#'
#' @description Computes the Phi coefficient for \code{2 x 2} tables.
#'
#' @param x A vector or a matrix.
#' @param y A vector that is ignored if x is a matrix and required if x is a vector.
#' @param \dots Extra parameters pass to the \code{\link{table}} function.
#'
#' @details
#' Phi is seldom applied for indexing a \code{2 x 2} table, because the
#' researcher will typically want to contrast the two proportions as
#' an increment or ratio, not with a correlation coefficient.
#' Alternatives to Phi are the Pearson's C; Tschuprow's T, and Cramer's V.
#'
#' @keywords Association, Nominal, Ordinal, Tests
#' @references
#' Friendly, Michael (2000) \emph{Visualizing Categorical Data}. SAS Institute Inc., p. 63.
#'
#'
#' @examples
#' # Admission to Berkeley graduate programs:
#' Berkeley <- data.frame(
#' expand.grid(GENDER=c("Male", "Female"),
#' ADMIT=c("Admitted", "Rejected")),
#' Freq=c(1198,557,1493,1278))
#'
#' tab = as.table(rbind(c(1198,557), c(1493,1278)))
#' calc.Phi(tab)
#'
#' @export
#' @rdname calc.Phi
`calc.Phi` <- function(x, y = NULL, ...) UseMethod("calc.Phi")


#' @rdname calc.Phi
#' @export
`calc.Phi.default` <- function(x, y = NULL, ...) {
  DNAME <- deparse(substitute(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) { if (min(dim(x)) == 1) x <- as.vector(x)	}
  if (!is.matrix(x)) {
    if (is.null(y)) { stop("y must be a non-null vector") } else {
      if (length(x) != length(y)) { stop("x and y must have the same length") }
    }
    YDNAME <- deparse(substitute(y))
    ok <- complete.cases(x,y)
    x <- factor(x[ok])
    y <- factor(y[ok])
    if ((nlevels(x) < 2L) || (nlevels(y) < 2L)) { stop("'x' and 'y' must have at least 2 levels") }
    x <- table(x,y, ...)
    names(dimnames(x)) <- c(DNAME, YDNAME)
  }
  if (all(dim(x) == 2)) {
    result <- as.numeric(sqrt(suppressWarnings(stats::chisq.test(x, correct = FALSE)$statistic) / sum(x)))
    #phicoef <- (prod(diag(result$observed)) - (result$observed[2,1]*result$observed[1,2]))/sqrt(prod(result$observed))
  } else {
    result <- as.numeric(sqrt(suppressWarnings(stats::chisq.test(x, correct = FALSE)$statistic) / sum(x)))
   # phicoef <- sqrt(result[[1]][[1]]/sum(x))
  }
  return(result)
} ### end -- phi function
NULL



#' @title Pearson's Contingency Coefficient for Tables
#'
#' @description Compute Pearson's contingency coefficient for an RxC contingency table. \deqn{C = \sqrt{\frac{T}{N + T}}}, where T = the chi-square test statistic and N = the total sample size.
#'
#' @param x A vector or a matrix.
#' @param y A vector that is ignored if x is a matrix and required if x is a vector.
#' @param \dots Extra parameters pass to the \code{\link{table}} function.
#'
#' @details If we have N observations with two variables where each observation can be classified into one of R mutually exclusive categories for variable one and one of C mutually exclusive categories for variable two, then a cross-tabulation of the data results in a two-way contingency table (also referred to as an RxC contingency table).
#' A common question with regards to a two-way contingency table is whether we have independence. By independence, we mean that the row and column variables are unassociated (i.e., knowing the value of the row variable will not help us predict the value of column variable and likewise knowing the value of the column variable will not help us predict the value of the row variable).
#' One criticism of this statistic is that it does not give a meaningful description of the degree of dependence (or strength of association). That is, it is useful for determining whether there is dependence. However, since the strength of that association also depends on the degrees of freedom as well as the value of the test statistic, it is not easy to interpert the strength of association.
#'
#' @note The Cramer contingency coefficient is more commonly used than the Pearson contingency coefficient.
#' @keywords Tests
#'
#' @references
#' Agresti, Alan (1996) \emph{Introduction to categorical data
#'  analysis}. NY: John Wiley and Sons.
#' Blaikie, N. 2003. \emph{Analyzing Quantative Data}. London: SAGE.
#'
#' @examples
#' # some data:
#' male <- c(33, 76, 6)
#' female <- c(47, 153, 25)
#' mat <- cbind( male, female )
#' rownames(mat) <- c( 'good', 'satisfactory', 'bad')
#'
#' calc.CC(mat)
#'
#' @export
`calc.CC` <- function(x, y = NULL, ...) UseMethod("calc.CC")

#' @rdname calc.CC
#' @export
`calc.CC.default` <- function(x, y = NULL, ...) {
  DNAME <- deparse(substitute(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) { if (min(dim(x)) == 1) x <- as.vector(x)	}
  if (!is.matrix(x)) {
    if (is.null(y)) { stop("y must be a non-null vector") } else {
      if (length(x) != length(y)) { stop("x and y must have the same length") }
    }
    YDNAME <- deparse(substitute(y))
    ok <- complete.cases(x,y)
    x <- factor(x[ok])
    y <- factor(y[ok])
    if ((nlevels(x) < 2L) || (nlevels(y) < 2L)) { stop("'x' and 'y' must have at least 2 levels") }
    x <- table(x,y, ...)
    names(dimnames(x)) <- c(DNAME, YDNAME)
  }
  if (all(dim(x) == 2)) {
    result <- suppressWarnings(stats::chisq.test(x, correct = TRUE))
    cont.coef <- (prod(diag(result$obs)) - (result$obs[2,1]*result$obs[1,2]))/sqrt(prod(result$obs))
  } else {
    result <- suppressWarnings(stats::chisq.test(x, correct = FALSE))
    cont.coef <- sqrt(result[[1]][[1]]/(sum(x) + result[[1]][[1]]))
  }
  return(cont.coef)
} ### end -- Contingency function
NULL




#' @title Tschuprow's T for Tables
#'
#' @description Computes the Tschuprow's T coefficient of association for tables.
#'
#' @param x A vector or a matrix.
#' @param y A vector that is ignored if x is a matrix and required if x is a vector.
#' @param \dots Extra parameters pass to the \code{\link{table}} function.
#'
#' @details Tschuprow's T has the disadvantage of producing an overcorrection. Although kept from being > 1, the correlation coefficient often cannot reach the permissible maximum value of 1. This problem is likely to occur if \code{R} is much greater than \code{C} (or the other way around) in a large \code{R x C} table.
#' @keywords Tests
#'
#' @references
#' Tschuprow, A. A. (1939) \emph{Principles of the Mathematical Theory of Correlation}. Translated by M. Kantorowitsch. W. Hodge & Co.
#'
#' @examples
#' # some data:
#' male <- c(33, 76, 6);
#' female <- c(47, 153, 25);
#' mat <- cbind( male, female );
#' rownames(mat) <- c( 'good', 'satisfactory', 'bad');
#'
#' calc.TT(mat);
#'
#' # long format
#' long = untable(mat);
#'
#' calc.TT(long$Var1, long$Var2)
#' @export
#' @rdname calc.TT
`calc.TT` <- function(x, y = NULL) UseMethod("calc.TT")


#' @rdname calc.TT
#' @export
`calc.TT.default` <- function(x, y = NULL, ...){

  if(!is.null(y)) x <- table(x, y, ...)
  # http://en.wikipedia.org/wiki/Tschuprow's_T
  # Hartung S. 451
  as.numeric( sqrt(suppressWarnings(stats::chisq.test(x, correct = FALSE)$statistic)/ (sum(x) * sqrt(prod(dim(x)-1)))))

}### end -- tschuprowT function
NULL




#' @title Bartels Rank Test of Randomness
#'
#' @description Performs Bartels rank test of randomness. The default method for testing the null hypothesis of randomness is \code{two.sided}. By using the alternative \code{left.sided}, the null hypothesis is tested against a trend. By using the alternative \code{right.sided}, the null hypothesis of randomness is tested against a systematic oscillation in the observed data.
#'
#' @param x A numeric vector of data values.
#' @param alternative A method for hypothesis testing, must be one of "two.sided" (default), "left.sided" or "right.sided".
#' @param pvalue A method for asymptotic aproximation used to compute the p-value.
#' @details Missing values are by default removed.
#' The RVN test statistic is
#' \deqn{RVN=\frac{\sum_{i=1}^{n-1}(R_i-R_{i+1})^2}{\sum_{i=1}^{n}\left(R_i-(n+1)/2\right)^2}}{RVN=\sum(R_i-R_{i+1})^2 / \sum(R_i-(n+1)/2)^2}
#' where \eqn{R_i=rank(X_i), i=1,\dots, n}{R_i=rank(X_i), i=1,...,n}. It is known that \eqn{(RVN-2)/\sigma} is asymptotically standard normal, where \eqn{\sigma^2=\frac{4(n-2)(5n^2-2n-9)}{5n(n+1)(n-1)^2}}{\sigma^2=[4(n-2)(5n^2-2n-9)]/[5n(n+1)(n-1)^2]}.
#'
#' @return
#' \item{statistic}{ The value of the RVN statistic test and the theoretical mean value and variance of the RVN statistic test.}
#' \item{n}{ the sample size, after the remotion of consecutive duplicate values.}
#' \item{p.value}{the asymptotic p-value.}
#' \item{method}{a character string indicating the test performed.}
#' \item{data.name}{a character string giving the name of the data.}
#' \item{alternative}{a character string describing the alternative.}
#'
#' @keywords Tests
#' @references
#' Bartels, R. (1982). The Rank Version of von Neumann's Ratio Test for Randomness, \emph{Journal of the American Statistical Association}, \bold{77}(377), 40-46.
#'
#' Gibbons, J.D. and Chakraborti, S. (2003). \emph{Nonparametric Statistical Inference}, 4th ed. (pp. 97-98). URL: \url{http://books.google.pt/books?id=dPhtioXwI9cC&lpg=PA97&ots=ZGaQCmuEUq}
#'
#' @examples
#' # Example 5.1 in Gibbons and Chakraborti (2003), p.98.
#' # Annual data on total number of tourists to the United States for 1970-1982.
#'  years <- 1970:1982
#' tourists <- c(12362, 12739, 13057, 13955, 14123,  15698, 17523,
#'  18610, 19842, 20310, 22500, 23080, 21916)
#'
#'  # See it graphically
#'  qplot(factor(years), tourists)+ geom_point()
#'
#' # Test the null against a trend
#'  bartels.rank(tourists, alternative="left.sided", pvalue="beta")
#' @export
`bartels.rank` <- function(x, alternative="two.sided", pvalue="normal") UseMethod("bartels.rank")


#' @rdname bartels.rank
#' @export
`bartels.rank.default` <- function(x, alternative="two.sided", pvalue="normal"){
  dname <- deparse(substitute(x))
  # Remove NAs
  x <- na.omit(x)
  stopifnot(is.numeric(x))
  n <- length(x)
  if (alternative == "t"){alternative <- "two.sided"}
  if (alternative == "l"){alternative <- "left.sided"}
  if (alternative == "r"){alternative <- "right.sided"}
  if (alternative != "two.sided" & alternative != "left.sided" & alternative != "right.sided")
  {stop("must give a valid alternative")}
  if (n < 10){stop("sample size must be greater than 9")}
  # unique
  rk <- rank(x)
  d <- diff(rk)
  #d.rank <- n*(n^2-1)/12
  d.rank <- sum(rk^2)-n*(mean(rk)^2)
  RVN <- sum(d^2)/d.rank
  mu <- 2
  vr <- (4*(n-2)*(5*n^2-2*n-9))/(5*n*(n+1)*(n-1)^2)
  # Computes the p-value
  if (pvalue == "auto"){pvalue<-ifelse(n<=100,"beta","normal")}
  if (pvalue == "beta"){
    btp <- (5*n*(n+1)*(n-1)^2)/(2*(n-2)*(5*n^2-2*n-9))-1/2
    pv0 <- stats::pbeta(RVN/4,shape1=btp,shape2=btp)
  }
  if (pvalue=="normal"){
    pv0 <- stats::pnorm((RVN - mu) / sqrt(vr))
  }
  if (alternative=="two.sided"){
    pv <- 2*min(pv0,1-pv0)
    alternative<-"nonrandomness"
  }
  if (alternative=="left.sided"){
    pv <- pv0
    alternative<-"trend"
  }
  if (alternative=="right.sided"){
    pv <- 1-pv0
    alternative<-"systematic oscillation"
  }
test <- (RVN - mu) / sqrt(vr)
rval <- list(statistic = c(statistic=test), nm=sum(d^2), rvn=RVN, mu=mu, var=vr, p.value = pv,
method = "Bartels Ratio Test", data.name = dname, parameter=c(n=n), n=n, alternative=alternative)
  class(rval) <- "htest"
  return(rval)
} ### end -- bartels.rank function
NULL





#' @title The Uncertainty Coefficient
#'
#' @description The uncertainty coefficient \code{U(C|R)} measures the proportion of uncertainty (entropy) in the column variable \code{Y} that is explained by the row variable \code{X}.
#'
#' @param x A numeric vector, a factor, matrix or data frame.
#' @param y A vector that is ignored if x is a matrix and required if x is a vector.
#' @param direction The direction of the calculation, either \code{"symmetric"} (default), \code{"row"}, or \code{"column"}. \code{"row"} calculates \code{uncertainty(R|C)} (column dependent relationship).
#' @param conf.level The confidence level of the interval. If set to NA (which is the default) no confidence interval will be calculated.
#' @param p.zero.correction Slightly nudge zero values so that their logarithm can be calculated.
#' @param \dots Further arguments are passed to the function \code{\link{table}}, allowing i.e. to set useNA. This refers only to the vector interface.
#'
#' @details The uncertainty coefficient is computed as
#' \deqn{U(C|R) = \frac{H(X) + H(Y) - H(XY)}{H(Y)} } and ranges from [0, 1].
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com},
#'  strongly based on code from Antti Arppe \email{antti.arppe@@helsinki.fi} and Andri Signorell \email{andri@@signorell.net}.
#'
#' @references
#' Theil, H. (1972), \emph{Statistical Decomposition Analysis}, Amsterdam: North-Holland Publishing Company.
#'
#' @keywords Multivariate
#'
#' @examples
#'  if (interactive()) {
#' # example from Goodman Kruskal (1954)
#' m <- as.table(cbind(c(1768,946,115), c(807,1387,438), c(189,746,288), c(47,53,16)));
#' dimnames(m) <- list(paste("A", 1:3), paste("B", 1:4));
#' print(m)
#'
#' calc.UC(m); # default is direction = "symmetric"
#'
#' calc.UC(m, conf.level=0.95); # direction "symmetric"
#'
#' calc.UC(m, direction="column");
#' }
#' @export
#' @rdname calc.UC
`calc.UC` <- function(x, y = NULL, direction = c("symmetric", "row", "column"), conf.level = NA, p.zero.correction = 1/sum(x)^2, ...) UseMethod("calc.UC")

#' @rdname calc.UC
#' @export
`calc.UC.default` <- function(x, y = NULL, direction = c("symmetric", "row", "column"), conf.level = NA, p.zero.correction = 1/sum(x)^2, ... ) {
  # Theil's UC (1970)
  # slightly nudge zero values so that their logarithm can be calculated (cf. Theil 1970: x->0 => xlogx->0)
  if(!is.null(y)) x <- table(x, y, ...)
  x[x == 0] <- p.zero.correction
  n <- sum(x)
  rsum <- apply(x, 1, sum)
  csum <- apply(x, 2, sum)
  hx <- -sum((apply(x, 1, sum) * log(apply(x, 1, sum)/n))/n)
  hy <- -sum((apply(x, 2, sum) * log(apply(x, 2, sum)/n))/n)
  hxy <- -sum(apply(x, c(1, 2), sum) * log(apply(x, c(1, 2), sum)/n)/n)
  switch( match.arg( arg = direction, choices = c("symmetric", "row", "column") )
          , "symmetric" = { res <- 2 * (hx + hy - hxy)/(hx + hy) }
          , "row" = { res <- (hx + hy - hxy)/hx }
          , "column" = { res <- (hx + hy - hxy)/hy }
  )

  if(!is.na(conf.level)){
    var.uc.RC <- var.uc.CR <- 0
    for(i in 1:nrow(x))
      for(j in 1:ncol(x))
      { var.uc.RC <- var.uc.RC + x[i,j]*(hx*log(x[i,j]/csum[j])+((hy-hxy)*log(rsum[i]/n)))^2/(n^2*hx^4);
      var.uc.CR <- var.uc.CR + x[i,j]*(hy*log(x[i,j]/rsum[i])+((hx-hxy)*log(csum[j]/n)))^2/(n^2*hy^4);
      }
    switch( match.arg( arg = direction, choices = c("symmetric", "row", "column") )
            , "symmetric" = {
              sigma2 <- 4*sum(x * (hxy * log(rsum %o% csum/n^2) - (hx+hy)*log(x/n))^2 ) /
                (n^2*(hx+hy)^4)
            }
            , "row" = { sigma2 <- var.uc.RC }
            , "column" = { sigma2 <- var.uc.CR }
    )
    pr2 <- 1 - (1 - conf.level)/2
    ci <- stats::qnorm(pr2) * sqrt(sigma2) * c(-1, 1) + res
    retval <- c(
      "Est. Mean" = res,
      "CI lower"=max(ci[1], -1),
      "CI upper"=min(ci[2], 1))
  }
  return(retval)
}### end -- uncertainty function
NULL







#' @title Cronbach's Alpha for a matrix or data frame
#' @name cronbach
#' @description This function calculates the Cronbach's alpha
#'  value of a data frame or matrix.
#'
#'
#'
#' @param df \code{data.frame} or matrix with more than 2 columns.
#' @return The Cronbach's alpha value for \code{df}.
#'
#' @importFrom stats na.omit var
#' @export
`cronbach` <- function(df) {
  df <- stats::na.omit(df)
  if (is.null(ncol(df)) || ncol(df) < 2) {
    warning("Too less columns in this factor to calculate alpha value!", call. = F)
    return(NULL)
  }
  return(dim(df)[2] / (dim(df)[2] - 1) * (1 - sum(apply(df, 2, var)) / stats::var(rowSums(df))))
}
NULL







#' Returns significance level.
#'
#' Returns the significance level as stars, or NA if a non-numeric value is
#' passed in.
#'
#' @param x p-value.
#' @export
star <- function(x) {
  if(is.numeric(x)) {
    return(cut(x,
               breaks=c(-Inf,.001,.01,.05,.1,Inf),
               labels=c('***','**','*','.',''))
    )
  } else {
    return(NA)
  }
}
NULL

