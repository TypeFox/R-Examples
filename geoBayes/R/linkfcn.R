##' Link function for the exponential family.
##'
##' \code{linkfcn} maps the mean of the response variable \code{mu} to
##' the linear predictor \code{z}. \code{linkinv} is its inverse.
##' 
##' Note that the logit link for the
##' binomial family is defined as the quantile of the logistic
##' distribution with scale 0.6458.
##'
##' For the Gaussian family, if the link parameter is positive, then
##' the extended link is used, defined by \deqn{z =
##' \frac{sign(\mu)|\mu|^\nu - 1}{\nu}}{z = (sign(mu)*abs(mu)^nu -
##' 1)/nu} In the other case, the link function is the same as for the
##' Poisson and gamma families.
##'  
##' For the Poisson and gamma families, the Box-Cox transformation is
##' used, defined by \deqn{z = \frac{\mu^\nu - 1}{\nu}}{z = (mu^nu -
##' 1)/nu}
##'
##' For the GEV binomial family, the link function is defined by
##' \deqn{\mu = 1 - \exp\{-\max(0, 1 + \nu z)^{\frac{1}{\nu}}\}}{mu =
##' 1 - exp[-max(0, 1 + nu z)^(1/nu)]} for any real \eqn{\nu}{nu}. At
##' \eqn{\nu = 0}{nu = 0} it reduces to the complementary log-log
##' link.
##'
##' The Wallace binomial family is a fast approximation to the robit
##' family. It is defined as \deqn{\mu = 
##' \Phi(\mbox{sign}(z) c(\nu) \sqrt{\nu \log(1 + z^2/\nu)})}{mu = 
##' Phi(sign(z) c(nu) sqrt{nu log(1 + z^2/nu)})}
##' where \eqn{c(\nu) = (8\nu+1)/(8\nu+3)}{c(nu) = (8*nu+1)/(8*nu+3)}
##'
##' @title Calculate the link function for exponential families
##' @param mu Numeric. The mean of the response variable. 
##' @param z Numeric. The linear predictor. 
##' @param linkp The link function parameter. A scalar but for the
##' binomial family is also allowed to have the character values
##' "logit" or "probit".
##' @param family The distribution of the response variable.
##' @return A numeric array of the same dimension as the function's
##' first argument.
##' @seealso \code{\link{comparebinlinks}}
##' @examples \dontrun{
##' mu <- seq(0.1, 0.9, 0.1)
##' linkfcn(mu, 7, "binomial")       # robit(7) link function
##' linkfcn(mu, "logit", "binomial") # logit link function
##'
##' mu <- seq(-3, 3, 1)
##' linkfcn(mu, 0.5, "gaussian")     # sqrt transformation
##' linkinv(linkfcn(mu, 0.5, "gaussian"), 0.5, "gaussian")
##' curve(linkfcn(x, 0.5, "gaussian"), -3, 3)
##' }
##' @family linkfcn
##' @name linkfcn
##' @rdname linkfcn
##' @importFrom stats qlogis qnorm qt qf plogis pnorm pt 
##' @export 
linkfcn <- function (mu, linkp,
                     family = c("gaussian", "binomial", "poisson", "Gamma",
                       "GEV.binomial", "GEVD.binomial",
                                "Wallace.binomial")) {
  family <- match.arg(family)
  if (length(linkp) != 1) stop ("The linkp argument must be scalar")
  if (is.character(linkp) | is.factor(linkp)) {
    if (family == "binomial") {
      if (linkp == "logit") {
        z <- qlogis(mu, scale = 0.6458)
      } else if (linkp == "probit") {
        z <- qnorm(mu)
      } else {
        stop ("Cannot recognise character link for binomial")
      }
    } else stop ("Character link is only allowed for binomial")
  } else if (!is.numeric(linkp)) {
    stop ("Argument linkp must be numeric, or in the case of
the binomial can also be the character \"logit\" or \"probit\"")
  } else { # Numeric link parameter
    nu <- as.double(linkp)
    if (family == "gaussian") {
      if (nu == 0) {
        z <- ifelse(mu > 0, log(mu), -Inf)
      } else if (nu > 0) {
        z <- (sign(mu)*abs(mu)^nu - 1)/nu
      } else if (nu < 0) {
        z <- ifelse(mu > 0, (mu^nu - 1)/nu, -Inf)
      }
    } else if (family == "binomial") {
      if(nu <= 0) {
        stop ("The robit link parameter must be positive")
      } else if (nu < 1) {
        z <- sign(mu-.5)*sqrt(qf(abs(2*mu-1), 1, nu))
      } else {
        z <- qt(mu, nu)
      }
    } else if (family == "poisson" | family == "Gamma") {
      if (nu == 0) {
        z <- ifelse(mu > 0, log(mu), -Inf)
      } else {
        z <- ifelse(mu > 0, (mu^nu - 1)/nu, -Inf)
      }
    } else if (family == "GEV.binomial") {
      z <- -log1p(-mu)
      if (nu == 0) {
        z <- ifelse(mu > 0, ifelse(mu < 1, log(z), Inf), -Inf)
      } else {
        z <- ifelse(mu > 0, ifelse(mu < 1, (-1 + z^nu)/nu, Inf), -Inf)
      }
    } else if (family == "GEVD.binomial") {
      z <- -log1p(-mu)
      if (nu == 0) {
        z <- ifelse(mu > 0, ifelse(mu < 1, -log(z), -Inf), Inf)
      } else {
        z <- ifelse(mu > 0, ifelse(mu < 1, (-1 + z^(-nu))/nu, -Inf), Inf)
      }
    } else if (family == "Wallace.binomial") {
      if (nu <= 0)
        stop ("The link parameter for Wallace.binomial must be positive.")
      z <- 8*nu
      z <- qnorm(mu)*(z+3)/(z+1)
      z <- sign(mu-.5)*sqrt(nu*expm1(z*z/nu))
    } else {
      stop ("Unrecognised family")
    }
  }
  z
}

##' @family linkfcn
##' @rdname linkfcn
##' @export 
linkinv <- function (z, linkp, 
                     family = c("gaussian", "binomial", "poisson", "Gamma",
                       "GEV.binomial", "GEVD.binomial",
                       "Wallace.binomial")) {
  family <- match.arg(family)
  if (length(linkp) != 1) stop ("The linkp argument must be scalar")
  if (is.character(linkp) | is.factor(linkp)) {
    if (family == "binomial") {
      if (linkp == "logit") {
        mu <- plogis(z, scale = 0.6458)
      } else if (linkp == "probit") {
        mu <- pnorm(z)
      } else {
        stop ("Cannot recognise character link for binomial")
      }
    } else stop ("Character link is only allowed for binomial")
  } else if (!is.numeric(linkp)) {
    stop ("Argument linkp must be numeric, or in the case of
the binomial can also be the character \"logit\" or \"probit\"")
  } else { # Numeric link parameter
    nu <- as.double(linkp)
    if (family == "gaussian") {
      if (nu == 0) {
        mu <- exp(z)
      } else if (nu > 0) {
        mu <- nu*z + 1
        mu <- sign(mu)*abs(mu)^(1/nu)
      } else { # nu < 0
        mu <- nu*z + 1
        ifelse(mu > 0, mu^(1/nu), Inf)
      }
    } else if (family == "binomial") {
      if(nu <= 0) {
        stop ("The robit link parameter must be positive")
      } else {
        mu <- pt(z, nu)
      }
    } else if (family == "poisson" | family == "Gamma") {
      if (nu == 0) {
        mu <- exp(z)
      } else {
        mu <- nu*z
        mu <- ifelse(mu > -1, exp(log1p(mu)/nu), 0)
      }
    } else if (family == "GEV.binomial") {
      if (nu == 0) {
        mu <- 1 - exp(-exp(z))
      } else {
        mu <- 1 + nu*z
        mu <- ifelse(mu > 0, -expm1(-mu^(1/nu)), as.double(nu < 0))
      }
    } else if (family == "GEVD.binomial") {
      if (nu == 0) {
        mu <- 1 - exp(-exp(-z))
      } else {
        mu <- 1 + nu*z
        mu <- ifelse(mu > 0, -expm1(-mu^(-1/nu)), as.double(nu > 0))
      }
    } else if (family == "Wallace.binomial") {
      if (nu <= 0)
        stop ("The link parameter for Wallace.binomial must be positive.")
      mu <- 8*nu
      mu <- (mu+1)/(mu+3)*sqrt(nu*log1p(z*z/nu))
      mu <- ifelse(z > 0, pnorm(mu), pnorm(mu, lower.tail = FALSE))
    } else {
      stop ("Unrecognised family")
    }
  }
  mu
}
