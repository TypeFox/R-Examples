#' Power-time modelling.
#'
#' Model the Power-time (Pt) relationship for a set of data. This is done via
#' nonlinear least squares regression of four models: an inverse model; an
#' exponential model; a bivariate power function model; and a three parameter
#' inverse model. An S3 object of class "Ptmodels" is returned, which currently
#' has methods for \link[base]{print}, \link[stats]{coef}, \link[base]{summary},
#' and \link[stats]{predict}. If inputs do not conform well to the models, a
#' warning message is generated. This function will make use of
#' \code{minpack.lm::nlsLM} if available.
#'
#' @param P a numeric vector of maximal mean power values for time periods given
#'   in the \code{tsec} argument.
#' @param tsec a numeric vector of time values that (positionally) correspond to
#'   elements in \code{P}.
#'
#' @return returns an S3 object of class "Ptmodels".
#'
#' @seealso \code{\link{predict.Ptmodels}}
#'
#' @references R. Hugh Morton (1996) A 3-parameter critical power model,
#'   Ergonomics, 39:4, 611-619,
#'   \href{http://dx.doi.org/10.1080/00140139608964484}{DOI}.
#'
#' @examples
#' data(Pt_prof)  # Example power-time profile.
#'
#' P    <- Pt_prof$pwr
#' tsec <- Pt_prof$time
#'
#' mdls <- Pt_model(P, tsec)  # Model.
#' print(mdls)
#'
#' coef(mdls)
#' summary(mdls)
#'
#' @export
Pt_model <- function(P, tsec) {
  if (length(P) != length(tsec))
    stop("Mismatching values.", call. = FALSE)
  #  ---------------------------------------------------------------------------
  if (length(P) == 2) {
    m <- CP_model_inv(P, tsec)
    params <- data.frame(CP = m["CP"], "W'" = m["W"], row.names = "Inverse model")
    return(params)
  }
  #  ---------------------------------------------------------------------------
  models <- Pt_nls(requireNamespace("minpack.lm", quietly = TRUE), P, tsec)
  # Assemble S3 object ---------------------------------------------------------
  out <- list(
    models = models,
    table  = Pt_table(models),
    Pfn    = Pt_fn(models, "P"),
    tfn    = Pt_fn(models, "tsec")
  )
  class(out) <- "Ptmodels"
  out
}

# Helper functions -------------------------------------------------------------
CP_model_inv <- function(P, tsec) {
  m   <- lm(P ~ {1 / tsec})
  out <- m$coefficients
  names(out) <- c("CP", "W")
  out
}

Pt_nls <- function(LM = FALSE, P, tsec) {  # Generate nls model objects.
  algor    <- ifelse(LM, "LM", "port")
  parseval <- function(x) eval(parse(text = x))
  maxit    <- function() "list(maxiter = 1000)"
  margs    <- list(
    inv  = list(formula = "P ~ {(W / tsec) + CP}", control = maxit(),
                start = "c(W  = 20000, CP = 300)",
                upper = "c(W  = 60000, CP = 600)", lower = "c(W  = 0, CP = 0)"),
    exp  = list(formula = "P ~ {a * exp(1) ^ (k * tsec) + CP}", control = maxit(),
                start = "c(a = Diff(range(P)), k = -0.005, CP = min(P))",
                upper = "c(a = 1000, k = 0, CP = 1000)", lower = "c(a = 0, k = -1, CP = 0)"),
    pwr  = list(formula = "P ~ {k * (tsec ^ n) + A}", control = maxit(),
                start = "c(k = 1500, n = -0.5, A = 0)"),
    thrp = list(formula = "P ~ {(W / (tsec - k)) + CP}", control = maxit(),
                start = "c(W = 20000, CP = 300, k = -20)",
                upper = "c(W = Inf, CP = 600, k = 0)", lower = "c(W = 0, CP = 0, k = -1000)")
  )
  lapply(margs, function(a)
    do.call(ifelse(LM, minpack.lm::nlsLM, nls), args = lapply(a, parseval)))
}

Pt_table <- function(m) { # Generate table from model objects.
  coeff <- function(i) unname(round(coef(m[[i]]), 3))

  f   <- c(
    paste(coeff("inv"), collapse = " / x + "),
    paste0(coeff("exp")[1], " * e^", coeff("exp")[2], "x", " + ", coeff("exp")[3]),
    paste0(coeff("pwr")[[1]], " * x^", coeff("pwr")[[2]],
           ifelse(coeff("pwr")[[3]] < 0, " - ", " + "),
           abs(coeff("pwr")[[3]])),
    paste0(coeff("thrp")[1], " / (x",
           ifelse(coeff("thrp")[3] < 0, " + ", " - " ),
           abs(coeff("thrp")[3]), ") + ", coeff("thrp")[2])
  )

  RSE <- vapply(m, function(x) summary(x)$sigma, numeric(1))
  tab <- data.frame(formula = f, RSE  = RSE, fit  = NA)
  # Annotate.
  bestfitrow <- which.min(tab$RSE)  # Residual Standard Error.
  tab$fit[bestfitrow] <- "**"
  tab$fit[is.na(tab$fit)] <- " "
  rownames(tab) <- c("Inverse", "Exponential", "Power", "Three-param")
  # Warning message if models go haywire.
  CP  <- c(coef(m$inv)[[2]], coef(m$exp)[[3]], coef(m$thrp)[[2]])
  if (any(CP < 0))
    warning("Inappropriate data for these models, revise inputs.", call. = FALSE)

  tab
}

Pt_fn <- function(m, y = "P") { # Prediction functions.
  switch(
    y,
    "P" = list(
      inv = function(x)
        unname({coef(m$inv)["W"] / x + coef(m$inv)["CP"]}),
      exp = function(x)
        unname({coef(m$exp)["a"] * exp(1) ^ (coef(m$exp)["k"] * x) + coef(m$exp)["CP"]}),
      pwr = function(x)
        unname({coef(m$pwr)["k"] * (x ^ coef(m$pwr)["n"]) + coef(m$pwr)["A"]}),
      thrp = function(x)
        unname({(coef(m$thrp)["W"] / (x - coef(m$thrp)["k"])) + coef(m$thrp)["CP"]})
    ),
    "tsec" = list(
      inv = function(x)
        unname({coef(m$inv)["W"] / (x - coef(m$inv)["CP"])}),
      exp = function(x)
        unname({(1 / coef(m$exp)["k"]) * log((x - coef(m$exp)["CP"]) / coef(m$exp)["a"])}),
      pwr = function(x)
        unname({((x - coef(m$pwr)["A"]) / coef(m$pwr)["k"]) ^ (1 / coef(m$pwr)["n"])}),
      thrp = function(x)
        unname({coef(m$thrp)["W"] / (x - coef(m$thrp)["CP"]) + coef(m$thrp)["k"]})
    )
  )
}
