quasibin <- function(formula, data, link = c("logit", "cloglog"), phi = NULL, tol = 0.001){
# utility function to remove white spaces at the beginning and at the end of character strings
  tr <- function(string) gsub("^[[:space:]]+|[[:space:]]+$", "", string)
## check call validity
  CALL <- match.call(expand.dots = FALSE)
  link <- match.arg(link)
  f <- formula
  if(length(f) != 3)
    stop(paste(tr(deparse(formula)), collapse = " "), "is not a valid formula.")
  else
    if(substring(deparse(f), 1, 5) != "cbind")
      stop(paste(tr(deparse(formula)), collapse = " "), " is not a valid formula.\n",
           "The response must be a matrix of the form cbind(success, failure)")
  mf <- model.frame(formula = f, data = data)
  Y <- model.response(mf)
  n <- rowSums(Y)
  fam <- eval(parse(text = paste("binomial(link =", link,")")))
  if(is.null(phi)){
    phi <- 1e-04
    X2 <- 0
    delta <- X2 + 2 * tol
    w <- rep(1, length(n))
    dfr <- data.frame(data, w = w)
    fm <- glm(formula = f, family = fam, weights = w, data = dfr)
    ok <- TRUE
    while(ok){
      X2 <- sum(residuals(fm, type = "pearson")^2)
      delta <- X2 - df.residual(fm)
      if(delta <= tol)
        ok <- FALSE
      else{
        phi <- phi * sum(residuals(fm, type = "pearson")^2) / df.residual(fm)
        w <- 1 / (1 + (n - 1) * phi)
        dfr <- data.frame(data, w = w)
        fm <- glm(formula = f, family = fam, weights = w, data = dfr)
        }
      }
    }
  #results
  w <- 1 / (1 + (n - 1) * phi)
  dfr <- data.frame(data, w = w)
  fm <- glm(formula = formula, family = fam, weights = w, data = dfr)
  # outputs
  new(Class = "glimQL", CALL = CALL, fm = fm, phi = phi)
  }
