setClass("initMethod",
         representation(step1 = "FLXcontrol",
                        step2 = "FLXcontrol"))

initMethod <- function(name = c("tol.em", "cem.em", "sem.em"),
                       step1 = list(tolerance = 10^-2), step2 = list(), control = list(), nrep = 3L) {
  name <- match.arg(name)
  z <- new("initMethod",
           step1 = as(c(step1, control), "FLXcontrol"),
           step2 = as(c(step2, control), "FLXcontrol"))
  z@step1@nrep <- as.integer(nrep)
  z@step2@nrep <- 1L
  z@step1@classify <- switch(name,
                             cem.em = "CEM",
                             sem.em = "SEM",
                             tol.em = "weighted")
  z
}

initFlexmix <- function(..., k, init = list(), control = list(), nrep = 3L, verbose = TRUE, drop = TRUE, unique = FALSE)
{
  MYCALL <- match.call()
  if (missing(k)) stop("'k' is missing.")
  if (!missing(control) & is(init, "initMethod")) warning("'control' argument ignored.")
  init <- do.call("initMethod", c(init, list(control = control, nrep = nrep)))
  MYCALL1 <- lapply(k, function(K) {
    MYCALL[["k"]] <- as.numeric(K)
    MYCALL
  })
  names(MYCALL1) <- paste(k)
  STEP1 <- stepFlexmix(..., k = k, verbose = verbose, drop = FALSE, unique = FALSE,
                       nrep = init@step1@nrep, control = init@step1)
  models <- lapply(k, function(K) {
    if (length(k) > 1 && verbose) cat("* ")
    new("flexmix",
        flexmix(..., control = init@step2,
                cluster = posterior(getModel(STEP1, paste(K)))),
        k0 = as.integer(K), call = MYCALL1[[paste(K)]])
  })
  if (length(k) > 1 && verbose) cat("\n")
  names(models) <- paste(k)
  if (drop & length(models) == 1) {
    return(models[[1]])
  } else {
    z <- new("stepFlexmix", models = models, k = as.integer(k), logLiks = STEP1@logLiks, nrep = STEP1@nrep, call = MYCALL)
    if (unique) 
      z <- unique(z)
    return(z)
  }
}
