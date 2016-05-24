FillEnvelope <-
function(Envelope, Alpha, Global) {

  # Sims contains simulated values (one line per simulation)
  Sims <- attr(Envelope, "simfuns")

  # Local or Global
  if (Global) {
    Env <- GlobalEnvelope(Sims, Alpha)
    attr(Envelope, "desc")[4] <- "lower global envelope of %s from simulations"
    attr(Envelope, "desc")[5] <- "upper global envelope of %s from simulations"
    attr(Envelope, "einfo")$global <- TRUE
  } else {
    Env <- apply(as.data.frame(Sims)[, -1], 1, stats::quantile, probs=c(Alpha/2, 1-Alpha/2), na.rm=TRUE)
    # quantile returns 0 when a value is NA. Force NA to appear in the result: add NA or 0.
    Env <- t(t(Env) + apply(as.data.frame(Sims), 1, sum)*0)
    attr(Envelope, "einfo")$nrank <- attr(Envelope, "einfo")$nsim*Alpha/2
  }
  # Computed values are injected into the envelope object
  Envelope$lo <- Env[1,]
  Envelope$hi <- Env[2,]
  # Prepare summary
  attr(Envelope, "einfo")$Alpha <- Alpha
  class(Envelope) <- c("dbmssEnvelope", class(Envelope))
  return(Envelope)
}