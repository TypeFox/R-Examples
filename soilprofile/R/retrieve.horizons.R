retrieve.horizons <-
function(profile_data, a, col, lwd, horizon.border=horizon.border) {
  horizon_retrieved <- profile_data[[a]][[1]]
  polygon(horizon_retrieved, col=col, lwd=lwd, border=horizon.border)
}
