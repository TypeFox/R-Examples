# return the sun of square of differences between the two models

.fitSSM <- function(par, temperatures, growth.rate) {
  growth.rate2 <- .SSM(273.15+temperatures, par)[[1]]*1E5
  return(sum((growth.rate-growth.rate2)^2, na.rm = TRUE))
}

