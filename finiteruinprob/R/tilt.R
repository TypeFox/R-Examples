## Model parameters under the exponentially tilted measure (with tilting
## parameter theta)
tilt <- function(model, theta) {
    stopifnot(names(model) %in% c('premium', 'freq', 'variance', 'rates'))
    list(premium  = model$premium - theta * model$variance,
         freq     = model$freq * mgfhypoexp(theta, model$rates, 0L),
         variance = model$variance,
         rates    = model$rates - theta)
}
