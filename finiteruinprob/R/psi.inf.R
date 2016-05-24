## Use this library to compute the exact probability of ruin in the infinite
## time horizon
psi.inf <- function(model) {
    stopifnot(names(model) %in% c('premium', 'freq', 'variance', 'rates'))
    sdprisk::hypoexpRuinprob(riskproc(claiminfo(hypoexp = list(rates = model$rates)),
                                      premium  = model$premium,
                                      freq     = model$freq,
                                      variance = model$variance))[['psi']]
}
