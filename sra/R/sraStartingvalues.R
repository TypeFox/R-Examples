sraStartingvalues <-
function (parameter, sradata, rand = 0) 
{
    "sraStartingvalues.mu0" <- function(sradata) {
        return(mean(sradata[sradata[, "gen"] == 1, "mean"]) + 
            rand * rnorm(1, 0, abs(mean(sradata[sradata[, "gen"] == 
                1, "mean"]))))
    }
    "sraStartingvalues.logvarA0" <- function(sradata) {
        return(log(0.2 * mean(sradata[, "var"] * exp(rand * rnorm(1, 
            0, 0.2 * mean(sradata[, "var"]))))))
    }
    "sraStartingvalues.logvarE0" <- function(sradata) {
        return(log(0.8 * mean(sradata[, "var"] * exp(rand * rnorm(1, 
            0, 0.8 * mean(sradata[, "var"]))))))
    }
    if (parameter == "mu0") {
        return(sraStartingvalues.mu0(sradata))
    }
    if (parameter == "logvarA0") {
        return(sraStartingvalues.logvarA0(sradata))
    }
    if (parameter == "logvarE0") {
        return(sraStartingvalues.logvarE0(sradata))
    }
    if (parameter == "logvarME") {
        return(log(exp(sraStartingvalues.logvarE0(sradata))/2))
    }
    if (parameter == "logIA0") {
        return(log(exp(sraStartingvalues.logvarA0(sradata))/(sraStartingvalues.mu0(sradata)^2)))
    }
    if (parameter == "logIE0") {
        return(log(exp(sraStartingvalues.logvarE0(sradata))/(sraStartingvalues.mu0(sradata)^2)))
    }
    if (parameter == "logith20") {
        a <- exp(sraStartingvalues.logvarA0(sradata))
        e <- exp(sraStartingvalues.logvarE0(sradata))
        h2 <- a/(a + e)
        return(log(h2/(1 - h2)))
    }
    if (parameter == "logvarP0") {
        return(log(exp(sraStartingvalues.logvarA0(sradata)) + 
            exp(sraStartingvalues.logvarE0(sradata))))
    }
    if (parameter == "o") {
        return(sraStartingvalues.mu0(sradata))
    }
    if (parameter == "s") {
        return(1 + rand * rnorm(1, 0, 1))
    }
    if (parameter == "logepsilon") {
        return(-10 + rand * rnorm(1, 0, 10))
    }
    if (parameter == "logminusepsilon") {
        return(-10 + rand * rnorm(1, 0, 10))
    }
    if (parameter == "logvarepsilon") {
        return(0 + rand * rnorm(1, 0, 1))
    }
    if (parameter == "kc") {
        return(0 + rand * rnorm(1, 0, 1))
    }
    if (parameter == "kg") {
        return(0 + rand * rnorm(1, 0, 1))
    }
    if (parameter == "logvarM") {
        return(-20 + rand * rnorm(1, 0, 10))
    }
    if (parameter == "logNe") {
        return(log(100) + rand * rnorm(1, 0, 2))
    }
    if (parameter == "relativekA0") {
        return(0 + rand * rnorm(1, 0, 0.5))
    }
    if (parameter == "relativekE0") {
        return(0 + rand * rnorm(1, 0, 0.5))
    }
    if (parameter == "kA1") {
        return(1 + rand * rnorm(1, 0, 0.5))
    }
    if (parameter == "kE1") {
        return(1 + rand * rnorm(1, 0, 0.5))
    }
    if (parameter == "kA2") {
        return(0 + rand * rnorm(1, 0, 0.5))
    }
    if (parameter == "kE2") {
        return(0 + rand * rnorm(1, 0, 0.5))
    }
    if (parameter == "kA3") {
        return(0 + rand * rnorm(1, 0, 0.5))
    }
    if (parameter == "kE3") {
        return(0 + rand * rnorm(1, 0, 0.5))
    }
    if (parameter == "logrelativekA0") {
        return(-1 + rand * rnorm(1, 0, 1))
    }
    if (parameter == "logrelativekE0") {
        return(-1 + rand * rnorm(1, 0, 1))
    }
    if (parameter == "logkA1") {
        return(0 + rand * rnorm(1, 0, 1))
    }
    if (parameter == "logkE1") {
        return(0 + rand * rnorm(1, 0, 1))
    }
    if (parameter == "logkA2") {
        return(-3 + rand * rnorm(1, 0, 2))
    }
    if (parameter == "logkE2") {
        return(-3 + rand * rnorm(1, 0, 2))
    }
    if (parameter == "logkA3") {
        return(-3 + rand * rnorm(1, 0, 2))
    }
    if (parameter == "logkE3") {
        return(-3 + rand * rnorm(1, 0, 2))
    }
    stop("Unknown parameter ", parameter, ".")
}
