## The example normal model for bdpopt

## Create an instance of the special normal model.
## theta.mu and theta.tau are mean and precision parameters of the priors for the parameters of the efficacy Emax model. Correspondingly for eta.
## n.II is the vector of phase II sample sizes.
## d.II is the vector of phase II doses. Length must be equal to length of n.II.
## sigmaE and sigmaS are the known standard deviations of the efficacy and safety responses.
create.normal.model <-
    function(theta.mu, theta.tau,
             eta.mu, eta.tau,
             n.II, d.II,
             YE.II, YS.II,
             sigmaE, sigmaS,
             k.III, path.to.package = NA) {        
    
    ## Contruct the specific data list for the normal model
    k.II <- length(n.II)
   
    data <- list(k.II = k.II,                 
                 theta.mu = theta.mu, theta.tau = theta.tau,
                 eta.mu = eta.mu, eta.tau = eta.tau,
                 n.II = n.II, d.II = d.II,
                 YE.II = YE.II, YS.II = YS.II,
                 sigmaE = sigmaE, sigmaS = sigmaS,
                 k.III = k.III)

    path.to.package <-
        if (is.na(path.to.package)) {
            if (!("package:bdpopt" %in% search()))
                stop("package 'bdpopt' must be attached if 'path.to.package' is not provided")
            
            path.package("bdpopt")
            
        } else {
            if ( !(is.character(path.to.package) &&
                   length(path.to.package) == 1) )
                stop("'path.to.package' must be an atomic character vector containing a single string")
            
            path.to.package
    }

    sim.model(paste0(path.to.package, "/extdata/normal_model_jags_model.R"), data)
}

## Create an instance of the special normal model directly from a data file.
create.normal.model.from.file <- function(path.to.package = NA) {
    path.to.package <-
        if (is.na(path.to.package)) {
            if (!("package:bdpopt" %in% search()))
                stop("package 'bdpopt' must be attached if 'path.to.package' is not provided")
            
            path.package("bdpopt")
            
        } else {
            if ( !(is.character(path.to.package) &&
                   length(path.to.package) == 1) )
                stop("'path.to.package' must be an atomic character vector containing a single string")
            
            path.to.package
    }    
    
    sim.model(paste0(path.to.package, "/extdata/normal_model_jags_model.R"),
              paste0(path.to.package, "/extdata/normal_model_jags_data.R"))
}

## Create utility function for use together with the normal model.
## sig.level is the one-sided significance level for the RA significance tests.
create.utility.function <- function(model,
                                    n.min, sig.level, safety.max,
                                    cE, cS, p,
                                    fixed.cost, cost.per.sample) {

    ## Check that the model has been constructed and that sig.level is of acceptable type
    if (!inherits(model, "sim.model"))
        stop("'model' must be of class \"sim.model\"")

    if ( !(is.numeric(sig.level) &&
           length(sig.level) == 1 &&
           0 < sig.level && sig.level < 1) )
        stop("'sig.level' must be a number strictly between 0 and 1")
    
    z.level <- qnorm(1 - sig.level)    
    sigmaE <- model$data$sigmaE
    sigmaS <- model$data$sigmaS
    k.III <- model$data$k.III
    
    ## Model of the decision of the regulatory authority (RA).
    ## Approve for marketing if and only if,
    ## (i) the sample size is at least n.min, and,
    ## (ii) all trials show acceptable results for both efficacy and safety.  
    RA.decision <- function(YE, YS, n) {
        if (n >= n.min &&
            YE > (z.level * sigmaE / sqrt(n)) &&
            YS < safety.max - (z.level * sigmaS / sqrt(n)))            
            1
        else
            0
    }          

    gain <- function(YE, YS, muE, muS) {     
        p * (cE * mean(YE) + cS * mean(YS)) + (1 - p) * (cE * mean(muE) + cS * mean(muS))
    }
    
    trial.cost <- function(n) {
        fixed.cost + cost.per.sample * k.III * n
    }
    
    function(YE.III, YS.III, muE.III, muS.III, n.III) {    
        RA.decision(YE.III, YS.III, n.III[1]) * gain(YE.III, YS.III, muE.III, muS.III) - trial.cost(n.III[1])
    }
}
