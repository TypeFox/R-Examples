epi.detectsize <- function(N, prev, se, sp, interpretation = "series", covar = c(0,0), conf.level = 0.95, finite.correction = TRUE){
    
    # Covar is a vector of length two. First element is covariance for D+ group, second element is covariance for D- group.
    # See Dohoo, Martin and Stryhn page 103.
    alpha <- (1 - conf.level)
        
    # Work out sensitivity and specificity:
    if (length(se) > 1 & interpretation == "series") {
        Ses <- se[1] * se[2] + covar[1]
        Sps <- 1 - (1 - sp[1]) * (1 - sp[2]) - covar[2]
        use <- Ses
        usp <- Sps
    }
    if (length(se) > 1 & interpretation == "parallel") {
        Sep <- 1 - (1 - se[1]) * (1 - se[2]) - covar[1]
        Spp <- sp[1] * sp[2] + covar[2]
        use <- Sep
        usp <- Spp
    }
    if (length(se) == 1) {
        use <- se
        usp <- sp
    }
    if (length(N) == 1) {
        units = round((1 - alpha^(1/(N[1] * prev[1] * use))) * (N[1] - (N[1] * prev[1] * use - 1)/2), digits = 0)
        units.corrected <- round(units/(1 + (units/N[1])), digits = 0)
        if (finite.correction == TRUE) {
            performance <- as.data.frame(cbind(sens = use, spec = usp))
            sample.size <- units.corrected
            rval <- list(performance = performance, sample.size = sample.size)
        }
        
        if (finite.correction == FALSE) {
            performance <- as.data.frame(cbind(sens = use, spec = usp))
            sample.size <- units
            rval <- list(performance = performance, sample.size = sample.size)
        }
    }
    if (length(N) == 2) {
        units <- round((1 - alpha^(1/(N[2] * prev[2] * use))) * (N[2] - (N[2] * prev[2] * use - 1)/2), digits = 0)
        pd <- prev[1] * (1 - alpha)
        clusters <- round((1 - alpha^(1/(N[1] * pd))) * (N[1] - 
            (N[1] * pd - 1)/2), digits = 0)
        total <- units * clusters
        units.corrected <- round(units/(1 + (units/N[2])), digits = 0)
        clusters.corrected <- round(clusters/(1 + (clusters/N[1])), 
            digits = 0)
        total.corrected <- units.corrected * clusters.corrected
        if (finite.correction == TRUE) {
            performance <- as.data.frame(cbind(sens = use, spec = usp))
            sample.size <- as.data.frame(cbind(clusters = clusters.corrected, 
                units = units.corrected, total = total.corrected))
            rval <- list(performance = performance, sample.size = sample.size)
        }
        if (finite.correction == FALSE) {
            performance <- as.data.frame(cbind(sens = use, spec = usp))
            sample.size <- as.data.frame(cbind(clusters = clusters, 
                units = units, total = total))
            rval <- list(performance = performance, sample.size = sample.size)
        }
    }
    return(rval)
}
