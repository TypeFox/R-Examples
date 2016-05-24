PD_internal = function(demand = demand, tau = 0.5, eps = 0.005) {
    # CHECK IF msm LIBRARY EXISTS
 #   library(msm)
    # M - NUMBER OF CONSTRAINED DEMAND INSTANCES; N - NUMBER OF UNCONSTRAINED DEMAND INSTANCES
    M <- length(demand[names(demand) == "1"])
    if (M == 0) {stop("Warning: All demand instances are unconstrained")}
    else {
        N <- length(demand[names(demand) == "0"])
        # INITIALIZATION: COMPUTE THE MEAN AND THE STANDARD DEVIATION USING THE UNCONSTRAINED DEMAND VALUES
        mean_v <- sum(demand[names(demand) == "0"])/N
        sd_v <- sqrt(sum((demand[names(demand) == "0"] - mean_v)^2)/N)
        history <- as.data.frame(matrix(c(mean_v, sd_v), 1,2))
        # SPECIFY WHERE THE DISTRIBUTIONS ARE LEFT TRUNCATED
        bl_const <- demand[names(demand) == "1"]
        # INITIALIZE NUMBER OF ITERATIONS
        niter <- 1
        repeat {
            # EXPECTATION: COMPUTE EXPECTED VALUES OF THE TRUNCATED NORMAL DISTRIBUTIONS
            tmp <- numeric()
            for (ind_constr in 1:length(bl_const)) {
                tmp <- c(tmp, qtnorm(tau, mean = mean_v, sd = sd_v, lower = bl_const[ind_constr], upper = Inf, lower.tail = F))
            } # end for ind_constr
            demand[names(demand) == "1"] <- tmp
            # MAXIMIZATION: COMPUTE THE NEW VALUES FOR THE MEAN AND THE STANDARD DEVIATION
            new_mean_v <- mean(demand)
            new_sd_v <- sd(demand)
            if (sum(c(abs(new_mean_v - mean_v), abs(new_sd_v - sd_v)) > eps) == 0 ){
                niter <- niter + 1
                history <- rbind(history, c(new_mean_v, new_sd_v))
                break
                }
            else {
                niter <- niter + 1
                history <- rbind(history, c(new_mean_v, new_sd_v))
                mean_v <- new_mean_v
                sd_v <- new_sd_v
            }
        } # end repeat
        names(history) <- c("Mean", "StDev")
        return(list(param = round(c(new_mean_v, new_sd_v), 2), niter = niter, demand = round(demand, 2), history = round(history, 2)))
    } # end if
} # end PD_internal function
