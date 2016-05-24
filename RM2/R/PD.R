PD = function (demand = demand, tau = 0.5, eps = 0.005) {
    res <- try(RM2:::PD_internal(demand, tau, eps), silent = TRUE)
    if (class(res) == 'try-error') {print("Warning: All demand instances are unconstrained")} else {return(res)}
} # end PD function
