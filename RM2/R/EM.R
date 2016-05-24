EM = function (demand = demand, eps = 0.005) {
    res <- try(RM2:::EM_internal(demand, eps), silent = TRUE)
    if (class(res) == 'try-error') {print("Warning: All demand instances are unconstrained")} else {return(res)}
} # end EM function
