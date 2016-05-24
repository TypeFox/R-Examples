check.distribution.args <-
function (distribution, param.list, check.params = TRUE) 
{
    idist <- charmatch(distribution, .Distribution.abb, nomatch = 0)
    if (idist == 0) 
        stop(paste("Unknown or ambiguous distribution abbreviation. ", 
            "See the help file for 'EnvStats::Distribution.df' for", 
            "a list of distribution names."))
    dist.abb <- .Distribution.abb[idist]
    dist.name <- .Distribution.name[idist]
    if (dist.name == "Stable") 
        stop(paste("No method for Stable Distribution.  The only function", 
            "for the Stable Distribution is 'rstab', which produces", 
            "random numbers from this distribution."))
    dist.type <- EnvStats::Distribution.df[idist, "Type"]
    n.dist.params <- EnvStats::Distribution.df[idist, "Number.parameters"]
    dist.params.names <- unlist(EnvStats::Distribution.df[idist, 
        paste("Parameter", 1:n.dist.params, sep = ".")])
    if (check.params) {
        if (!is.list(param.list)) 
            stop("'param.list' must be a list.")
        if (any(dist.abb == c("beta", "chisq", "f", "t")) && 
            length(param.list) < n.dist.params) 
            param.list <- c(param.list, list(ncp = 0))
        names.param.list <- names(param.list)
        if (length(param.list) != n.dist.params || is.null(names.param.list) || 
            any(names.param.list == "")) 
            stop(paste("You must supply the names and values of all of the ", 
                "distribution parameters via the argument ", 
                "'param.list'.  Use the form\n\n\t", "param.list = list(name1=value1, name2 = value2)\n\n\t", 
                "for example, for a distribution with two parameters.  ", 
                "You may abbreviate distribution parameter names as ", 
                "long as the abbreviation uniquely identifies each ", 
                "distribution parameter relative to the other ", 
                "distribution parameters.  See the help file for ", 
                "'EnvStats::Distribution.df' for a list of distribution ", 
                "parameter names.", sep = ""))
        param.index <- pmatch(names.param.list, dist.params.names)
        if (sum(!is.na(param.index)) != n.dist.params) 
            stop(paste("Unknown or ambiguous argument name(s) for the", 
                "distribution parameter(s).  See the help file for", 
                "'EnvStats::Distribution.df' for", "a list of distribution parameter names."))
        if (!all(sapply(param.list, length) == 1)) 
            stop("All distribution parameters must be scalars.")
        if (any(sapply(param.list, is.na))) 
            stop("All distribution parameters must be non-missing.")
        param.list <- param.list[order(param.index)]
        names(param.list) <- dist.params.names
        for (i in 1:n.dist.params) assign(dist.params.names[i], 
            param.list[[i]])
        out.of.bounds.vec <- logical(n.dist.params)
        for (i in 1:n.dist.params) {
            param <- param.list[[i]]
            out.of.bounds.vec[i] <- (param < eval(parse(text = EnvStats::Distribution.df[idist, 
                paste("Parameter", i, "Min", sep = ".")]))) || 
                (param > eval(parse(text = EnvStats::Distribution.df[idist, 
                  paste("Parameter", i, "Max", sep = ".")])))
        }
        if (any(out.of.bounds.vec)) 
            stop(paste("Illegal values for the following distribution parameter(s):\n\t\t", 
                paste(dist.params.names[out.of.bounds.vec], collapse = ", "), 
                "\n\t", "See the help file for 'EnvStats::Distribution.df' for ", 
                "a list of legal values for distribution parameters.", 
                sep = ""))
    }
    ret.list <- list(dist.abb = dist.abb, dist.name = dist.name, 
        dist.type = dist.type, n.dist.params = n.dist.params, 
        dist.params.names = dist.params.names)
    if (check.params) 
        ret.list <- c(ret.list, list(param.list = param.list))
    ret.list
}
