.charToDistFunc = function(distribution, type = "q") {                                                           ####   .CHARTODISTFUNC-FUNCTION
    fun = NULL
    if (identical("beta", distribution)) 
        fun = eval(parse(text = paste(type, "beta", sep = "")))
    if (identical("cauchy", distribution)) 
        fun = eval(parse(text = paste(type, "cauchy", sep = "")))
    if (identical("chi-squared", distribution)) 
        fun = eval(parse(text = paste(type, "chisq", sep = "")))
    if (identical("exponential", distribution)) 
        fun = eval(parse(text = paste(type, "exp", sep = "")))
    if (identical("f", distribution)) 
        fun = eval(parse(text = paste(type, "f", sep = "")))
    if (identical("geometric", distribution)) 
        fun = eval(parse(text = paste(type, "geom", sep = "")))
    if (identical("log-normal", distribution) || identical("lognormal", distribution))         ####
        fun = eval(parse(text = paste(type, "lnorm", sep = "")))
    if (identical("log-normal3", distribution) || identical("lognormal3", distribution))       ####
        fun = eval(parse(text = paste(type, "lnorm3", sep = "")))                              ####
    if (identical("logistic", distribution)) 
        fun = eval(parse(text = paste(type, "logis", sep = "")))
    if (identical("negative binomial", distribution)) 
        fun = eval(parse(text = paste(type, "nbinom", sep = "")))
    if (identical("normal", distribution)) 
        fun = eval(parse(text = paste(type, "norm", sep = "")))
    if (identical("poisson", distribution)) 
        fun = eval(parse(text = paste(type, "pois", sep = "")))
    if (identical("t", distribution)) 
        fun = eval(parse(text = paste(type, "t", sep = "")))
    if (identical("weibull", distribution)) 
        fun = eval(parse(text = paste(type, "weibull", sep = "")))
    if (identical("weibull3", distribution))                                                   ####
        fun = eval(parse(text = paste(type, "weibull3", sep = "")))                            ####
    if (identical("gamma", distribution)) 
        fun = eval(parse(text = paste(type, "gamma", sep = "")))
    if (identical("gamma3", distribution)) 
        fun = eval(parse(text = paste(type, "gamma3", sep = "")))
    return(fun)
}
.lfrm = function(wholeList, filterList) {
    if (!is.list(wholeList)) 
        stop(paste(deparse(substitute(wholeList)), "is not a list!"))
    if (length(wholeList) == 0) 
        return(wholeList)
    if (!is.list(filterList)) 
        stop(paste(deparse(substitute(filterList)), "is not a list!"))
    if (length(filterList) == 0) 
        return(wholeList)
    logVec = lapply(names(wholeList), "%in%", names(filterList))
    filteredList = wholeList[!unlist(logVec)]
    return(filteredList)
}
.lfkp = function(wholeList, filterList) {
    if (!is.list(wholeList)) 
        stop(paste(deparse(substitute(wholeList)), "is not a list!"))
    if (length(wholeList) == 0) 
        return(wholeList)
    if (!is.list(filterList)) 
        stop(paste(deparse(substitute(filterList)), "is not a list!"))
    if (length(filterList) == 0) 
        return(filterList)
    logVec = lapply(names(wholeList), "%in%", names(filterList))
    filteredList = wholeList[unlist(logVec)]
    return(filteredList)
} 
