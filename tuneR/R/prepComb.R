prepComb <- function(object, zero = 0, where = c("both", "start", "end")){
    if(!(is(object, "Wave") || is(object, "WaveMC"))) 
        stop("'object' needs to be of class 'Wave' or 'WaveMC'")
    validObject(object)
    where <- match.arg(where)
    starts <- 1
    ends <- length(object)  
    samples <- if(is(object, "Wave")) object@left else object@.Data[,1]  
    if(is.na(zero)) zero <- mean(samples)
    if(!is.numeric(zero) || length(zero) != 1)
        stop("zero must be either NA or a numeric value of length 1")
    up <- which(diff(sign(samples - zero)) > 0)
    if(length(up) > 1){
        if(where %in% c("both", "start"))
            starts <- up[1] + 1
        if(where %in% c("both", "end"))
            ends <- rev(up)[1]
    }else warning("returned object is unchanged\nat least two zero level crossings from negative to positive are required")
    return(object[starts:ends])
}
