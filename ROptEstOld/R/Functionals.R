
setMethod("var", signature(x = "Gumbel"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(var(as(x,"AbscontDistribution"),...))
    else{  b <- scale(x)
            return(b^2 * pi^2/6)
    }})
## http://mathworld.wolfram.com/GumbelDistribution.html

