
setMethod("skewness", signature(x = "Gumbel"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp))  
        return(skewness(as(x,"AbscontDistribution"),...))
    else{
         return( -12 * sqrt(6) * APERYCONSTANT / pi^3 ) 
# http://mathworld.wolfram.com/GumbelDistribution.html         
    }
})

