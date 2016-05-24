###################################################################################
#kurtosis  --- code due to G. Jay Kerns, gkerns@ysu.edu
###################################################################################



setMethod("kurtosis", signature(x = "Gumbel"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(kurtosis(as(x,"AbscontDistribution"),...))
    else{
         return(12/5)
# http://mathworld.wolfram.com/GumbelDistribution.html         
    }
})

### source http://en.wikipedia.org/wiki/Generalized_extreme_value_distribution
###        http://en.wikipedia.org/wiki/Gumbel_distribution
###        http://en.wikipedia.org/wiki/Riemann_zeta_function 
