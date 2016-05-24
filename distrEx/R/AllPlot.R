# plot
setMethod("plot", signature(x="MultivariateDistribution",y="missing"), 
    function(x, y,...){ 
        warning("'plot' not yet implemented for objects",
                " of class ", class(x))
    })

setMethod("plot", signature(x="UnivariateCondDistribution",y="missing"), 
    function(x, y,...){ 
        warning("'plot' not yet implemented for objects",
                " of class ", class(x))
    })


