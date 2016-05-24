setMethod("getAIC",
    signature(object = "cold"),
    function (object) 
    {
        return(object@aic)
    }
)

