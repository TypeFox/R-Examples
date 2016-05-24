setMethod("getAIC",
    signature(object = "bild"),
    function (object) 
    {
        return(object@aic)
    }
)

