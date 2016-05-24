setMethod("getLogLik",
    signature(object = "bild"),
    function (object) 
    {return(object@log.likelihood)})