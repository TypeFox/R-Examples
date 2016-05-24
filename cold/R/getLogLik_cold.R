setMethod("getLogLik",
    signature(object = "cold"),
    function (object) 
    {return(object@log.likelihood)})