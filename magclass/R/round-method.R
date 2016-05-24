setMethod("round",
    signature(x = "magpie"),
    function (x, digits=0) 
    {      
        x@.Data <- round(x@.Data,digits=digits)
        return(x)
    }
)