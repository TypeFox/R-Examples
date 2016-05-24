
setMethod("setNames",
    signature(object = "magpie"),
    function (object, nm) 
    {
        getNames(object) <- nm
        return(object)
    }
)
