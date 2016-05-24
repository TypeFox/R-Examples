setMethod("as.array",
    signature(x = "magpie"),
    function (x) 
    {      
      return(as(x,"array"))
    }
)