setGeneric("lineX12",
    function(x, param1,param2,...) { standardGeneric("lineX12")} )

setMethod(
    f='lineX12',
    signature=signature(x = "x12Output"),definition=function(x,param1,param2,...)
    {
      plot(1)
    })
setMethod(
    f='lineX12',
    signature=signature(x = "x12Single"),definition=function(x,param1,param2,...)
    {
      lineX12(x@x12Output,param1=param1,param2=param2,...)
    })


#x <- s1@x12Output
setGeneric("specPlotX12",
    function(x, param1,param2,...) { standardGeneric("specPlotX12")} )

setGeneric("acfPlotX12",
    function(x, param1,param2,...) { standardGeneric("acfPlotX12")} )

setGeneric("seasFacX12",
    function(x, param1,param2,...) { standardGeneric("seasFacX12")} )
