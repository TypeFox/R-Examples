
#classes and methods in the following
#Method for coercing x of class relimplm to list 
#with as(x,"list")
setAs("relimplm","list",as.list.relimplm)
setMethod("show",signature(object="relimplm"),function(object) print.relimplm(object))
setMethod("show",signature(object="relimplmbooteval"),function(object) print.relimplmbooteval(object))
setMethod("show",signature(object="relimplmbootMI"),function(object) print.relimplmbooteval(object))
setMethod("show",signature(object="relimplmboot"),function(object) 
    {
    cat("Objects of class relimplmboot should not be printed.", "\n", 
        "To see their structure, use the function str().", "\n")
    })
# all subsequent methods done as S3-methods only
#setMethod("print",signature(x="relimplm"),function(x) print.relimplm(x))
#setMethod("print",signature(x="relimplmbooteval"),function(x) print.relimplmbooteval(x))
#setMethod("print",signature(x="relimplmbootMI"),function(x) print.relimplmbootMI(x))
#setMethod("print",signature(x="relimplmboot"),function(x) 
#    {
#    cat("Objects of class relimplmboot should not be printed.", "\n", 
#        "To see their structure, use the function str().", "\n")
#    })

#setMethod("summary",signature(object="relimplmbootMI"),function(object, ...) summary.relimplmbootMI(object, ...))
#setMethod("plot",signature(x="relimplm",y="missing"),function(x,y,...) plot.relimplm(x,...))
#setMethod("plot",signature(x="relimplmbooteval",y="missing"),function(x,y,...) plot.relimplmbooteval(x,...))
#setMethod("plot",signature(x="relimplmbootMI",y="missing"),function(x,y,...) plot.relimplmbootMI(x,...))
