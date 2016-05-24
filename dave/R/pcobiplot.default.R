pcobiplot.default <-
function(veg,method,y=1,...) {
    o.pcobiplot<- pcocoor(veg,method,y)
    class(o.pcobiplot)<- "pcobiplot"
    o.pcobiplot$call<- match.call()
    cat("Call:\n")
    print(o.pcobiplot$call) 
    o.pcobiplot
  }
