vvelocity.default <-
function(pveg,timescale,y,...) {
     o.vvelocity<- vvelocity2(pveg,timescale,y)
     o.vvelocity$call<- match.call()
     cat("Call:\n") 
     class(o.vvelocity) <- "vvelocity"
     print(o.vvelocity$call)
     o.vvelocity
     }
