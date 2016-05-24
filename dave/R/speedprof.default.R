speedprof.default <-
function(veg,timescale,orders,y=1,adjust,...) {
     o.speedprof<- speedprof2(veg,timescale,orders,y=1,adjust)
     o.speedprof$call<- match.call()
     cat("Call:\n") 
     class(o.speedprof) <- "speedprof"
     print(o.speedprof$call)
     o.speedprof
     }
