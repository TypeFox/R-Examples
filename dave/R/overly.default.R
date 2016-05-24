overly.default <-
function(veg,Plot.no,y,sint,...) {
     o.overly<- overly2(veg,Plot.no,y,sint)
     o.overly$call<- match.call()
     cat("Call:\n") 
     class(o.overly) <- "overly"
     print(o.overly$call)
     o.overly
     nt<- o.overly$n.tsteps
     cat("Number of time steps in new time series:  ",nt,"\n")
     cat("Time span of the new time series:          0 -",o.overly$tsteps[nt],"\n")
     o.overly
     }
