SNPtm.default <-
function(trange,tsl,x6,r6,...) {
     o.SNPtm<- SNPtm2(trange,tsl,x6,r6)
     o.SNPtm$call<- match.call()
     cat("Call:\n") 
     class(o.SNPtm) <- "SNPtm"
     print(o.SNPtm$call)
     cat("\n")
     cat("Time range:       0 -",o.SNPtm$n.time.steps,"\n")
     cat("Time step length:",o.SNPtm$time.step.length,"\n")
     cat("Species:       Initial state:   Final state:  Growth rate:")
     cat("\n")
     x<- o.SNPtm$sim.data
     for (i in 1:6){
          cat(sprintf(" %-20s %6.1f  %13.1f %13.3f",o.SNPtm$veg.types[7-i],x[1,i],x[o.SNPtm$n.time.steps,i],o.SNPtm$growth.rates[i]))
          cat("\n")
          }
     o.SNPtm
     }
