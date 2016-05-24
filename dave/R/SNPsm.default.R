SNPsm.default <-
function(trange,tsl,diff,r6,...) {
     o.SNPsm<- SNPsm2(trange,tsl,diff,r6)
     o.SNPsm$call<- match.call()
     cat("Call:\n") 
     class(o.SNPsm) <- "SNPsm"
     print(o.SNPsm$call)
     o.SNPsm
     }
