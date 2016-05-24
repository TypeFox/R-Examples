pcaser.default <-
function(veg,plotlabels,y,...) {
     o.pcaser<- pcaser2(veg,plotlabels,y)
     o.pcaser$call<- match.call()
     cat("Call:\n") 
     class(o.pcaser) <- "pcaser"
     print(o.pcaser$call)
     o.pcaser
     cat("\n")
     cat("Eigenvalues (first 5):",round(o.pcaser$Eigv[1:5],digits=3),"\n") 
     o.pcaser
     }
