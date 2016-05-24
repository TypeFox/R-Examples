Mtabs.default <-
function(veg,method="raw",y.r=0.5,y.s=0.25,k.r=3,k.s=5,ndiffs=10,...) {
     o.Mtabs<- mtab(veg,method,y.r,y.s,k.r,k.s,ndiffs)
     o.Mtabs$call<- match.call()
     cat("Call:\n") 
     class(o.Mtabs) <- "Mtabs"
     print(o.Mtabs$call)
     o.Mtabs
     cat("\n")
     cat("CA eigenvalues, %:      ",round(o.Mtabs$CAeig.rel*100,digits=2)[1:3],"\n") 
     o.Mtabs
     if (method == "raw")   cat("AOC eigenvalues, %:     ",o.Mtabs$AOCeig.rel*100[1:3],"\n") 
     if (method == "sort")  cat("AOC eigenvalues, %:     ",o.Mtabs$AOCeig.rel*100[1:3],"\n") 
     if (method == "ca")    cat("AOC eigenvalues, %:     ",o.Mtabs$AOCeig.rel*100[1:3],"\n") 
     if (method == "clust") cat("AOC eigenvalues, %:     ",o.Mtabs$AOCeig.rel*100[1:3],"\n") 
     if (method == "aoc")   cat("AOC eigenvalues, %:     ",round(o.Mtabs$AOCeig.rel*100,digits=2)[1:3],"\n") 
     if (method == "mulva") cat("AOC eigenvalues, %:     ",round(o.Mtabs$AOCeig.rel*100,digits=2)[1:3],"\n") 
     o.Mtabs
     cat("\n")
     if (method == "raw")   cat("Mean square contingency coefficient:",o.Mtabs$MSCC,"\n") 
     if (method == "sort")  cat("Mean square contingency coefficient:",o.Mtabs$MSCC,"\n") 
     if (method == "ca")    cat("Mean square contingency coefficient:",o.Mtabs$MSCC,"\n") 
     if (method == "clust") cat("Mean square contingency coefficient:",o.Mtabs$MSCC,"\n") 
     if (method == "aoc")   cat("Mean square contingency coefficient:",round(o.Mtabs$MSCC,digits=5),"\n") 
     if (method == "mulva") cat("Mean square contingency coefficient:",round(o.Mtabs$MSCC,digits=5),"\n") 
     o.Mtabs
     }
