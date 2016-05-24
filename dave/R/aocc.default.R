aocc.default <-
function(veg,o.rgr,o.sgr,...) {
     o.aocc<- aoc(veg,o.rgr,o.sgr)
     o.aocc$call<- match.call()
     cat("Call:\n") 
     class(o.aocc) <- "aocc"
     print(o.aocc$call)
     o.aocc
     cat("\n")
     cat("Mean square contingency coefficient:",round(o.aocc$MSCC,digits=5),"\n") 
     o.aocc
     cat("\n")
     cat("Chi squared:            ",round(o.aocc$chisquare,digits=3),"\n") 
     o.aocc
     cat("Eigenvalues:            ",round(o.aocc$eig,digits=3),"\n") 
     o.aocc
     cat("Canonical correlations: ",round(o.aocc$cancorr,digits=3),"\n") 
     o.aocc
     }
