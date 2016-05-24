print.srank <-
function(x,...) {
    o.srank<- x
    method<- o.srank$method
    nspecies<- length(o.srank$species)
    if (method == "indval")  {
      cat("    Rank    No. Species                             Indval   error p","\n")
      for (i in 1:nspecies) {
        cat(sprintf("%6.0f  %6.0f  %-30s  %10.3f   %-10.3g",o.srank$rank[i],o.srank$species.no[i],o.srank$species[i],o.srank$Indval[i],o.srank$error.probability[i]))
        cat("\n")
      }
    }
    if (method == "jancey") {
      cat("    Rank    No. Species                            F-value     error p","\n")
      for (i in 1:nspecies) {
        cat(sprintf("%6.0f  %6.0f  %-30s  %10.3f   %-10.3g",o.srank$rank[i],o.srank$species.no[i],o.srank$species[i],o.srank$F_value[i],o.srank$error.probability[i]))
        cat("\n")
      }
    }
 }
