sandwich.error <- function(g){
#  if (g$method!="nls")
#    stop("sandwich.error is currently only available for method='nls.'")
#  sandwich.nls.polar=coeftest(g$fit,vcov = sandwich)
#  sandwich.est.nls.polar=c(sandwich.nls.polar[,1])
#  sandwich.se.nls.polar=c(sandwich.nls.polar[,2])
#  tcrit.nls.polar <- qt(0.975,g$fit.statistics["d.f."])
#  sandwich.low.nls.polar = sandwich.est.nls.polar -tcrit.nls.polar*sandwich.se.nls.polar
#  sandwich.high.nls.polar = sandwich.est.nls.polar +tcrit.nls.polar*sandwich.se.nls.polar
#return(list("low"=sandwich.low.nls.polar,"high"=sandwich.high.nls.polar))
NULL
}
