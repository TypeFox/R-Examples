"print.sum.psgc" <-
function(x,...) {
cat("\n #### MCMC details         ####","\n")
cat("\n number of saved samples:",x$nsamp,"\n")
cat("\n average effective sample size:",mean(x$ESS),"\n")
cat("\n effective sample sizes","\n" )
print(round(x$ESS,2))

cat("\n #### Parameter estimation ####","\n")
cat(" Posterior quantiles of correlation coefficients:","\n")
print(round(x$QC,2))

cat("\n Posterior quantiles of regression coefficients:","\n")
print(round(x$QR,2))


                                }

