logMargDenNR = function(ll) {
#
# purpose:  compute log marginal density using Newton-Raftery
#  importance sampling estimator: 1/ (1/g sum_g exp(-log like) )
#    where log like is the likelihood of the model evaluated as the 
#    posterior draws (x).
#
# arguments:
#  ll -- vector of log-likelihood values evaluated at posterior draws
#
# output:
#  estimated log-marginal density

med=median(ll)
return(med-log(mean(exp(-ll+med))))
}

