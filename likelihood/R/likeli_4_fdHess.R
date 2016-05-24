###################################################
# likeli_4_fdHess
#
# Allows us to use likeli as a function in fdHess in the nlme package.
# This can handle ragged arrays in the parameters being analyzed.
# par_2_analyze is the numeric vector corresponding to par in anneal.
# It will be passed to likeli_4_fdHess by fdHess.  It is a "flattened"
# version of par; that is, if par[[1]] was 3 elements long, then:
# par_2_analyze[1] = par[[1]][1]
# par_2_analyze[2] = par[[1]][2]
# par_2_analyze[3] = par[[1]][3]
# par_2_analyze[4] = par[[2]][1] etc.
# All other elements are as they are in anneal.  This will transfer the
# values in par_2_analyze to par for transfer to likeli.
#
# Author:  Lora Murphy, Cary Institute of Ecosystem Studies
# murphyl@caryinstitute.org
##################################################
likeli_4_fdHess<-function(par_2_analyze, model, par, var, source_data, pdf) {
  k = 1
  for (i in 1:length(par)) {
    for (j in 1:length(par[[i]])) {
      par[[i]][j] = par_2_analyze[k]
      k = k + 1
    }
  }
  likeli(model,par,var,source_data,pdf)
}