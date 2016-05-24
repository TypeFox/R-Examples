###################################################
# likeli_4_optim
#
# Allows us to use likeli as a function in optim in the stats package.
# This can handle ragged arrays in the parameters being analyzed.
# par_2_analyze is the numeric vector corresponding to par in anneal.
# It will be passed to likeli_4_optim by optim.  
# This does not support ragged arrays in par (but should it?)
#
# Author:  Lora Murphy, Cary Institute of Ecosystem Studies
# murphyl@caryinstitute.org
##################################################
likeli_4_optim<-function(par_2_analyze, model, par_names, var, source_data, pdf) {
  par <- as.list(par_2_analyze)
  names(par) <- par_names
  likeli(model,par,var,source_data,pdf)
}