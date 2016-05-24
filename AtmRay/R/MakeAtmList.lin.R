MakeAtmList.lin = function(z0 = 0, c0 = 343, gc = 0, wx0 = 0, gwx = 0, wy0 = 0, gwy = 0, rho0 = 1.2929 * exp(-z0/6800), grho = -0.0001901323 * exp(-z0/6800)){
  # it makes the most sense for z0 to be a scalar, but it won't crash if it's a vector.  Anything else can be a vector or scalar.  Makes a list of ATMs.
  # replace gc = 0 with gc = epsilon
  gc[gc == 0] = 10^-9

  # make every possible combination out of the input variables.  
  out = meshgridn(list(z0, c0, gc, wx0, gwx, wy0, gwy, rho0, grho))

  ATMlist = list()
  for(i in 1:length(out[[1]])){
    ATMlist[[i]] = list(z0 = out[[1]][i], c0 = out[[2]][i], gc = out[[3]][i], wx0 = out[[4]][i], gwx = out[[5]][i], wy0 = out[[6]][i], gwy = out[[7]][i], rho0 = out[[8]][i], grho = out[[9]][i])
  }

  return(ATMlist)
}

