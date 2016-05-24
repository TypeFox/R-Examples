EMMax.sd <-
function(y,z,mu, sdtype, sd.fixed=0.05) {
  switch(sdtype,
    sd.const = EMMax.sd.const(y,z,mu),   # estimate constant sigma
    sd.free  = EMMax.sd.free(y,z,mu),    # estimate free sigma's
    sd.fixed = rep(sd.fixed,length(mu)) )
}
