DeconSeis = function(GH, inst, L, fl = 0.1, fh = NaN, bitweight = NULL, dec = rep(1, length(GH$JSTR))){
  # GH is seismogram
  # inst is vector of indices in L corresponding to traces in GH (0 for no deconvolution)
  # L is list of DPZs of instruments used
  for(i in which(inst != 0)){
    # make sure that sample rate is correct
    if(GH$dt[i] != (L[[inst[i]]]$dt * dec[i])){
      stop('GH$dt must equal sample rates in L times dec')
    }
    GH$JSTR[[i]] = DeconTrace(GH$JSTR[[i]], L[[inst[i]]], fl = fl, fh = fh, bitweight = bitweight[i], dec = dec[i])
    GH$units[i] = 'm/s'
  }
  GH$process = c(GH$process, paste('DeconRec', fl, fh))
  invisible(GH)
}
