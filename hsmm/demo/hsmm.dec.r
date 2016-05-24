# Computation of the smoothing probabilities:
fit.sm <- hsmm.smooth(sim$obs, od = "norm", rd = "log", 
                      pi.par = pipar, tpm.par = tpmpar, 
                      od.par = odpar, rd.par = rdpar)

# Executing the Viterbi algorithm:
fit.vi <- hsmm.viterbi(sim$obs, od = "norm", rd = "log", 
                       pi.par = pipar, tpm.par = tpmpar, 
                       od.par = odpar, rd.par = rdpar)


# The first 15 smoohting probabilities:
round(fit.sm$sm[, 1:15], 2)
# The first 15 values of the resulting path:
fit.sm$path[1:15]

# The first 15 values of the resulting path (Viterbi):
fit.vi$path[1:15]

# For comparison, the real/simulated path (first 15 values):
sim$path[1:15]

