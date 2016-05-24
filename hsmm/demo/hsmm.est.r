# Executing the EM algorithm:
fit    <- hsmm(sim$obs, od = "norm", rd = "log", 
               pi.par = pipar, tpm.par = tpmpar, 
               od.par = odpar, rd.par = rdpar)

# The log-likelihood:
fit$logl

# Ehe estimated parameters:
fit$para

# For comparison, the estimated parameters seperately together with the true parameter values
# ar given below.
# Transition probability matrix:
tpmpar
fit$para$tpm
# Observation distribution:
odpar
fit$para$od
# Runlength distribution:
rdpar
fit$para$rd

