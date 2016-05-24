# Simulation of sequences of observations and hidden states from a 
# 3-state HSMM with a logarithmic runlength distribution and a 
# conditional Gaussian distributions.

### Setting up the parameter values:
# Initial probabilities of the semi-Markov chain:
pipar  <- rep(1/3, 3)
# Transition probabilites:
# (Note: For two states, the matrix degenerates, taking 0 for the 
# diagonal and 1 for the off-diagonal elements.)
tpmpar <- matrix(c(0, 0.5, 0.5,
                   0.7, 0, 0.3,
                   0.8, 0.2, 0), 3, byrow = TRUE)
# Runlength distibution:
rdpar  <- list(p = c(0.98, 0.98, 0.99))
# Observation distribution:
odpar  <- list(mean = c(-1.5, 0, 1.5), var = c(0.5, 0.6, 0.8))

# Invoking the simulation:
sim    <- hsmm.sim(n = 2000, od = "norm", rd = "log", 
                   pi.par = pipar, tpm.par = tpmpar, 
                   rd.par = rdpar, od.par = odpar, seed = 3539)


# The first 15 simulated observations:
round(sim$obs[1:15], 3)

# The first 15 simulated states:
sim$path[1:15]
