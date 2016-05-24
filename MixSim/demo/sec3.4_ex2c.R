### Example of Section 3.4 plot (c)
set.seed(1232)
Q <- MixSim(BarOmega = 0.05, K = 6, p = 4)
pdplot(Pi = Q$Pi, Mu = Q$Mu, S = Q$S)

