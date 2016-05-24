### Example of Section 3.4 plot (b)
set.seed(1234)
Q <- MixSim(BarOmega = 0.001, K = 6, p = 4)
pdplot(Pi = Q$Pi, Mu = Q$Mu, S = Q$S)