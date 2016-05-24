### Test main functions
set.seed(1234)
repeat {
   Q <- MixSim(BarOmega = 0.01, MaxOmega = 0.05, K = 5, p = 7)
   if (Q$fail == 0) break
}

overlap(Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
pdplot(Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
