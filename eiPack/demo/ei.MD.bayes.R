
data(redistrict) 
data(tuneA)
data(tuneB)
tuneB <- array(tuneB[[1]], dim = c(3, 2, 150))

out2 <- ei.MD.bayes(cbind(dem, rep, novote) ~ cbind(black, white, hispanic),
             data = redistrict, lambda1 = 4, lambda2 = 2,
                    tune.list = list(tuneA, tuneB), sample = 1000,
                    thin = 200, burnin = 0, verbose = 1000)
summary(out2)


data(senc)
out3 <- ei.MD.bayes(cbind(dem, rep, non) ~ cbind(white, black, natam),
                       covariate = ~ I(white/total),
                       data = senc,
                       sample = 1000,
                       thin = 100, burnin = 100000, verbose = 1000, ret.beta='r')
summary(out3)
