
data(redistrict)
data(tuneA)
data(tuneB)
tuneB <- array(tuneB[[1]], dim = c(3, 2, 150))

out <- ei.MD.bayes(cbind(dem, rep, novote) ~ cbind(black, white, hispanic),
             data = redistrict, lambda1 = 4, lambda2 = 2,
             tune.list = list(tuneA, tuneB), sample = 1000,
             thin = 200, burnin = 0, verbose = 1000)

lambda <- lambda.MD(out, c("dem", "rep"))
density.plot(lambda)

