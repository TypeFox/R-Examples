data(senc)

out <- ei.reg.bayes(cbind(dem, rep, non) ~ cbind(black, white, natam), 
              data = senc)

lambda <- lambda.reg.bayes(out, c("dem", "rep"), ret.mcmc = TRUE)
density.plot(lambda)

