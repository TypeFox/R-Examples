data(senc)

out <- ei.reg(cbind(dem, rep, non) ~ cbind(black, white, natam), 
              data = senc)

lambda <- lambda.reg(out, c("dem", "rep"))
density.plot(lambda)

