data(senc)

out <- ei.reg(cbind(dem, rep, non) ~ cbind(black, white, natam), 
              data = senc)
summary(out)


