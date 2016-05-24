data(redistrict)

out <- ei.reg.bayes(cbind(dem, rep, novote) ~ cbind(black, white, hispanic), 
                    data = redistrict)
summary(out)


