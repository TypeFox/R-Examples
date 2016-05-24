library(gnm)
set.seed(1)

RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               MultHomog(origin, destination), family = poisson,
               data = occupationalStatus)

print(RChomog$deviance, digits=10)
print(RChomog$df)


###  Fit an association model with homogeneous row-column effects
set.seed(4)
### Set diagonal elements to NA (rather than fitting exactly)
dat <- as.data.frame(friend)
id <- with(dat, r == c)
dat[id,] <- NA
rc2 <- gnm(Freq ~ r + c + instances(MultHomog(r, c), 2),
           family = poisson, data = dat, iterStart = 0)

print(rc2$deviance, digits=10)
print(rc2$df)
