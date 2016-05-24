### defines the inverse function of digamma called igamma for simplicity

## an extensive grid of x-values
.xg <- sort(c(10^(-70:-1),qexp(unique(pmin(seq(0,1,length=5e3)+1e-10,1-1e-10))),qcauchy(seq(0.999,1-1e-10,length=5e3))))
.dxg <- digamma(.xg)
igamma <- approxfun(.dxg,.xg)
rm(.xg,.dxg)

