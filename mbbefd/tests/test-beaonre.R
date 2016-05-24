library(mbbefd)

data(beaonre)
x <- beaonre$ClaimCost/beaonre$SumInsured


# #Nelder Mead
# mledist(x[x!=1], "gbeta", start=list(shape0=3.373523e-03, 
#     shape1=2.913619e+02, shape2=7.198632e+00), control=list(trace=1,REPORT=1))
# #L-BFGS-B
# mledist(x[x!=1], "gbeta", lower=0, start=list(shape0=3.373523e-03, 
#     shape1=2.913619e+02, shape2=7.198632e+00), control=list(trace=1, REPORT=1, fnscale=1e-6))
# 
# fitDR(x, "oigbeta", method="mle", control=list(trace=1, REPORT=1))
# fitDR(x, "MBBEFD", method="mle", control=list(trace=1, REPORT=1))


dlist <- c("oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD")
dlist <- c("oistpareto", "oibeta", "mbbefd", "MBBEFD")
flist <- lapply(dlist, function(d) {print(d);
  fitDR(x, d, method="mle")})
names(flist) <- dlist


cdfcomp(flist, do.points=FALSE, leg=dlist)
ppcomp(flist, leg=dlist, fitpch=".", addlegend = FALSE)
legend("bottomright", fill=c("red", "green", "blue", "cyan"), leg=dlist)

qqcomp(flist, leg=dlist, use.ppoints=TRUE)
