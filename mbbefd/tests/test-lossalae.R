library(mbbefd)

data(lossalaefull)
x <- lossalaefull$Loss/lossalaefull$Limit


# fitDR(x, "oigbeta", method="mle", control=list(trace=1, REPORT=1))

dlist <- c("oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD")
flist <- lapply(dlist, function(d) {print(d);
  fitDR(x, d, method="mle")})
names(flist) <- dlist


cdfcomp(flist, do.points=FALSE, leg=dlist, xlogscale = TRUE)
ppcomp(flist, leg=dlist, fitpch=".", addlegend = FALSE)
legend("bottomright", fill=c("red", "green", "blue", "cyan", "purple"), leg=dlist)

qqcomp(flist, leg=dlist, use.ppoints=TRUE)
