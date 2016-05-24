# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
attribute.prob.bin<- function(x, ...){
# retreives data from a verify object.

# assign("obar.i", x$obar.i)
# assign("thres", x$thres)
# assign("prob.y", x$prob.y)
# assign("obar", x$obar)
# assign("class", "prob.bin")
# assign("obs", x$obs)
# assign("pred", x$pred)
# assign("bins", x$bins)
# assign("x", x$y.i)

# do.call("attribute.default", list(x, obar.i, prob.y, obar, class, obs=obs, pred = pred, thres = thres, bins = bins,...))

res <- attribute.default(x$y.i, obar.i=x$obar.i, prob.y=x$prob.y, obar=x$obar,
    class="prob.bin", obs=x$obs, pred=x$pred, thres=x$thres, bins=x$bins, ...)
return(res)
}
