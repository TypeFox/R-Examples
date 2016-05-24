
library(RcppXts)

X  <- xts(1:4, order.by=Sys.time()+0:3)
X2 <- xts(1:4, order.by=Sys.time()+4:7)


stopifnot( xtsIsOrdered(X) )
stopifnot( xtsCoredata(X) == coredata(X) )
stopifnot( xtsIs(X) )
stopifnot( all.equal(coredata(xtsTry(as.zoo(X))), coredata(X) ) )
stopifnot( all.equal(index(xtsTry(as.zoo(X))), index(X) ) )
xtsRbind(X, X2, FALSE)
xtsRbind(X, X, TRUE)
xtsNaCheck(X)
xtsLag(X, 2L, TRUE)
