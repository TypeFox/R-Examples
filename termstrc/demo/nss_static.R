oldpar <- par(no.readonly = TRUE)

data(govbonds)

ns_res <- estim_nss(govbonds, c("GERMANY", "AUSTRIA", "FRANCE"),matrange = c(0,30), method = "ns", tauconstr = list(c(0.2,5,0.1),c(0.2,5,0.1), c(0.2,5,0.1)))

print(ns_res)
plot(ns_res)
summary(ns_res)

## Plot startparameter grid search results
par(mfrow=c(1,3))
plot(ns_res$spsearch$GERMANY,main="GERMANY")
plot(ns_res$spsearch$AUSTRIA,main="AUSTRIA")
plot(ns_res$spsearch$FRANCE,main="FRANCE")

## Plot all yield curves in one figure
par(mfrow=c(1,1))
plot(ns_res,multiple=TRUE)

par(oldpar)

