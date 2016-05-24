oldpar <- par(no.readonly = TRUE)

data(datadyncouponbonds)

## Diebold/Li estimation
dl_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "dl", lambda = 1/3)

## 3d yield curve plot
plot(dl_res)

## Estimated parameters
plot(param(dl_res))
summary(param(dl_res))

## Estimate Nelson/Siegel model
ns_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", tauconstr = list(c(0.2, 7, 0.2)), optimtype = "allglobal")

## Estimated parameters
plot(param(ns_res))
summary(param(ns_res))

## Estimate Svensson model
sv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "sv",tauconstr = list(c(0.2,7,0.2,0.5)))

## Plot start parameter grid search for t=1

## Estimated parameters
plot(param(sv_res))
summary(param(sv_res))

## Estimate Adjusted Svensson model
asv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "asv",tauconstr = list(c(0.2,10,0.2)))

## Estimated parameters
plot(param(asv_res))
summary(param(asv_res))


## Factor contributions at t=1
par(mfrow=c(2,2))
fcontrib(param(dl_res), index = 1, method="dl")
fcontrib(param(ns_res), index = 1, method="ns")
fcontrib(param(sv_res), index = 1, method="sv")
fcontrib(param(asv_res), index = 1, method="asv")



## Compare GOF
allgof <- cbind(summary(dl_res)$gof, summary(ns_res)$gof, summary(sv_res)$gof, summary(asv_res)$gof)
colnames(allgof) <- c("Diebold/Li", "Nelson/Siegel", "Svensson", "Adj. Svensson")
print(allgof)

par(oldpar)
