oldpar <- par(no.readonly = TRUE)

## Load yields matrix from csv
data(zyields)
x <- zyields

maturities <- c(1/12,3/12,6/12,9/12,1:12)
yields <- as.matrix(x[,2:ncol(x)])
dates <- as.Date(x[,1],format="%d.%m.%Y")

## Call class constructor
datazeroyields <- zeroyields(maturities, yields, dates)

## Estimate Diebold/Li model
dl_res <- estim_nss(datazeroyields, "dl", lambda = 1/2)
print(dl_res)
summary(dl_res)
plot(param(dl_res))

## Estimate Nelson/Siegel model
ns_res <- estim_nss(datazeroyields, "ns", tauconstr = c(0.2, 6, 0.1))
print(ns_res)
summary(ns_res)
plot(param(ns_res))

## Plot start parameter grid search for t=1
plot(ns_res$spsearch[[1]])

## Estimate Svensson model
sv_res <- estim_nss(datazeroyields, "sv", tauconstr =  c(0.2, 5, 0.1, 0.5))
print(sv_res)
summary(sv_res)
plot(param(sv_res))

## Plot start parameter search for t=1
plot(sv_res$spsearch[[1]])

## Plot yield curves in 3D
plot(sv_res)

## Estimate Adjusted Svensson model
asv_res <- estim_nss(datazeroyields, "asv",  tauconstr =  c(0.2, 7, 0.1))
plot(param(asv_res))

## Compare GOF
allgof <- cbind(summary(dl_res)$gof, summary(ns_res)$gof, summary(sv_res)$gof, summary(asv_res)$gof)
colnames(allgof) <- c("Diebold/Li", "Nelson/Siegel", "Svensson", "Adj. Svensson")
print(allgof)

par(oldpar)









