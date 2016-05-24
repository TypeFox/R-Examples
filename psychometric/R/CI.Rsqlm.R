"CI.Rsqlm" <-
function (obj, level=.95)
 {
 l <- level
 rsq <- summary(obj)$r.squared
 k <-  summary(obj)$df[1] - 1
 n <- obj$df + k + 1
 mat <- CI.Rsq (rsq, n, k, level=l)
 return(mat)
 }

