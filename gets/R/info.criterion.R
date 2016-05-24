info.criterion <-
function(logl, n = NULL, k = NULL,
  method = c("sc", "aic", "aicc", "hq"))
{
#check arguments:
if (!is.numeric(logl))
  stop(" Non-numeric argument to logl")
if(is.null(n))
  stop(" Value n is NULL")
if(is.null(k))
  stop(" Value k is NULL")

#initiate:
method <- match.arg(method)
#info.val <- NULL

#Schwarch criterion
if(method == "sc") info.val <- I(-2)*logl/n + k*log(n)/n

#Akaike criterion
if(method == "aic") info.val <- -2*logl/n + 2*k/n

#AICC of Hurvich and Tsai (1989):
if(method == "aicc") info.val <- -2*logl/n + 2*k/(n-k-1)

#Hannan-Quinn criterion
if(method == "hq") info.val <- -2*logl/n + 2*k*log(log(n))/n

#out values:
out <- list()
out$method <- method
out$n <- n
out$k <- k
out$value <- info.val

return(out)

}
