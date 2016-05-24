iprop <-
function (
                  	irates,
                  	ci.fun=NULL,
                  	ci.level=NULL )
{

  object = irates
  covar.code = c(object$covar.code, object$full.sample.code)
  covar.lab = c(object$covar.lab, object$full.sample.lab)
  n.cov = length(covar.code)

  event.code <- object$event.code
  event.lab = object$event.lab
  n.event <- length(event.code)
  	
  ## Calculation of z
  if ( is.null(ci.level) ) ci.level <- object$ci.level
  if ( is.null(ci.fun) ) ci.fun <- "score"
  else if(!(ci.fun %in% c("lin", "score"))) stop(paste("ci.fun", ci.fun, "is not implemented in iprop"))

  z <- qnorm(1-ci.level/2)

  empty.data <- as.data.frame(matrix(NA, nrow=n.cov, ncol=n.event))
  names(empty.data) <- event.code
  row.names(empty.data) <- covar.code
 
  ip <- var <- lower <- upper <- empty.data
 
  r = object$n[covar.code,]
  n = as.vector(rep(object$N, n.event)[covar.code])
  if(any(covar.code == object$full.sample.code)) n[which(covar.code == object$full.sample.code)] = sum(object$N)
  n = as.vector(n)

  ip = as.data.frame(r[covar.code, event.code] / n, row.names = covar.code)
  names(ip) = event.code

  var = ip*(1-ip)
 
      if ( ci.fun == "lin" ){
        lower <- ip-z*sqrt(var/n)
        upper <- ip+z*sqrt(var/n)
      }
      if ( ci.fun == "score" ){
        lower <- (2*r[covar.code,event.code]+z^2-z*sqrt(z^2+4*r[covar.code,event.code]*(1-ip)))/(2*(n+z^2))
        upper <- (2*r[covar.code,event.code]+z^2+z*sqrt(z^2+4*r[covar.code,event.code]*(1-ip)))/(2*(n+z^2))
        }
  
  res <- list (
               ip=ip,
               var=var,
               conf.lower=lower,
               conf.upper=upper,
               event.code = event.code,
               event.lab = event.lab,
               covar.code = object$covar.code,
               covar.lab = object$covar.lab,
               full.sample.code = object$full.sample.code,
               full.sample.lab = object$full.sample.lab,
               ci.level=ci.level,
               ci.fun=ci.fun )

  class(res) <- "iprop"
               
  return (res)
}

