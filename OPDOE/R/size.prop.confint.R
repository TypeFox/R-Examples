size.prop.confint <- function(p=NULL,delta,alpha)
{
 q=qnorm(1-alpha/2)
 if (is.null(p)) {n=q^2/(4*delta^2)}
 else {n=p*(1-p)*q^2/delta^2}
 n=ceiling(n)
 return(n)
}
 
