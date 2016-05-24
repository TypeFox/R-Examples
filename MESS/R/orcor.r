orcor <- function(object, type=c("default", "youden"), lambda=1) {

    type <- match.arg(type)

    rho <- object
    l <- 2*asin(rho)/pi
    
    or <- switch(type,
                 default=((1+rho)/(1-rho))^(1/lambda),
                 youden=((l+1)/(1-l))^2)
    or
}
 
