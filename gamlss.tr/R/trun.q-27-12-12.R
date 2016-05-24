trun.q <- function(par, family = "NO", type = c("left", "right", "both"), varying = FALSE, ...)
  {
   type <- match.arg(type)
  fname <- family
distype <- eval(call(family))$type
  if (mode(family) != "character" && mode(family) != "name")
  fname <- as.character(substitute(family))
   qfun <- paste("q",fname,sep="")
   pfun <- paste("p",fname,sep="")
 invcdf <- eval(parse(text=qfun))
    cdf <- eval(parse(text=pfun))
 if (!varying)
 {
   if (type=="both" && length(par)!= 2)  stop(paste("the length of par should be 2 \n"))
   if (type!="both" && length(par)!= 1)  stop(paste("the length of par should be 1 \n"))
#-- 	    
fun <- if (type=="left")  
       function(p, lower.tail = TRUE, log.p = FALSE, ...)
      {
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
      pp <- cdf(par,...)
      cc <- invcdf((pp+p*(1-pp)),...)
      cc
      }
     else if (type=="right")
       function(p, lower.tail = TRUE, log.p = FALSE, ...)
      {
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
    pp <- if (distype == "Discrete") cdf(par-1,...)
          else cdf(par,...) # added Friday, February 26, 2010              
      cc <- invcdf(p*pp,...)
      cc
      }
     else if (type=="both")    
      function(p, lower.tail = TRUE, log.p = FALSE, ...)
      {
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
      pp1 <- cdf(par[1],...)
      pp2 <- if (distype == "Discrete") cdf(par[2]-1,...)
             else cdf(par[2],...)
      cc <- invcdf((p*(pp2-pp1)+pp1),...)
      cc
      }
  }
     else # this is for varying truncation only 
  {
  if (type=="both" && dim(par)[2]!= 2)  stop(paste("the rows of par should be 2 \n"))
#--
 fun <- if (type=="left")  
       function(p, lower.tail = TRUE, log.p = FALSE, ...)
      {
    if (length(p)!=length(par)) 
          stop(paste("length of p must be equal to ", length(par), "\n","" ))     	
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
      pp <- cdf(par,...)
      cc <- invcdf((pp+p*(1-pp)),...)
      cc
      }
     else if (type=="right")
       function(p, lower.tail = TRUE, log.p = FALSE, ...)
      {
    if (length(p)!=length(par)) 
          stop(paste("length of p must be equal to ", length(par), "\n" )) 	
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
    pp <- if (distype == "Discrete") cdf(par-1,...)
          else cdf(par,...) # added Friday, February 26, 2010              
      cc <- invcdf(p*pp,...)
      cc
      }
     else if (type=="both")    
      function(p, lower.tail = TRUE, log.p = FALSE, ...)
      {
    if (dim(par)[1]!= length(p))  
          stop(paste("the length of p should be the equal to ", dim(par)[1],  "\n"))  	
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
      pp1 <- cdf(par[,1],...)
      pp2 <- if (distype == "Discrete") cdf(par[,2]-1,...)
             else cdf(par[,2],...)
      cc <- invcdf((p*(pp2-pp1)+pp1),...)
      cc
      }	
  }	
  fun
 }
