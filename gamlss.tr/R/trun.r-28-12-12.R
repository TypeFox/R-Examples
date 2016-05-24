trun.r <- function(par, family = "NO", type = c("left", "right", "both"), varying=FALSE,...)
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
       function(n,...)
    {
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))   
     n <- ceiling(n)
     p <- runif(n)
    pp <- cdf(par,...)
     r <- invcdf((pp+p*(1-pp)),...)
     r
    }
     else if (type=="right")
     function(n,...)
    {
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))   
     n <- ceiling(n)
     p <- runif(n)
    pp <- cdf(par,...)
     r <- invcdf(p*pp,...)
     r
    }
     else if (type=="both")    
      function(n,...)
    {
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))   
     n <- ceiling(n)
     p <- runif(n)
   pp1 <- cdf(par[1],...)
   pp2 <- cdf(par[2],...)
     r <- invcdf(p*(pp2-pp1)+pp1,...)
     r
    }
  }
   else # this is for varying truncation only 
   {
    if (type=="both" && dim(par)[2]!= 2)  stop(paste("the rows of par should be 2 \n"))
#--
    fun <- if (type=="left")  
      function(n,...)
      {     
        if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))   
        n <- ceiling(n)
        if (n!=length(par)) 
          stop(paste(" n must be equal to ", length(par), "\n","" ))
        p <- runif(n)
        pp <- cdf(par,...)
        r <- invcdf((pp+p*(1-pp)),...)
        r
      }
    else if (type=="right")
      function(n,...)
      {
        if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))   
        n <- ceiling(n)
        if (n!=length(par)) 
          stop(paste(" n must be equal to ", length(par), "\n","" )) 
        p <- runif(n)
        pp <- cdf(par,...)
        r <- invcdf(p*pp,...)
        r
      }
    else if (type=="both")    
      function(n,...)
      {
        if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))   
        n <- ceiling(n)
        if (dim(par)[1]!= n)  
          stop(paste("the n should be the equal to ", dim(par)[1],  "\n"))  
        p <- runif(n)
        pp1 <- cdf(par[,1],...)
        pp2 <- cdf(par[,2],...)
        r <- invcdf(p*(pp2-pp1)+pp1,...)
        r     
       }
   }  
fun
}
