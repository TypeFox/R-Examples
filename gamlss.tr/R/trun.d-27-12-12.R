# this is to create a truncating distribution with varying truncation according to  par
# not finished yet 
trun.d <-function(par, family = "NO", type = c("left", "right", "both"), varying = FALSE, ...)
  {
     type <- match.arg(type)
    fname <- family
    if (mode(family) != "character" && mode(family) != "name")
    fname <- as.character(substitute(family))
  distype <- eval(call(family))$type
     dfun <- paste("d",fname,sep="")
     pfun <- paste("p",fname,sep="")
      pdf <- eval(parse(text=dfun))
      cdf <- eval(parse(text=pfun))  
 if (!varying)
 {
  if (type=="both" && length(par)!= 2)  stop(paste("the length of par should be 2 \n"))
  if (type!="both" && length(par)!= 1)  stop(paste("the length of par should be 1 \n"))
  #--
  fun <- if (type=="left")  
    function(x, log = FALSE, ...)
    {
      if (distype=="Discrete" &&  any(x <= par))  
        stop(paste("x must be greater than ", par, "\n", ""))
      if (distype!="Discrete" && any(x < par))
        stop(paste("x must be greater or equal than ", par, "\n", ""))
      dfun <- pdf(x,log = TRUE,...)-log(1-cdf(par,...))
      dfun <- if (log == TRUE) dfun else exp(dfun)
      dfun
    }
  else if (type=="right")
    function(x, log = FALSE, ...)
    {
      if (distype=="Discrete" &&  any(x >= par))  
        stop(paste("x must be less than ", par, "\n", ""))
      if (distype!="Discrete" && any(x > par))
        stop(paste("x must be less or equal than ", par, "\n", ""))   
      dfun <-  if (distype=="Discrete") pdf(x, log = TRUE,...)-log(cdf(par-1,...))
               else pdf(x, log = TRUE,...)-log(cdf(par,...))
      dfun <- if (log == TRUE) dfun else exp(dfun)
      dfun
    } 
  else if (type=="both")    
    function(x, log = FALSE, ...)
    {
      if (distype=="Discrete" &&  (any(x <= par[1]) || any(x >= par[2])) )  
        stop(paste("x must be greater than", par[1], "and less than", par[2], 
                   "\n", ""))
      if (distype!="Discrete" && (any(x < par[1]) || any(x > par[2])) )
        stop(paste("x must be greater than", par[1], "and less or equal to", par[2], 
                   "\n", ""))  
      dfun <- if (distype=="Discrete") pdf(x, log = TRUE,...) - log(cdf(par[2]-1,...)-cdf(par[1],...)) 
              else pdf(x, log = TRUE,...) - log(cdf(par[2],...)-cdf(par[1],...))      
      dfun <- if (log == TRUE) dfun else exp(dfun)
      dfun
    }   
  }
  else # this is for varying truncation only 
  {
  #  browser()
    if (type=="both" && dim(par)[2]!= 2)  stop(paste("the rows of par should be 2 \n"))
    #--
    fun <- if (type=="left")  
      function(x, log = FALSE, ...)
      {
        if (length(x)!=length(par)) 
          stop(paste("length of x must be equal to ", length(par), "\n","" ))
        if (distype=="Discrete" &&  any(x <= par))  
          stop(paste("x must be greater than the truncation parameters", "\n", ""))
        if (distype!="Discrete" && any(x < par))
          stop(paste("x must be greater or equal than ", par, "\n", ""))
        dfun <- pdf(x,log = TRUE,...)-log(1-cdf(par,...))
        dfun <- if (log == TRUE) dfun else exp(dfun)
        dfun
      }
    else if (type=="right")
      function(x, log = FALSE, ...)
      {
        if (length(x)!=length(par)) 
          stop(paste("length of x must be equal to ", length(par), "\n" ))
        if (distype=="Discrete" &&  any(x >= par))  
          stop(paste("x must be less than ", par, "\n", ""))
        if (distype!="Discrete" && any(x > par))
          stop(paste("x must be less or equal than ", par, "\n", ""))
        dfun <- if (distype=="Discrete") pdf(x, log = TRUE,...)-log(cdf(par-1,...))
                else pdf(x, log = TRUE,...)-log(cdf(par,...))
        dfun <- if (log == TRUE) dfun else exp(dfun)
        dfun
      } 
    else if (type=="both")    
      function(x, log = FALSE, ...)
      {
        if (dim(par)[1]!= length(x))  
          stop(paste("the length of x should be the equal to ", dim(par)[1],  "\n"))
        if (distype=="Discrete" &&  (any(x <= par[,1]) || any(x >= par[,2])) )  
          stop(paste("x must be between the lower and uppper of the par argument", 
                     "\n", ""))
        if (distype!="Discrete" && (any(x < par[,1]) || any(x > par[,2])) )
          stop(paste("x must be between the lower and uppper of the par argument",
                     "\n", "")) 
        dfun <-   if (distype=="Discrete") pdf(x, log = TRUE,...) - log(cdf(par[,2]-1,...)-cdf(par[,1],...)) 
                  else pdf(x, log = TRUE,...) - log(cdf(par[,2],...)-cdf(par[,1],...))      
        dfun <- if (log == TRUE) dfun else exp(dfun)
        dfun
      }  
  }
  fun    
  }
