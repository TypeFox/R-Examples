# Power polynomial
# This is an experimental non-linear additive function 
# created Monday, February 2, 2004 at 11:09
# Mikis Stasinopoulos
#----------------------------------------------------------------------------------------
pp <-function(x, start=list(), shift = NULL, scale = NULL ) 
{
# 1st function ------------------
fp.scale <-function(x)
   {
    if(min(x) <= 0) 
    {        z <- sort(x)[-1] - sort(x)[ - length(x)]
         shift <- min(z[z > 0]) - min(x)
    }
    else shift <- 0
    range <- max(x) - min(x)
    scale <- 10^(sign(log10(range)) * trunc(abs(log10(range))))
    list(shift=shift, scale=scale)
   }
#-end of 1st function------------
    scall <- deparse(sys.call())
    if(ncol(as.matrix(x)) > 1)
    stop(paste("The default smoother is bivariate; you gave a matrix as an argument in ",
            scall, "\n"))
    if(!is.null(levels(x))) 
        {
        if(inherits(x, "ordered"))
            x <- as.numeric(x)
        else stop("unordered factors cannot be used as smoothing variables")
        }
    npoly <- length(start)
    if(npoly==0) stop("no starting values for the power parameter(s) is set")     
    if(!any(npoly %in%c(1,2,3))) 
       stop("the number of power polynomials can be 1,2 or 3")
    if(is.null(scale)|is.null(shift)) 
        {
           out <- fp.scale(x)
         shift <- out$shift
         scale <- out$scale
        }
    if(shift > 0) 
       warning("The origing of the x variable has been shifted by ", shift, "\n" )
    x <- x + shift
    x <- x/scale 
    xvar <- rep(0, length(x))
    attr(xvar, "npoly") <- npoly
    attr(xvar, "start") <- start
    attr(xvar, "values") <- x
    real.call <- substitute(gamlss.pp(data[[scall]], z, w))
    attr(xvar, "call") <- real.call
    attr(xvar, "class") <- "smooth"
    xvar
}
# the definition of the backfitting additive function
gamlss.pp <-function(x, y, w)
{
     xvar <- as.vector(attr(x,"values"))
    npoly <- as.vector(attr(x, "npoly"))
    start <- as.list(attr(x,"start"))

if (npoly==1) 
    { 
    wy <-sqrt(w)*y
    fit <- nls(wy~cbind(1*sqrt(w),(xvar^p1)*sqrt(w)),  start=c(p1=start[1]), algorithm = "plinear")
     fv <- coef(fit)[2]+coef(fit)[3]*xvar^coef(fit)[1]
    }
 if (npoly==2)
    {
    wy <-sqrt(w)*y
    fit <- nls(wy~cbind(1*sqrt(w),(xvar^p1)*sqrt(w),(xvar^p2)*sqrt(w) ),  start=c(p1=start[1], p2=start[2]), algorithm = "plinear")
     fv <- coef(fit)[3]+coef(fit)[4]*xvar^coef(fit)[1]+coef(fit)[5]*xvar^coef(fit)[2] 
   } 
  if (npoly==3)
    {
    wy <-sqrt(w)*y
    fit <- nls(wy~cbind(1*sqrt(w),(xvar^p1)*sqrt(w),(xvar^p2)*sqrt(w),(xvar^p3)*sqrt(w)),  
              start=c(p1=start[1], p2=start[2], p3=start[3]), algorithm = "plinear")
     fv <- coef(fit)[4]+coef(fit)[5]*xvar^coef(fit)[1]+coef(fit)[5]*xvar^coef(fit)[2]+coef(fit)[7]*xvar^coef(fit)[3] 
   }  
    #new.start<- as.vector(coef(fit)[1])
 
   # eval(expression(attr(data[["pp(x, start = list(p1 = 0.5))"]],"start")),sys.parent())
   #  eval(expression(attr(mobject$smooth.frame[["pp(x, start = list(p1 = 0.5))"]],"start")),sys.parent(2))
   # eval(expression(attr(mobject$smooth.frame[["pp(x, start = list(p1 = 0.5))"]],"start")<-  eval(expression(new.start),sys.parent(2)) ),sys.parent(3))

   # attr(data[["pp(x, start = list(p1 = 0.5))"]],"start") <<- 3
   #  as.vector(coef(fit)[1])

residuals <- y-fv 
     if(npoly==1) df <- 2 
     if(npoly==2) df <- 4 
     if(npoly==3) df <- 6 
list(x = xvar, fitted.values=fv, residuals=residuals,
     nl.df =df, lambda=NA, 
     coefSmo=coef(fit), var=fv )    # var=fv has to fixed
}
      
      
      
      
  
