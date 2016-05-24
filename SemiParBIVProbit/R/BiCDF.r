BiCDF <- function (u1, u2, family, par1, par2 = NULL, test = TRUE){

  epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)
  max.p   <- 0.9999999

if(test == TRUE){

    if(is.null(u1) == TRUE || is.null(u2) == TRUE) stop("The margins are not set or have length zero.")
    if(any(u1 > 1) || any(u1 < 0)) stop("First margin must be in [0,1].")
    if(any(u2 > 1) || any(u2 < 0)) stop("Second margin must be [0,1].")
    if(length(u1) != length(u2)) stop("The lengths of two margins do not match.")
    if( length(par1) < 1 ) stop("The dependence parameter must be set.")
    if( length(par1) > 1 && length(par1) != length(u2) ) stop("The dependence parameter has length different from the number of observations.")
    
    if(  family == 1  && any( abs(par1) >= 1) )  stop("The parameter of Gaussian must be in (-1,1).")
    if(  family == 55 && any( abs(par1) >= 1) )  stop("The parameter of AMH must be in (-1,1).")
    if(  family == 56 && any( abs(par1) >= 1) )  stop("The parameter of FGM must be in (-1,1).")
    if( (family == 2 || family == 4) && any(par1 <= 0) )  stop("The parameter of Clayton0/180 must be positive.")
    if( (family == 3 || family == 5) && any(par1 >= 0) )  stop("The parameter of Clayton90/270 must be negative.")
    if( (family == 10 || family == 12) && any(par1 < 1))  stop("The parameter of Gumbel0/180 must be in [1,oo).")
    if( (family == 11 || family == 13) && any(par1 > -1)) stop("The parameter of Gumbel90/270 must be in (-oo,-1].")
    if(  family == 14 && any(par1 == 0))                  stop("The parameter of Frank must be different from 0.") 
    if( (family == 6 || family == 8) && any(par1 <= 1))   stop("The parameter of Joe0/180 must be in (1,oo).")
    if( (family == 7 || family == 9) && any(par1 >= -1))  stop("The parameter of Joe90/270 must be in (-oo,-1).")
       
}       
       
    if(family == 1)                     res <- pbinorm( qnorm(u1), qnorm(u2), cov12 = par1)
    if(family %in% c(2,6,10,14,55,56))  res <- BCDF(u1, u2, family, par1)                            # 0
    if(family %in% c(4,8,12))           res <- u1 + u2 - 1 + BCDF(1 - u1, 1 - u2, family, par1)      # 180
    if(family %in% c(3,7,11))           res <- u2 - BCDF(1 - u1, u2, family, -par1)                  # 90
    if(family %in% c(5,9,13))           res <- u1 - BCDF(u1, 1 - u2, family, -par1)                  # 270

    res <- ifelse(res > max.p, max.p, res) 
    res <- ifelse(res < epsilon, epsilon, res) 

    res
      
}


