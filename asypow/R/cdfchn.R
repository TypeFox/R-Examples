cdfchn <- function(which, p=0.0, x=0.0, df=0.0, pnonc=0.0)
#**********************************************************************
#
#
#  Changed for version 1.2 to eliminate use of compiled fortran, and replace it with
#  standard R function.
#
#               Cumulative Distribution Function
#               Non-central Chi-Square
#
#
#                              Function
#
#
#     Calculates any one parameter of the non-central chi-square
#     distribution given values for the others.
#
#
#                              Arguments
#
#
#     WHICH --> Integer indicating which of the next three argument
#               values is to be calculated from the others.
#               Input range: 1..4
#               iwhich = 1 : Calculate P from X and DF
#               iwhich = 2 : Calculate X from P,DF and PNONC
#               iwhich = 3 : Calculate DF from P,X and PNONC
#               iwhich = 3 : Calculate PNONC from P,X and DF
#                    INTEGER WHICH
#
#     P <--> The integral from 0 to X of the non-central chi-square
#            distribution.
#            Input range: [0, 1-1E-16).
#                    DOUBLE PRECISION P
#
#     X <--> Upper limit of integration of the non-central
#            chi-square distribution.
#            Input range: [0, +infinity).
#            Search range: [0,1E300]
#                    DOUBLE PRECISION X
#
#     DF <--> Degrees of freedom of the non-central
#             chi-square distribution.
#             Input range: (0, +infinity).
#             Search range: [ 1E-300, 1E300]
#                    DOUBLE PRECISION DF
#
#     PNONC <--> Non-centrality parameter of the non-central
#                chi-square distribution.
#                Input range: [0, +infinity).
#                Search range: [0,1E4]
#                    DOUBLE PRECISION PNONC
#
#                              Method
#
#
#     Formula  26.4.25   of   Abramowitz   and   Stegun,  Handbook  of
#     Mathematical  Functions (1966) is used to compute the cumulative
#     distribution function.
#
#     Computation of other parameters involve a seach for a value that
#     produces  the desired  value  of P.   The search relies  on  the
#     monotinicity of P with the other parameter.
#
#
#                            WARNING
#
#     The computation time  required for this  routine is proportional
#     to the noncentrality  parameter  (PNONC).  Very large  values of
#     this parameter can consume immense  computer resources.  This is
#     why the search range is bounded by 10,000. 
#     When changing to internal R functions this limit is lowered to about
#     1417.  
#
#**********************************************************************
{
  request <- c("p","x","df","pnonc")
  stopifnot(which %in% 1:4)
  stopifnot( (p>=0.0), (p<=1.0-1E-16) )
  stopifnot( x >= 0.0 )
  stopifnot( df >= 0.0 )
  stopifnot( pnonc >= 0.0 )

# **** cdfchn should be loaded before the next statement ****
# well, unnecessary now!

#  ans <- .Fortran("cdfchn",
#		  as.integer(which),
#		  as.double(p),
#		  as.double(1.0-p),
#		  as.double(x),
#                  as.double(df),
#		  as.double(pnonc),
#		  status=as.integer(0),
#		  bound = as.double(0), 
#                  PACKAGE="asypow")
#  status <- ans$status
#  bound  <- ans$bound
#  if (status != 0) {
#    if (status == -1) stop("which should be 1,2,3 or 4")else 
#    if (status == -2) stop("p should be in range [0, 1-1E-16)") else 
#    if (status == -3) stop("x should be in range [0, +infinity)") else 
#    if (status == -4) stop("df should be in range [0, +infinity)") else 
#    if (status == -5) stop("pnonc should be in range [0, +infinity)") else 
#    if (status == 1) 
#	stop(paste(request[which],"lower than lowest search bound", bound)) else 
#    if (status == 2)
#	stop(paste(request[which],"larger than largest search bound",bound)) else
#        stop("SMALL, X, BIG not monotone in INVR")
#  } else {
#    if (which == 1) return(list(p=ans[[2]])) else 
#    if (which == 2) return(list(x=ans[[4]])) else 
#    if (which == 3) return(list(df=ans[[5]])) else 
#    if (which == 4) return(list(pnonc=ans[[6]]))
#  }
     if(which==1){ return( list(p=pchisq(x, df, pnonc))) } else {
       if(which==2){ return(list(x= qchisq(p, df, pnonc) ))}  else {
         if(which==3) { #calculate df:
              temp <- uniroot( function(df) pchisq(x, df, pnonc)-p, 
                               interval=c(0, Inf) )$root
              return( list(df=temp) ) } else {
         if(which==4) { #calculate pnonc: 
              temp <- uniroot( function(pnonc) pchisq(x, df, pnonc)-p, 
                               interval=c(0, 1417) )$root
              return( list(pnonc=temp) )   }}}}
}









