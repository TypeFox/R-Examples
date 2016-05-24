pctab <- function( TT, margin=length( dim( TT ) ), dec=1 )
  {
  nd <- length( dim( TT ) )
  sw <- (1:nd)[-margin[1]]
  rt <-
  sweep( addmargins( TT,
                     margin,
                     list( list( All=sum,
                                   N=function( x ) sum( x )^2/100 ) ) ),
         sw,
         apply( TT,
                sw,
                sum )/100,
         "/" )
  if( dec ) print( round( rt, dec ) )
  invisible( rt )
  }
