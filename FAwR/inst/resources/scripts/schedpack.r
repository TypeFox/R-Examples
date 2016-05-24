#	$Id: schedpack.r 4519 2010-12-19 23:32:13Z hamannj $	

## this is a set of functions to help with the R book
## for harvest scheduling.

##    No.   Row name   St   Activity     Lower bound   Upper bound    Marginal
## ------ ------------ -- ------------- ------------- ------------- -------------
##      1 hvol1        NS             0             0             =       78.3347
##      2 hvol2        NS             0             0             =      -14.3923
##      3 hvol3        NS             0             0             =      -104.408

get.row.report <- function( lp ) {
  
  ## get the number of columns
  nrows <- lpx_get_num_rows( lp )
  
  ## build the table for the rows
  r <- NULL
  for( i in 1:nrows ) {
    r <- rbind( r,
               as.data.frame( cbind( lpx_get_row_name( lp, i ),
                                    lpx_get_row_stat( lp, i ),
                                    ##glpk_strerror(lpx_get_row_stat(lp,i)),
                                    lpx_get_row_prim( lp, i ),
                                    lpx_get_row_lb( lp, i ),
                                    lpx_get_row_ub( lp, i ),
                                    lpx_get_row_dual( lp, i ),
                                    ##lpx_get_row_b_ind( lp, i ),
                                    glpk_strerror(lpx_get_row_type(lp,i))
                                    )
                             )
               )


  }
  
  
  names(r) <- c("name","status","prim","lb","ub","dual","strerr" )
  rownames( r ) <- 1:nrows

  r
}

  

#    No. Column name  St   Activity     Lower bound   Upper bound    Marginal
# ------ ------------ -- ------------- ------------- ------------- -------------
#      1 x(1)         B        109.153             0               
#      2 x(2)         NL             0             0                    -76.3559
#      3 x(3)         B            200             0               


get.col.report <- function( lp ) {

  ncols <- lpx_get_num_cols( lp )

  ## build the table for the columns
  c <- NULL
  for( i in 1:ncols ) {
    c <- rbind( c,
               as.data.frame( cbind( I(lpx_get_col_name( lp, i )),
                                    I(lpx_get_col_stat( lp, i ) ),
                                    I(lpx_get_col_prim( lp, i ) ),
                                    I(lpx_get_col_lb( lp, i ) ),
                                    I(lpx_get_col_ub( lp, i ) ),
                                    I(lpx_get_col_dual( lp, i ) ),
                                    I(lpx_get_col_b_ind( lp, i ) ),
                                    I(lpx_get_col_type( lp, i ) )
                                    )
                             )
               )
  }
  
  
  names(c) <- c("name","status","activity","lb","ub","dual","b_ind","type" )
  rownames( c ) <- 1:ncols

  c
}

