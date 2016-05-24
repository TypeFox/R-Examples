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
#  nrows <- Get_num_rows( lp )
  nrows <-  getNumRowsGLPK(lp)  
  ## build the table for the rows
  r <- NULL
  i <- NULL
  for( i in 1:nrows ) {
    r <- rbind( r,
               as.data.frame( cbind( getRowNameGLPK(lp, i),
                                     getRowStatGLPK( lp, i ),
                                     getRowPrimGLPK( lp, i ),
                                     getRowLowBndGLPK( lp, i ),
                                     getRowUppBndGLPK( lp, i ),
                                     getRowDualGLPK( lp, i )
                                    )
                             )
               )


  }
  
  
#  names(r) <- c("name","status","prim","lb","ub","dual","strerr" )
  names(r) <- c("name","status","prim","lb","ub","dual")  
  rownames( r ) <- 1:nrows

  r
}

  

#    No. Column name  St   Activity     Lower bound   Upper bound    Marginal
# ------ ------------ -- ------------- ------------- ------------- -------------
#      1 x(1)         B        109.153             0               
#      2 x(2)         NL             0             0                    -76.3559
#      3 x(3)         B            200             0               


get.col.report <- function( lp ) {

  ncols <- getNumColsGLPK( lp )

  ## build the table for the columns
  c <- NULL
  i <- NULL
  for( i in 1:ncols ) {
    c <- rbind( c,
               as.data.frame( cbind( I(getColNameGLPK( lp, i ) ),
                                     I(getColStatGLPK( lp, i ) ),
                                     I(getColPrimGLPK( lp, i ) ),
                                     I(getObjCoefGLPK( lp, i ) ),
                                     I(getColLowBndGLPK( lp, i ) ),
                                     I(getColUppBndGLPK( lp, i ) ),
                                     I(getColDualGLPK( lp, i ) ),
#                                     I(getCbindGLPK( lp, i ) ),
# Error in routine getCbindGLPK reported to Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de> March 24th                                      
                                     I(getColTypeGLPK( lp, i ) )
                                    )
                             )
               )
print(i)
print(c)
    }
  
  
#  names(c) <- c("name","status","activity","coef", "lb","ub","dual","b_ind","type" )
  names(c) <- c("name","status","activity","coef","lb","ub","dual","type" )
  
#  colnames( c ) <- 1:ncols
#  Error in `colnames<-`(`*tmp*`, value = 1:48) : 
#    'names' attribute [48] must be the same length as the vector [7]

  c
}

