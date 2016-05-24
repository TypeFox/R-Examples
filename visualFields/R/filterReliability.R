filterReliability <- function( vf, relCriteria = c( 0.2, 0.2, 0.2 ) ) {

# returns a report of visual fields that didn't pass reliability criteria or
# has NA
  idx_na <- which( is.na( vf$sfp ) |  is.na( vf$sfn ) | is.na( vf$sfl ) )
  idx_exc <- which( vf$sfp > relCriteria[1] | vf$sfn > relCriteria[2] | vf$sfl > relCriteria[3]  )
  
  return( union( idx_na, idx_exc ) )
}
