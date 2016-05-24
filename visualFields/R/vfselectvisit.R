vfselectvisit <- function( vf, sel = "last", numTests = 1,
                           beginDate = NA, endDate = NA ) {

# check input
  if( is.numeric( sel ) & !is.vector( sel ) ) stop( "wrong type of selection" )
  if( !is.numeric( sel ) && ( sel != "last" & sel != "first" & sel != "date" ) ) stop( "wrong type of selection" )

# first all sort data
  if( !is.numeric( sel ) && sel == "last") {
    vf <- vfsort( vf, decreasing = TRUE )
    sel <- "first"
  } else {
    vf <- vfsort( vf )
  }

# set dates if they are NA and sel equals "date"
  if( !is.numeric( sel ) && ( sel == "date" & is.na( beginDate ) ) ) beginDate <- "1900-01-01"
  if( !is.numeric( sel ) && ( sel == "date" & is.na( endDate ) ) )   endDate   <- format(Sys.time(), "%Y-%m-%d")
  beginDate <- as.Date( beginDate )
  endDate   <- as.Date( endDate )
# get unique records
  uid      <- NULL
  uid$id   <- vf$id
  uid$seye <- vf$seye
  uid      <- as.data.frame( uid )
  uid      <- unique( uid )

  vf_select <- NULL
  if( is.numeric( sel ) ) idx <- sel
# perform subselection
  for( i in 1:nrow( uid ) ) {
    vfidx <- which( vf$id == uid$id[i] & vf$seye == uid$seye[i] )
# too much to select
    if( !is.numeric( sel ) && sel == "first" ) {
      if( length( vfidx ) < numTests ) next
      idx <- c( 1:numTests )
    } else if( !is.numeric( sel ) && sel == "date" ) {
      idx <- which( vf$tdate[vfidx] >= beginDate & vf$tdate[vfidx] <= endDate )
    }
# perform actual selection
    vfidx <- vfidx[idx]

    if( length( vfidx ) > 0 ) vf_select <- rbind( vf_select, vf[vfidx,] )
  }

  if( is.null( vf_select ) ) return( NA )
  return( vfsort( vf_select ) )

}