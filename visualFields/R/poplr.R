poplr <- function( vf, nperm = 5000, type = "slr", truncVal = 1,
                   sl_test = NULL, typecomb = "fisher",
                   details = FALSE ) {
##############
# input checks
##############
# check that all rows in vf belong to the same subject, the same test, the same perimetry
# testing and the same algorithm
  if( length( unique( vf$tperimetry ) ) > 1 |
      length( unique( vf$tpattern   ) ) > 1 |
      length( unique( vf$talgorithm ) ) > 1 |
      length( unique( vf$id ) ) > 1         |
      length( unique( vf$seye ) ) > 1 ) {
    stop( "all visual fields should belong to the same subject tested with the same perimeter and algorithm on the same locations" )
  }
  if( nperm < 5 ) stop( "number of permutations was lower than 5" )
  if( nperm > 1000000 ) stop( "please don't! Don't use more than a million permutations!" )
  if( ( type != "slr" & !is.null( sl_test ) ) ) stop( "tests about slopes being larger than a value are only valid for slr analysis" )
# truncation must be between zero and one
  if( truncVal <= 0 | truncVal > 1 ) stop("truncation must be between 0 and 1")

# permutation matrix

  porder <- make.permSpace( c( 1:nrow( vf ) ), nperm, return.permIDs = TRUE )$permID
  porder <- rbind( c( 1:nrow( vf ) ), porder )

  res <- NULL
# get last VF in res$vfdata
  res$vfdata      <- vf[nrow( vf ),]
# get and remove blind spot
  evaltxt <- paste("vfsettings$", vf$tpattern[1], "$bs", sep = "")
  bs <- eval(parse(text = evaltxt)) + visualFields::vfsettings$locini - 1
  res$vfdata <- res$vfdata[-bs]

  res$nvisits  <- nrow( vf )
  res$nperm    <- nperm
  res$type     <- type
  res$typecomb <- typecomb
# get the p-value statitics of the permuation analysis ...
  pstat    <- poplr_pstat( vf, porder = porder, type = "slr", sl_test = sl_test )
# ... and the actual analysis
  cstat    <- poplr_cstat( pstat$pval, typecomb = typecomb, truncVal = truncVal )
  if( type == "slr" ) {
    res$sl   <- pstat$sl[1,]
    res$int  <- pstat$int[1,]
    res$se   <- pstat$se[1,]
  } else if( type == "rank" ) {
    res$rho  <- pstat$rho[1,]
  }
  res$pval      <- pstat$pval[1,]
  res$scomb_obs <- cstat$scomb_obs
  res$pcomb_obs <- cstat$pcomb_obs
  res$pcomb     <- cstat$pcomb
  res$scomb     <- cstat$scomb
#  return detail or just final results?
  if( details ) {
    res$pstat <- pstat
    res$cstat <- cstat
  }
    
  return( res )
}
