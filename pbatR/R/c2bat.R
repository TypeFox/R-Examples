## alright, seriously, why would you ever force a user to make a snp a number?
c2bat <- function( snps,
                   case.aa, case.Aa, case.AA,
                   control.aa, control.Aa, control.AA,
                   replicates=0,
                   statistic="armitage",
                   method="exact" ) {
  ## fix up the parameters and check the users input
  if( statistic=="armitage" ){
    statistic <- 1;
  }else if( statistic=="logrank" ){
    statistic <- 2;
  }else{
    stop( paste( "'statistic' must be either 'armitage' or 'exact'; you specified '", statistic, "', which is not supported.", sep="" ) );
  }

  if( method=="exact" ){
    method <- 1;
  }else if( method=="approximate" ){
    method <- 0;
  }

  if( replicates<0 )
    stop( "Replicates should be 0 for the first analysis, and then > 1000 for the selected SNPs afterward." );

  if( length(case.aa) != length(case.Aa) ||
      length(case.aa) != length(case.AA) ||
      length(case.aa) != length(control.aa) ||
      length(case.aa) != length(control.Aa) ||
      length(case.aa) != length(control.AA) ) {
    stop( "length of the vectors of genotypes for cases and controls must all be the same." );
  }

  ## now lets get it done.

  ## kill the cc.txt file
  ##file.remove( "cc.txt" ); ## leftovers...
  ## be civil, maybe back it up?
  if( file.exists( "cc.txt" ) ) {
    file.copy( from="cc.txt", to="cc.txt.bak", overwrite=TRUE );
    file.remove( "cc.txt" );
  }
  
  ## now loop through each of them
  ## -- maybe we should launch this from c++ in the future?
  for( i in 1:length(case.aa) )
    system( paste( pbat.get(), i, case.aa[i], case.Aa[i], case.AA[i], control.aa[i], control.Aa[i], control.AA[i], replicates, statistic, method ) );

  ## now read in the output
  res <- read.table( "cc.txt" );
  res[,1] <- snps;  ## allow the user to use names for snps if they like
  names(res) <- c( "snp", "case.aa", "case.Aa", "case.AA", "control.aa", "control.Aa", "control.AA", "pvalue.Monte", "pvalue", "noncentrality", "modelc2.or", "allelic.or");

  ## and return
  return( res[,1:12] );
}
