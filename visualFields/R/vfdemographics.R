vfdemographics <- function( vf ) {

  demog <- NULL
# from a set of visual fields it fits lines to characterize age effect on the visual-field sensitivities
# For this function all visual fields should correspond to the same device, algorithm, and pattern of
# test locations. If not, stop!!!
  if( length( unique( vf$tperimetry ) ) > 1 ) {
    stop("mixing different perimeters data")
  }
  if( length( unique( vf$talgorithm ) ) > 1 ) {
    stop("mixing different algorithm data")
  }
  if( length( unique( vf$tpattern ) ) > 1 ) {
    stop("mixing different patterns of locations")
  }

# get weights based on number of visits per subject
  idu <- NULL
  idu$id <- unique( vf$id )
  for( i in 1:length( idu$id ) ) {
    idu$weight[i] <- 1 / length( which( vf$id == idu$id[i] ) )
  }
  idweight <- NULL
  for( i in 1:length( vf$id ) ) {
    idweight[i] <- idu$weight[which( idu$id == vf$id[i] )]
  }

# get settings for the pattern of test locations
  locini   <- visualFields::vfsettings$locini
  texteval <- paste( "vfsettings$", vf$tpattern[1], sep = "" )
  settings <- eval( parse( text = texteval ) )

# calculates demographic statistics of subjects and visits of data frame
  nvisit <- NULL
  futime <- NULL
  for( i in 1:length( idu$id ) ) {
    idx <- which( vf$id == idu$id[i] )
    nvisit[i] <- length( idx )
    futime[i] <- difftime( vf$tdate[idx[length( idx )]], vf$tdate[idx[length( 1 )]], units = "days" )
  }
  demog$nsubjects    <- length( unique( idu$id) )
  demog$ntotalvisits <- nrow( vf )
# number of visits per subject and average follow up time
  demog$stats$nvisitspersubject <- as.numeric( summary( nvisit ) )
  demog$stats$futime <- as.numeric( summary( futime ) )
# age of subjects: weighted average
  demog$stats$age[1] <- min( vf$sage )
  demog$stats$age[2] <- wtd.quantile( vf$sage, probs = 0.25, weights = idweight, normwt = TRUE )
  demog$stats$age[3] <- wtd.quantile( vf$sage, probs = 0.5, weights = idweight, normwt = TRUE )
  demog$stats$age[4] <- weighted.mean( vf$sage, w = idweight )
  demog$stats$age[5] <- wtd.quantile( vf$sage, probs = 0.75, weights = idweight, normwt = TRUE )
  demog$stats$age[6] <- max( vf$sage )
# false positives
  demog$stats$sfp[1] <- min( vf$sfp )
  demog$stats$sfp[2] <- wtd.quantile( vf$sfp, probs = 0.25, weights = idweight, normwt = TRUE )
  demog$stats$sfp[3] <- wtd.quantile( vf$sfp, probs = 0.5, weights = idweight, normwt = TRUE )
  demog$stats$sfp[4] <- weighted.mean( vf$sfp, w = idweight )
  demog$stats$sfp[5] <- wtd.quantile( vf$sfp, probs = 0.75, weights = idweight, normwt = TRUE )
  demog$stats$sfp[6] <- max( vf$sfp )
# false negatives, fixation losses
  demog$stats$sfn[1] <- min( vf$sfn )
  demog$stats$sfn[2] <- wtd.quantile( vf$sfn, probs = 0.25, weights = idweight, normwt = TRUE )
  demog$stats$sfn[3] <- wtd.quantile( vf$sfn, probs = 0.5, weights = idweight, normwt = TRUE )
  demog$stats$sfn[4] <- weighted.mean( vf$sfn, w = idweight )
  demog$stats$sfn[5] <- wtd.quantile( vf$sfn, probs = 0.75, weights = idweight, normwt = TRUE )
  demog$stats$sfn[6] <- max( vf$sfn )
# fixation losses
  demog$stats$sfl[1] <- min( vf$sfl )
  demog$stats$sfl[2] <- wtd.quantile( vf$sfl, probs = 0.25, weights = idweight, normwt = TRUE )
  demog$stats$sfl[3] <- wtd.quantile( vf$sfl, probs = 0.5, weights = idweight, normwt = TRUE )
  demog$stats$sfl[4] <- weighted.mean( vf$sfl, w = idweight )
  demog$stats$sfl[5] <- wtd.quantile( vf$sfl, probs = 0.75, weights = idweight, normwt = TRUE )
  demog$stats$sfl[6] <- max( vf$sfl )
# as data frame
  demog$stats <- as.data.frame( demog$stats )
  nc <- ncol( demog$stats )
# stats of locations
  for( i in locini:( locini - 1 + settings$locnum ) ) {
    nc <- nc + 1
    demog$stats[1,nc] <- min( vf[,i] )
    demog$stats[2,nc] <- wtd.quantile( vf[,i], probs = 0.25, weights = idweight, normwt = TRUE )
    demog$stats[3,nc] <- wtd.quantile( vf[,i], probs = 0.5, weights = idweight, normwt = TRUE )
    demog$stats[4,nc] <- weighted.mean( vf[,i], w = idweight )
    demog$stats[5,nc] <- wtd.quantile( vf[,i], probs = 0.75, weights = idweight, normwt = TRUE )
    demog$stats[6,nc] <- max( vf[,i] )
    names( demog$stats )[nc] <- names( vf )[i]
  }
  row.names( demog$stats ) <- c( "Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max" )

  return( demog )

}
