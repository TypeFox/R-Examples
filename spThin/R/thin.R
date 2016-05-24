#' @export thin
#' @title Spatially thin species occurence data
#' 
#' @description
#' \code{thin} returns spatially thinned species occurence data sets.
#' A randomizaiton algorithm (\code{\link{thin.algorithm}}) is used to create
#' data set in which all occurnece locations are at least \code{thin.par}
#' distance apart. Spatial thinning helps to reduce the effect of uneven,
#' or biased, species occurence collections on spatial model outcomes.
#' 
#' @param loc.data A data.frame of occurence locations. It can include several
#'   columnns, but must include at minimum a column of latitude and a
#'   column of longitude values
#' @param lat.col Name of column of latitude values. Caps sensitive.
#' @param long.col Name of column of longitude values. Caps sensitive.
#' @param spec.col Name of column of species name. Caps sensitive.
#' @param thin.par Thinning parameter - the distance (in kilometers) that you want
#'   records to be separated by.
#' @param reps The number of times to repete the thinning process. Given the random
#'   process of removing nearest-neighbors there should be 'rep' number of different
#'   sets of coordinates.
#' @param locs.thinned.list.return TRUE/FALSE - If true, the `list` of 
#'   the data.frame of thinned locs resulting from each replication 
#'   is returned (see Returns below).
#' @param write.files TRUE/FALSE - If true, new *.csv files will be 
#'   written with the thinned locs data
#' @param max.files The maximum number of *csv files to be written based on the
#'   thinned data
#' @param out.dir Directory to write new *csv files to
#' @param out.base A file basename to give to the thinned datasets created
#' @param write.log.file TRUE/FALSE create/append log file of thinning run
#' @param log.file Text log file 
#' @param verbose TRUE/FALSE - If true, running details of the function are print at the console.
#' 
#' @return locs.thinned.dfs A list of data.frames, each data.frame
#'   the spatially thinned locations of the algorithm for a 
#'   single replication. This list will have `reps` elements.
#'   
#' @seealso \code{\link{thin.algorithm}}
#'
thin <- function( loc.data, lat.col="LAT", long.col="LONG", spec.col="SPEC",
                  thin.par, reps,
                  locs.thinned.list.return = FALSE,
                  write.files = TRUE, 
                  max.files = 5, 
                  out.dir, 
                  out.base = "thinned_data",
                  write.log.file = TRUE,
                  log.file = 'spatial_thin_log.txt',
                  verbose = TRUE ){ 
  
  ## Begin writing to log file
  log.begin <- paste("**********************************************","\n",
                     "Beginning Spatial Thinning.\n",
                     "Script Started at:",
                     date(), sep=" ")
  ## Print information to the console
  if( verbose ){ cat( log.begin ) }
  ## Write information to the log.file
  if( write.log.file ){ write( log.begin, file=log.file, append = TRUE ) }
  
  ## Copy loc.data to new data.frame names locs.df
  locs.df <- loc.data
  
  ## Get the species name used in the `locs.df`
  species <- unique( locs.df[[ which( names(locs.df) == spec.col ) ]] )
  ## Send a warning message if there are more than one species in the df
  if( length( species ) > 1 ){
    log.spec.warn.1 <- "There appear to be more than one species name in this *.csv file."
    warning( log.spec.warn.1 )
    if( write.log.file ){ write( log.spec.warn.1, file=log.file, append=TRUE ) }
    species <- species[1]
    log.spec.warn.2 <- paste( "Only using species name:", species )
    warning( log.spec.warn.2 )
    if( write.log.file ){ write( log.spec.warn.2, file=log.file, append=TRUE ) }
  }
  
  ## Determine the columns associated with Lat and Long
  lat <- which( names(locs.df) == lat.col )
  long <- which( names(locs.df) == long.col )
  
  ## Make a data.frame that contains only the Long and 
  ## Lat values, in that order (ie df$Long, df$Lat)
  locs.long.lat <- as.data.frame( cbind( locs.df[[long]], locs.df[[lat]] ))

  ## Note in the log file what thinning parameter is being used
  log.thin.par <- paste("\nThinning Parameter Used (in km):", thin.par)
  if( write.log.file ){ write( log.thin.par, file=log.file, append = TRUE ) }
  log.num.reps <- paste("Number of replicates of thinning script:", reps )
  if( write.log.file ){ write( log.num.reps, file=log.file, append = TRUE ) }
  
  ## Execute spatial thinning function `thin.pres.data.R`. This
  ## function returns a `list` of spatially thinned data.frames
  
  # Keep track of how much time it takes to run this algorithm
  thin.time <- system.time( 
    locs.thinned <- thin.algorithm( rec.df.orig=locs.long.lat,
                                    thin.par=thin.par, reps=reps )
  )
  
  ## Record in log file elapsed system time for running the script
  if( write.log.file ){ write( "\nElapsed time for thinning completion", file=log.file, append = TRUE ) }
  if( write.log.file ){ write( thin.time, file=log.file, append = TRUE ) }
    
  ## Look at the number of locs kept in each thinned dataset
  ## by determining the number of rows in each returned data.frame
  lat.long.thin.count <- unlist(lapply(locs.thinned, nrow ))
  
  ## Use `table` to deterine number of dfs for each
  ## locs count
  locs.thinned.tbl <- table(lat.long.thin.count)
  ## Print `locs.thinned.tbl` to console
  if( verbose ){ cat("\n")
                 print(locs.thinned.tbl) }
  ## Print `locs.thinned.tbl` to log file
  if( write.log.file ){ write("\nNumber of data.frames per locations retained\nloc.cnt df.freq",
                              file=log.file, append=TRUE) }
  if( write.log.file ){ write(names(locs.thinned.tbl),file=log.file, append=TRUE, 
                              ncolumns=length(names(locs.thinned.tbl)),sep="\t") }
  if( write.log.file ){ write(locs.thinned.tbl, file=log.file, append=TRUE, 
                              ncolumns=length(locs.thinned.tbl),sep="\t") }
  
#   ## Plot a histogram of lat.long.thin.count
#   hist(lat.long.thin.count)
  
  ## Find max number of records
  max.thin.recs <- max( lat.long.thin.count)
  ## Save to log and Print this out for user to see
  log.max.rec <- paste( "Maximum number of records after thinning:",
                        max.thin.recs)
  if( verbose ){ print( log.max.rec ) }
  if( write.log.file ){ write( log.max.rec, file=log.file, append=TRUE) }
  
  ## Determine which data.frames
  ## have max.no. records
  max.dfs <- which( lat.long.thin.count == max.thin.recs)
  max.dfs.length <- length(max.dfs)
  log.max.df.cnt <- paste( "Number of data.frames with max records:", 
                           max.dfs.length)
  if( verbose ){ print( log.max.df.cnt ) }
  if( write.log.file ){ write(log.max.df.cnt, file=log.file, append=TRUE) }
  
  ## Write files if `write.files==TRUE`
  if( write.files ){
    if( verbose ){ print( "Writing new *.csv files" ) }
    if( write.log.file ){ write("\n**New *.csv file creation:**", file=log.file, append=TRUE) }
    
    # Determine number of files to write - should be the min
    # of the max number requested and the `max.dfs.lenght`
    n.csv <- min( c(max.files, max.dfs.length) )
    
    ## Write the first `n.csv` max data.frames
    # Check that `out.dir` exists. If not, create this directory.
    if ( !file.exists( out.dir ) ) {
      log.dir <- paste('Created new output directory: ', out.dir, sep="") 
      dir.create( out.dir, recursive=TRUE )
    } else {
      log.dir <- paste('Writing new *.csv files to output directory: ',
                       out.dir, sep="")
    }
    warning( log.dir )
    if( write.log.file ){ write( log.dir, file=log.file, append=TRUE ) }
    # Check that `out.dir` terminates in a '/'
    if( !grepl( '/$', out.dir ) ){
      out.dir <- paste( out.dir, '/', sep='' )
    }
    
    ## Make csv file names for thinned datasets
    csv.files <- paste( out.dir, out.base, "_thin", rep(1:n.csv), 
                        ".csv", sep="")
    
    for ( df in 1:n.csv ){
      # Get long and lat values for thinned locs
      df.temp <- locs.thinned[[ max.dfs[df] ]]
      # Add column of species name
      df.temp <- cbind( rep(as.character(species), max.thin.recs),
                        df.temp )
      # Give columns names
      colnames(df.temp) <- c(spec.col, long.col, lat.col)
      # Change file name in case of possible overwrite 
      # and send a warning to the user
      while( file.exists( csv.files[df] ) ){
        # Change file name
        csv.files[df] <- sub( '.csv$', '_new.csv', csv.files[df] )
        log.csv.overwrite <- paste(csv.files[df],
                                   "' already exists. Renaming file 
                                   to avoid overwriting.")
        warning( log.csv.overwrite )
        if( write.log.file ){ write(log.csv.overwrite, file=log.file, append=TRUE ) }
      }
      # Write new *.csv file with new name
      write.csv( df.temp, file=csv.files[df], quote=FALSE,
                 row.names=FALSE)
      log.write.file <- paste( "Writing file:", csv.files[df] )
      if( verbose ){ print( log.write.file ) }
      if( write.log.file ){ write( log.write.file, file=log.file, append=TRUE ) }
    }
    
  } else {
    log.write.file <- "No files written for this run."
    if( verbose ){ print( log.write.file ) }
    if( write.log.file ){ write( log.write.file, file=log.file, append=TRUE) }
  }
  
  ## Return `locs.thinned.list` if that setting is TRUE
  if ( locs.thinned.list.return ){
    return( locs.thinned )
  }
}
