#' Method initialize for class ViSibook object.
#' @title Method \code{initialize-ViSibook}
#' @rdname initialize-ViSibook-methods
#' @aliases initialize,ViSibook-method
#' @exportMethod initialize
#' @docType methods
#' @param .Object a ViSibook object.
#' @param vars a vector storing names of actions.
#' @param label a vector. storing brief description of actions.
#' @param typeA Vector storing type of actions, "l" for long actions, "p" for punctuals.
#' @param showorder vector storing order in which actions will be plotted, is an actions is not to be plot its showorder is "NA".
#' @param deb Vector storing, for long actions, the punctual action names that corresponds to its start.
#' @param fin Vector storing, for long actions, the punctual action that corresponds to its end.
#' @param GZDeb Vector storing punctuals actions green zone starting time.
#' @param GZFin Vector storing punctual action green zone ending time.
#' @param Repetition Vector storing if the green zones should be repeated the time interval of repetition.
#' @param BZBeforeDeb Vector storing punctual black zone 1 starting time.
#' @param BZBeforeFin Vector storing punctual black zone 1 ending time.
#' @param BZAfterDeb Vector storing punctual black zone 2 starting time. 
#' @param BZAfterFin Vector storing punctual black zone 2 ending time. 
#' @param BZLong Vector storing the long action black zone time. 
#' @param BZLtype Vector storing the type of the black zone,
#' "time" if the action should be finish at a time, "span" if the action should be finish in a time.
#' @param NAMES Vector storing names of slots that are to be considered for \code{\link{plot-ViSigrid-method}}.
#' @return a ViSibook object
#' @seealso See \code{\link{plot-ViSigrid-method}} for examples.
setMethod("initialize", signature("ViSibook"),  
          function(.Object,vars,label,typeA,showorder,deb,
                   fin,GZDeb,GZFin,Repetition,BZBeforeDeb,
                   BZBeforeFin,BZAfterDeb,BZAfterFin,BZLong,BZLtype,NAMES){
  ## Minimum structure 
  #### Checking for null or  that will disturbe other script (build.ViSigrid and plot.ViSigrid)
  #lengths
  if (any( rep( length( vars ) , 5 ) != c( length( label ) , length( typeA ), length( showorder ), length( deb ) , length( fin ) ) ) ) {
    stop( " initialize ( ViSibook ) : length of vars,label,showorder,deb,fin are not equals \n" )
  }
  ### vars	
  if (any(c( is.na( vars ),is.null( vars ) ) ) ) { stop( " initialize ( ViSibook ) : vars can not be NA or NULL \n" ) }else{ 
  }
  ### label
  if (any(c( is.na( label ) , is.null( label ) ) ) ) { 
    label[ which( is.na(label) ) ] <- vars[ which( is.na( label ) ) ] 
    warning( " initialize ( ViSibook ) : NA or NULL in label remplaced by vars values \n" ) 
  }
  ### typeA
  if (c( any( is.na( typeA ) ) ) ) { stop( " initialize ( ViSibook ) : typeA can not be NA or NULL \n " ) }else{ 
    ####Checking typeA value "p or l"
    if (any(((typeA == "p" | typeA == "l" ) | is.na( typeA ) ) == FALSE ) ) {
      stop( " initialize ( ViSibook ) : Error typeA should be \"p\" or \"l\" \n ") 
    }
  }	
  #### Showorder 
  if (is.integer( showorder ) == FALSE ) {
    showorder <- as.integer(showorder)
  }
  if (length( unique( showorder[ which( is.na( showorder ) == FALSE ) ] ) ) != length( showorder[ which( is.na( showorder ) == FALSE ) ] ) ) {
  stop( " initialize ( ViSibook ) : Error in showorder has one or more duplicates  \n " ) 
  }    
  ### deb
  ####Checking deb def for each long action 
  if (any( is.na( deb[ which( typeA == "l" ) ] ) ) ) { stop( " initialize ( ViSibook ) : Not all deb are defined for type action long \n " ) }else{deb <- deb }
  #### Names in deb and fin Matching with vars
  temp <- vars[ which( typeA == "l" ) ][ which( unlist( lapply( deb[ which( typeA == "l" ) ] , function(x )(is.na( match( x, vars[ which( typeA == "p" ) ] ) ) ) ) ) == TRUE ) ]
  if (length( temp ) > 0 ) {stop( paste( " initialize ( ViSibook ) : Error ", " in deb type long action(s) ", temp ," do not match with any punctual name action in vars" ) ) }
  ### fin 
  ####Checking deb def for each long action 
  if (any( is.na( fin[ which( typeA == "l" ) ] ) ) ) { stop( " initialize ( ViSibook ) : Not all fin are defined for type action long \n " ) }else{fin <- fin }
  # Names in fin Matching with vars
  temp <- vars[ which( typeA == "l" ) ][ which( unlist( lapply( fin[ which( typeA == "l" ) ] , function(x )(is.na( match( x, vars[ which( typeA == "p" ) ] ) ) ) ) ) == TRUE ) ]
  if (length( temp ) > 0 ) {stop( paste( " initialize ( ViSibook ) : Error ", " in fin type long action(s) ", temp ," do not match with any punctual name action in vars" ) ) }
  methods::slot( .Object , "vars" ) <- vars
  methods::slot( .Object , "label" ) <- label
  methods::slot( .Object , "typeA" ) <- typeA
  methods::slot( .Object , "showorder" ) <- as.numeric( showorder )
  methods::slot( .Object , "deb" ) <- deb
  methods::slot( .Object , "fin" ) <- fin
  methods::slot( .Object , "NAMES" ) <- c("vars","label","typeA","showorder","deb","fin")
  #Options#
  # Green Zone...........................................................................................................
  if (any( is.na( c( GZDeb , GZFin ) ) == FALSE ) ) {
    if (any( rep( length( vars ) , 2 ) != c( length( GZDeb ) , length( GZFin ) ) ) ) {
      Repetition <- rep( NA , length( vars ) )
      warning( " initialize ( ViSibook ) : Length of GZDeb and/or GZFin are not equals to the length of vars \n " )
      warning(" initialize ( ViSibook ) : No green zone defined for punctuals actions \n ")
    }else{
      #### Checking fin' zones are > to def' zones
      if (any( (as.numeric(GZDeb) >= as.numeric(GZFin) ) , na.rm = TRUE ) ) {
        GZDeb[ which( (as.numeric(GZDeb) < as.numeric(GZFin) ) ) ] <- rep( NA , sum( (as.numeric(GZDeb) < as.numeric(GZFin) ), na.rm = TRUE ) )
        GZFin[ which( (as.numeric(GZDeb) < as.numeric(GZFin) ) ) ] <- rep( NA , sum( (as.numeric(GZDeb) < as.numeric(GZFin) ), na.rm = TRUE ) )
        warning( "  initialize ( ViSibook ) : when as.numeric(GZDeb) >= as.numeric(GZFin) values replaced by NA \n " )
      }
      #### Checking  "GZDeb" and "GZFin" either both empty either both defined
      if (any( is.na( GZDeb ) != is.na( GZFin ) ) ) { 
        GZDeb[ which( is.na( GZDeb ) != is.na( GZFin )  ) ] <- rep( NA,sum( is.na( GZDeb ) != is.na( GZFin )  )) 
        GZDeb[ which( is.na( GZDeb ) != is.na( GZFin )  ) ] <- rep( NA,sum( is.na( GZDeb ) != is.na( GZFin )  )) 
        warning( "initialize ( ViSibook ) : For action(s) ", vars[which( is.na( GZDeb ) != is.na( GZFin )  ) ] ," when only one of GZDeb and GZFin defined value(s) remplaced by NA \n " )
      }
      if (any( is.na( c( GZDeb , GZFin ) ) == FALSE ) ) {
        methods::slot( .Object , "GZDeb" ) <-  GZDeb 
        methods::slot( .Object , "GZFin" ) <-  GZFin 
        methods::slot( .Object , "NAMES" ) <- c( methods::slot( .Object , "NAMES" ) , "GZDeb" )
        methods::slot( .Object , "NAMES" ) <- c( methods::slot( .Object , "NAMES" ) , "GZFin" ) 
      }
      #### Checking If repetition defined and correctly
      if (any( is.na( Repetition ) == FALSE )  ) {
                 slot( .Object , "Repetition" ) <-  Repetition 
                 methods::slot( .Object , "NAMES" ) <- c( methods::slot( .Object , "NAMES" ) , "Repetition" )
        }
    }
}else{
    Repetition <- rep( NA , length( vars ) )
    warning("No green zone defined for punctuals actions \n ")
  }
  # Green zone fin..................................................................................................................
  #.................................................................................................................................
  # Black zone 1 (Before)...........................................................................................................
  if (any( is.na( c( BZBeforeDeb , BZBeforeFin ) ) == FALSE ) ) {
    ### Checking lengths
    if (any( rep( length( vars ) , 2 ) != c( length( BZBeforeDeb ) , length( BZBeforeFin ) ) ) ) {
      BZBeforeDeb <- rep( NA , length( vars ) )
      BZBeforeFin <- rep( NA , length( vars ) )
      warning( "initialize ( ViSibook ) Length of BZBeforeDeb and BZBeforeFin are not equal to the length of vars \n " )
    }else{
      ### Checking No repetition and Black zone define for the same Action
      if (any( is.na( Repetition[ which( is.na( BZBeforeDeb ) == FALSE ) ] ) == FALSE ) ) {
        BZBeforeDeb[ which( is.na( Repetition ) == FALSE ) ] <- rep( NA , sum( is.na( Repetition ) == FALSE )  )
      }
      if (any( is.na( Repetition[ which( is.na( BZBeforeFin) == FALSE ) ] ) == FALSE ) ) {
        BZBeforeFin[ which( is.na( Repetition ) == FALSE ) ] <- rep( NA , sum( is.na( Repetition ) == FALSE )  )
      }
      #### Checking  "BZBeforeDeb" and "BZBeforeFin" either both empty either both defined
      if (any( is.na( BZBeforeDeb ) != is.na( BZBeforeFin )  ) ) { 
        BZBeforeDeb[which( is.na( BZBeforeDeb ) != is.na( BZBeforeFin )  )] <- rep( NA,sum( is.na( BZBeforeDeb ) != is.na( BZBeforeFin )  )) 
        BZBeforeFin[which( is.na( BZBeforeDeb ) != is.na( BZBeforeFin )  )] <- rep( NA,sum( is.na( BZBeforeDeb ) != is.na( BZBeforeFin )  ))
        warning( " initialize ( ViSibook ) : For action, ", vars[which( is.na( BZBeforeDeb ) != is.na( BZBeforeFin )  )] ," only BZBeforeDeb or BZBeforeFin, value is remplaced by NA \n " )
      }
      ####	Checking Deb < Fin
      if (any( (as.numeric(BZBeforeDeb) >= as.numeric(BZBeforeFin) ) , na.rm = TRUE ) ) {
        BZBeforeDeb[ which( (as.numeric(BZBeforeDeb) < as.numeric(BZBeforeFin) ) ) ] <- rep( NA , sum( (as.numeric(BZBeforeDeb) < as.numeric(BZBeforeFin) ), na.rm = TRUE ) )
        BZBeforeFin[ which( (as.numeric(BZBeforeDeb) < as.numeric(BZBeforeFin) ) ) ] <- rep( NA , sum( (as.numeric(BZBeforeDeb) < as.numeric(BZBeforeFin) ), na.rm = TRUE ) )
        warning( "  initialize ( ViSibook ) : when as.numeric(BZBeforeDeb) >= as.numeric(BZBeforeFin) replaced by NA \n " )
      }
      if (any( is.na( c( BZBeforeDeb , BZBeforeFin ) ) == FALSE ) ) {
        methods::slot( .Object , "BZBeforeDeb" ) <-  BZBeforeDeb 
        methods::slot( .Object , "BZBeforeFin" ) <-  BZBeforeFin 
        methods::slot( .Object , "NAMES" ) <- c( methods::slot( .Object , "NAMES" ) , "BZBeforeDeb" )
        methods::slot( .Object , "NAMES" ) <- c( methods::slot( .Object , "NAMES" ) , "BZBeforeFin" ) 
      }
    }
  }
  # END Black zone 1 (Before)...........................................................................................................
  #.................................................................................................................................
  # Black zone 2 (After)...........................................................................................................
  if (any( is.na( c( BZAfterDeb , BZAfterFin ) ) == FALSE ) ) {
    ### Checking lengths
    if (any( rep( length( vars ) , 2 ) != c( length( BZAfterDeb ) , length( BZAfterFin ) ) ) ) {
      BZAfterDeb <- rep( NA , length( vars ) )
      BZAfterFin <- rep( NA , length( vars ) )
      warning( "initialize ( ViSibook ) Length of BZAfterDeb and BZAfterFin are not equal to the length of vars \n " )
    }else{
      ### Checking No repetition and Black zone define for the same Action
      if (any( is.na( Repetition[ which( is.na( BZAfterDeb ) == FALSE ) ] ) == FALSE ) ) {
        BZAfterDeb[ which( is.na( Repetition ) == FALSE ) ] <- rep( NA , sum( is.na( Repetition ) == FALSE )  )
      }
      if (any( is.na( Repetition[ which( is.na( BZAfterFin) == FALSE ) ] ) == FALSE ) ) {
        BZAfterFin[ which( is.na( Repetition ) == FALSE ) ] <- rep( NA , sum( is.na( Repetition ) == FALSE )  )
      }
      #### Checking  "BZAfterDeb" and "BZAfterFin" either both empty either both defined
      if (any( is.na( BZAfterDeb ) != is.na( BZAfterFin )  ) ) { 
        BZAfterDeb[which( is.na( BZAfterDeb ) != is.na( BZAfterFin )  )] <- rep( NA,sum( is.na( BZAfterDeb ) != is.na( BZAfterFin )  )) 
        BZAfterFin[which( is.na( BZAfterDeb ) != is.na( BZAfterFin )  )] <- rep( NA,sum( is.na( BZAfterDeb ) != is.na( BZAfterFin )  ))
        warning( " initialize ( ViSibook ) : For action, ", vars[which( is.na( BZAfterDeb ) != is.na( BZAfterFin )  )] ," only BZAfterDeb or BZAfterFin, value is remplaced by NA \n " )
      }
      ####	Checking Deb < Fin
      if (any( (as.numeric(BZAfterDeb) >= as.numeric(BZAfterFin) ) , na.rm = TRUE ) ) {
        BZAfterDeb[ which( (as.numeric(BZAfterDeb) < as.numeric(BZAfterFin) ) ) ] <- rep( NA , sum( (as.numeric(BZAfterDeb) < as.numeric(BZAfterFin) ), na.rm = TRUE ) )
        BZAfterFin[ which( (as.numeric(BZAfterDeb) < as.numeric(BZAfterFin) ) ) ] <- rep( NA , sum( (as.numeric(BZAfterDeb) < as.numeric(BZAfterFin) ), na.rm = TRUE ) )
        warning( "  initialize ( ViSibook ) : when as.numeric(BZAfterDeb) >= as.numeric(BZAfterFin) replaced by NA \n " )
      }
      
      if (any( is.na( c( BZAfterDeb , BZAfterFin ) ) == FALSE ) ) {
        methods::slot( .Object , "BZAfterDeb" ) <-  BZAfterDeb
        methods::slot( .Object , "BZAfterFin" ) <-  BZAfterFin
        methods::slot( .Object , "NAMES" ) <- c( methods::slot( .Object , "NAMES" ) , "BZAfterDeb" )
        methods::slot( .Object , "NAMES" ) <- c( methods::slot( .Object , "NAMES" ) , "BZAfterFin" ) 
      }
    }
  }
  # END Black zone 2 (After)...........................................................................................................
  #.................................................................................................................................
  # Black zone long ...........................................................................................................
  if (any( is.na( BZLong ) == FALSE ) ) {
    ### Checking lengths
    if (length( vars ) != length( BZLong )  ) {
      BZLong <- rep( NA , length( vars ) )
      warning( "initialize ( ViSibook ) Length of BZLong is not equal to the length of vars \n " )
    }else{
      ####Checking BZLtype value "span or time"
      if (any( (BZLtype == "span" | BZLtype == "time" | is.na( BZLtype ) ) == FALSE ) ) {
        BZLtype[ which( (BZLtype == "span" | BZLtype == "time" | is.na( BZLtype ) ) == FALSE ) ] <- rep( NA , sum( (BZLtype == "span" | BZLtype == "time" | is.na( BZLtype ) ) == FALSE , na.rm = TRUE ) )
        warning( "  initialize ( ViSibook ) : BZLtype should be \"span\" or \"time\", unrecognized values replaces by NA \n ") 
      }	
      if (any( is.na( BZLtype[ which( is.na( BZLong ) == FALSE ) ] ) ) ) {		
        BZLtype[ which( is.na( BZLong) == FALSE ) ][ which( is.na( BZLtype[ which( is.na( BZLong ) == FALSE ) ] ) )] <- rep( "time" , sum( is.na( BZLtype[ which( is.na( BZLong ) == FALSE ) ] ) ))
        warning( "  initialize ( ViSibook ) : when BZLong not NA and BZLtype NA, BZLtype is set to \"time\" \n " )
      }
      if (any( is.na( BZLong ) == FALSE ) ) {
        methods::slot( .Object , "BZLong" ) <-  BZLong 
        methods::slot( .Object , "NAMES" ) <- c( slot( .Object , "NAMES" ) , "BZLong" )
        methods::slot( .Object , "BZLtype" ) <- BZLtype
        methods::slot( .Object , "NAMES" ) <- c( slot( .Object , "NAMES" ) , "BZLtype" )
      }
    }
  }
  return( .Object ) 
}
)
