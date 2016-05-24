# -----------------------------------------------------------------------------
# FUNCTION: dummy.data.frame:
#
#   Given a _name_ of a column and a data frame, return the data.frame for that
#
#   x          : a data.frame, matrix or single variable or variable name
#
#   data       : An object such as a matrix or data.frame with colnames.
#                If provided, x is take as the the name of a column on the data
#
#   sep        : the seperator to be used between the variable name and the 
#                value.  Used in new variable construction.  Also, set as an
#                attribute
#
#   drop       : Drop unused levels?.  When x is a factor, whether to produce 
#                dummy variable for only the used levels. If x has unused 
#                levels and drop=T ( the default ), dummy variables will not be
#                created for the values of y and not the levels.
# 
#   constant   : Whether to return an identity vectors for variables that 
#                assume one value.
#                
#   fun        : Function to coerce the value in the final matrix.  
#                Default: 'as,integer'
#
#   verbose    : logical.  Whether to print(cat) the number of variables 
#                Default: FALSE
#   
#   NA         : Options: 
#                 encode as a seperate dummy variable default.
#                 create a row of NA for that observation --> model.frame(
#                   na.action=na.pass )
#                 omit
#
#  TODO:
#   - If continuous variables, allow for quantile encoding.  
#   - 
# -----------------------------------------------------------------------------

dummy <- function( x, data=NULL, sep="", drop=TRUE, fun=as.integer, verbose = FALSE ) { 


  # HANDLE IF DATA IS MISSING.  
    if( is.null(data) ) {
      name <- as.character( sys.call(1) )[2]   
      name <- sub( "^(.*\\$)", "", name )    # REMOVE prefix e.f
      name <- sub( "\\[.*\\]$", "", name )   # REMOVE suffix   
    } else {
      if( length(x) > 1 ) stop( "More than one variable provided to produce dummy variable." )  
      name <- x
      x    <- data[ , name]
    }


  # CHANGE TO FACTOR: KEEP LEVELS?
    if( drop == FALSE && class(x) == "factor" ) {
      x <- factor( x, levels=levels(x), exclude=NULL ) 
    } else {
      x<-factor( x, exclude=NULL )
    }
   

  # TRAP FOR ONE LEVEL :  
  #   model.matrix does not work on factor w/ one level.  Here we trap for the spacial case.
    if( length(levels(x))<2 ) {
      
      if( verbose ) warning( name, " has only 1 level. Producing dummy variable anyway." )

      return(          
        matrix( 
          rep(1,length(x)), 
          ncol=1, 
          dimnames=list( rownames(x), c( paste( name, sep, x[[1]], sep="" ) ) ) 
        )
      )

    }


  # GET THE MODEL MATRIX   
    mm <- model.matrix( ~ x - 1, model.frame( ~ x - 1 ),  contrasts=FALSE )  # vec
    colnames.mm <- colnames(mm) 

    if( verbose ) cat( " ", name, ":", ncol(mm), "dummy varibles created\n" ) 

    mm <- matrix( fun(mm), nrow=nrow(mm), ncol=ncol(mm), dimnames=list(NULL, colnames.mm) ) 


  # Replace the column names 'x'... with the true variable name and a seperator
    colnames(mm) <- sub( "^x", paste( name, sep, sep="" ), colnames(mm) )
    if(! is.null(row.names(data)) ) rownames(mm) <- rownames(data)

    return(mm)   

}

                                                           

# TESTING:  See Rd Files

