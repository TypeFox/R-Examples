# -----------------------------------------------------------------------------
# FUNCTION: dummy.data.frame
#  Produce a dummy data.frame, i.e. where all categorical ( non-continous ) 
#  variables are expanded. 
#  
#  all  : return all columns or only the categorical variables.
#  names: names of cols to expand as dummy variables.
#
#  TODO:
#   - matrix?
#   - na.action
#
# -----------------------------------------------------------------------------

dummy.data.frame <- function( data, names=NULL, omit.constants = TRUE, dummy.classes=getOption("dummy.classes"), all=TRUE, ... ) {

  # Initialize the data.frame
    df<-data.frame( row.names=row.names(data) )     
    new.attr <- list()  # Track location of dummy variables

    for( nm in names(data) ) {
      
# cat( nm )
      old.attr <- attr(df,'dummies')
      
      if(
        nm %in% names || 
        ( is.null(names) && ( dummy.classes == "ALL" || class(data[,nm]) %in% dummy.classes ) )
      ) {

        dummies <- dummy( nm, data, ... )

        # OMIT CONSTANT COLUMNS:
        #  Variables that are constant will return a matrix with one column
        if( ncol(dummies) == 1  & omit.constants ) {
          dummies <- matrix( nrow=nrow(data), ncol=0 ) 
        }
            
        if( ncol(dummies)>0 ) new.attr[[nm]] <- (ncol(df)+1):( ncol(df)+ncol(dummies) ) 

      } else {
        if( ! all ) next()
        dummies <- data[,nm, drop=FALSE ]
      }

      df <- cbind(df, dummies)

    }

    attr( df, 'dummies' ) <- new.attr
    return(df)

}
      
    
# TESTING:
# dummy.data.frame(iris)
# dummy.data.frame(iris all=FALSE)



