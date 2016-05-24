#########################################################
# CGB, 20100707, created
#########################################################

jython.call <- function( rJython, py.foo, ... ){

    foo.args <- list( ... )

    if( is.null( names( foo.args ) ) )
        which.dict <- rep( FALSE, length( foo.args ) )
    else
        which.dict <- names( foo.args ) != ""

    n.args.vect <- sum( !which.dict )
    n.args.dict <- sum(  which.dict )

    foo.args.dict <- toJSON( foo.args[  which.dict ] )
    foo.args.vect <- toJSON( foo.args[ !which.dict ] )

    # Passing data to jython

    #rJython$exec( paste( "_r_args_dict ='", foo.args.dict, "'", sep = "" ) )
    #rJython$exec( paste( "_r_args_vect ='", foo.args.vect, "'", sep = "" ) )
    #rJython$exec( "_r_args_dict = json.loads( _r_args_dict )" )
    #rJython$exec( "_r_args_vect = json.loads( _r_args_vect )" )

    # Creating the call

    jython.command <- c( 
        paste( "_r_args_dict ='", foo.args.dict, "'", sep = "" ),
        paste( "_r_args_vect ='", foo.args.vect, "'", sep = "" ),
        "_r_args_dict = json.loads( _r_args_dict )",
        "_r_args_vect = json.loads( _r_args_vect )",
        jython.command <- paste( "_r_return = ", py.foo, "(",
                                  ifelse( n.args.vect == 1, "_r_args_vect[0]", "*_r_args_vect" ),
                                  ifelse( n.args.dict == 0, ")", ", **_r_args_dict)" ), 
                                  sep = "" )
    )

    jython.exec( rJython, jython.command )
    rJython$exec( "_r_return = json.dumps( [ _r_return ] )" )
    ret <- fromJSON( .jstrVal( rJython$get ( "_r_return" )) )

    if( length( ret ) == 1 ) ret <- ret[[1]]

    ret
}

