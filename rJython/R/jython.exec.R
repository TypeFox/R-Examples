#########################################################
# CGB, 20100703
#########################################################

jython.exec <- function( rJython, python.code ){

    rJython$exec( "_r_error = None" )

    # Creating the call
    if( length( python.code ) == 1 )
        python.code <- strsplit(python.code, "\n|\r\n" ) 

    python.code <- as.character( sapply( python.code, function( x ) paste( "\t", x, sep = "" ) ) )

    python.code <- paste( "try:", paste( python.code, collapse = "\n" ), "except Exception, e:_r_error = e.__str__()\n", sep = "\n" )

    rJython$exec( python.code )
    rJython$exec( "_r_error = json.dumps( _r_error )" )

    python.exception <- fromJSON( .jstrVal( rJython$get ( "_r_error" ) ) ) 

    if( !is.null( python.exception ) )
        stop( python.exception )

    invisible( NULL )
}



