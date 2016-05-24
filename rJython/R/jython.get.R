#########################################################
# CGB, 20100708, created
#########################################################

jython.get <- function( rJython, py.var ){

    jython.command <- paste( "_r_return = json.dumps( [", py.var, "] )", sep = "" )
    jython.exec( rJython, jython.command )
    ret <- fromJSON( .jstrVal( rJython$get ( "_r_return" )) )
    if( length( ret ) == 1 ) ret <- ret[[1]]
    ret
}

