#########################################################
# CGB, 20100707, created
#########################################################

jython.assign <- function( rJython, var.name, value ){

    value <- toJSON( value )

    # Creating the call

    jython.command <- c( 
        paste( var.name , "='", value, "'",  sep = "" ),
        paste( var.name , "= json.loads(", var.name, ")", sep = "" )
    )

    jython.exec( rJython, jython.command )
    invisible( NULL )
}

