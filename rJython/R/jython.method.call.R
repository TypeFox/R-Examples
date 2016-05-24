#########################################################
# CGB, 20100703
#########################################################

jython.method.call <- function( rJython, py.object, py.method, ... ){
    jython.call( rJython, paste( py.object, py.method, sep = "." ) )
}

