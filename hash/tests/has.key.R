
library( hash ) 

#' Ensure that there are no keys on hash
h <- hash()
all( has.key( 'c' , h ) == FALSE )

