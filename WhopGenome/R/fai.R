
#
#
#
fai_open <- function( filename ) .Call("FAI_open", filename, PACKAGE="WhopGenome" )
fai_close <- function( faifh ) .Call("FAI_close", faifh, PACKAGE="WhopGenome")
fai_reopen <- function( faifh ) .Call("FAI_reopen", faifh, PACKAGE="WhopGenome")
fai_build <- function( filename ) .Call("FAI_build", filename, PACKAGE="WhopGenome")
fai_query3 <- function( faifh, regionstring, resultstring ) .Call("FAI_query3", faifh, regionstring, resultstring, PACKAGE="WhopGenome")
fai_query5 <- function( faifh, sequencename, beginpos, endpos, resultstring ) .Call("FAI_query5", faifh, sequencename, beginpos, endpos, resultstring, PACKAGE="WhopGenome")
