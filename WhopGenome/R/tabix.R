###
###
###		Access to any-format tabix-indexed files
###
###
###
###
###


##
##	Tabix : file
##
tabix_open <- function( filename )	.Call("tabix_open",filename,PACKAGE="WhopGenome" )
tabix_close <- function( tabfh )	.Call("tabix_close",tabfh,PACKAGE="WhopGenome" )
tabix_reopen <- function( tabfh )	.Call("tabix_reopen",tabfh,PACKAGE="WhopGenome" )

##
##	Tabix : regions
##
tabix_getregion <- function( tabfh )	.Call("tabix_getRegion",tabfh,PACKAGE="WhopGenome" )
tabix_restartregion <- function( tabfh )	.Call("tabix_restartRegion",tabfh,PACKAGE="WhopGenome" )
tabix_setregion <- function( tabfh, tid, beginpos, endpos ) .Call("tabix_setRegion",tabfh,tid,beginpos,endpos,PACKAGE="WhopGenome" )


##
##	Tabix : read
##
tabix_readraw <- function( tabfh ) return( .Call("tabix_readLine",tabfh,PACKAGE="WhopGenome") )
tabix_read <- function( tabfh )
{
	ln <- .Call("tabix_readLine",tabfh,PACKAGE="WhopGenome")
	ln <- strsplit(ln,"\t")[[1]]
	flds <- strsplit(ln[9],";")[[1]]
	flds <- unlist( strsplit(flds,"=") )
	return(ln)
}

##
##
##


##
##
##
tabix_build <- function( filename, sc, bc, ec, meta, lineskip ) .Call("tabix_build", filename, sc, bc, ec, meta, lineskip, PACKAGE="WhopGenome")
