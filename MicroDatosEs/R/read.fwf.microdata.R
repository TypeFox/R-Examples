###################################################################
# cjgb
# 20120305
# Reads a generic microdata file
###################################################################

read.fwf.microdata <- function( file, file.mdat.1, file.mdat.2, fileEncoding = "UTF-8" ){

	## read fixed file using mdat1 metadata file

	mdat.1 <- read.table( file.mdat.1, header = T, sep = "\t", fileEncoding = fileEncoding )
	dat <- read.fwf( file, mdat.1$width, header = F, col.names = mdat.1$var, fileEncoding = fileEncoding )

	rm( mdat.1 ); gc()

	# Replaces keys in raw data by actual column values

	mdat.2 <- read.table( file.mdat.2, header = T, sep = "\t", fileEncoding = fileEncoding )

	assign.labels <-  function( v, var.name, mdat ){
		tmp <- mdat[ mdat$var == var.name, ]

		if( nrow( tmp ) == 1 ){
			if( ! is.na( tmp$nulo ) && any( v == tmp$nulo, na.rm = T ) )
				v[ v == tmp$nulo ] <- NA
			if( tmp$tipo == "HHMM" )
				v <- v %/% 100 + ( v %% 100 ) / 60
			return( v )
		}

		indices <- match( v, tmp$llave )
		return( as.character( tmp$valor )[ indices ] )
	}

	as.data.frame( sapply( names( dat ), function( x ) assign.labels( dat[[x]], x, mdat.2 ), simplify = F ) )
}

