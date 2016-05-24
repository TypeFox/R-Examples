###
###
###		VCF file open/close/valid
###
###
###

######
#
#	is given variable a VCFhandle ?
#
vcf_valid <- function( vcffh )
{
	if( class(vcffh) != "externalptr" ){ return( FALSE );}
	if( names(attributes(vcffh)[1]) != "VCF.filename" ){ return( FALSE );}
	return( TRUE );
}


######
#
#	Open a VCF file
#
vcf_open <- function( filename )
{
	if( ! is.character(filename) || length(filename)!=1 )
		stop("vcf_open : parameter <filename> needs to be a single string!");
	res = .Call("VCF_open", filename,PACKAGE="WhopGenome" )
	return( res )
}



######
#
#	Close a VCF file
#
vcf_close <- function( vcf_filehandle )
{
	if( class(vcf_filehandle) != "externalptr" )
		stop("vcf_close : parameter <vcffh> must be a VCFhandle as returned by vcf_open!");
	res = .Call("VCF_close", vcf_filehandle,PACKAGE="WhopGenome" )
	return( res )
}


######
#
#	Open a VCF file
#
vcf_reopen <- function( vcffh )
{
	if( class(vcffh) != "externalptr" )
		stop("vcf_reopen : parameter <vcffh> must be a VCFhandle as returned by vcf_open!");
	res = .Call("VCF_reopen", vcffh,PACKAGE="WhopGenome" )
	return( res )
}

######
#
#	Build a tabix index for a (already compressed) VCF file
#
vcf_buildindex <- function( filename ) .Call("VCF_buildindex", filename,PACKAGE="WhopGenome")

