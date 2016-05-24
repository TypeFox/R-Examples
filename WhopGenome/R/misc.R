###
###
###
###
###
###


vcf_getfieldnames <- function( vcff )	return( .Call("VCF_getFieldNames",vcff,PACKAGE="WhopGenome") )

vcf_getnumcontigs <- function( vcff )	return( .Call("VCF_getNumContigs",vcff,PACKAGE="WhopGenome") )

vcf_getcontignames <- function( vcff )	return( .Call("VCF_getContigNames",vcff,PACKAGE="WhopGenome") )

vcf_getheaderline <- function( vcff, whichnum )	return( .Call("VCF_getHeaderLine",vcff,whichnum,PACKAGE="WhopGenome") )


###
#
#	both use			f->getFieldPtr(REF) cmp f->getFieldPtr(ALT)
#

vcf_isSNP <- function( vcff )	return( .Call("VCF_isSNP",vcff,PACKAGE="WhopGenome") )

vcf_isINDEL <- function( vcff )	return( .Call("VCF_isInDel",vcff,PACKAGE="WhopGenome") )

