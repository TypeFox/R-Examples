###
###
###		Samples
###
###
###


#
#	returns a string-vector of sample names used in the VCF file
#
vcf_getsamples <- function( vcffh )return( .Call("VCF_getSampleNames",vcffh,PACKAGE="WhopGenome") )


#
#	Selects from which samples data should be returned in the read_matrix functions
#
vcf_selectsamples <- function( vcffh, sampleslist )	return( .Call("VCF_selectSamples",vcffh,sampleslist,PACKAGE="WhopGenome") )

#
#
#
vcf_getselectedsamples <- function( vcffh ) return( .Call("VCF_getSelectedSamples",vcffh,PACKAGE="WhopGenome") )
