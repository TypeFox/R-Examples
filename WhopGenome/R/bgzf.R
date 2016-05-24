
#
#
#
bgzf_compress <- function( infilename, outfilename ) .Call("BGZF_compress", infilename, outfilename, PACKAGE="WhopGenome")
