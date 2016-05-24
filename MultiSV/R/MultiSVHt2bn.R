##' @export

###############################################################################
#
# MultiSVExe
#
# Compute Bin Counts from hits file
#
# Copyright 2014, Khurram Maqbool <khurram.maqbool@outlook.com>
#
# This software is distributed under the terms of the GNU General
# Public License Version 2, June 1991.  


`ComputeBinCounts.default` <- function(RDBinSampleFile, RDBinChrSizeFile, RDBinSampleName , 
								RDBinWindowSize, OutFolder)
{
Exec <- gtExEs("Hits2MultiSVBin")
RDBinOut <- paste(OutFolder,"BinRDCount.",RDBinWindowSize,"/", sep="")
system(paste(Exec,
               " --Sample ", RDBinSampleFile,
               " --ChrFile ", RDBinChrSizeFile,
               " --samplename ", RDBinSampleFile,
               "  --window-size ", RDBinWindowSize,
               " --outputdirectory ", RDBinOut,
               sep = "")
)
}
