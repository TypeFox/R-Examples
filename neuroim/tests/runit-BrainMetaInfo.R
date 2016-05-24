# TODO: Add comment
# 
# Author: brad
###############################################################################


test.asNIfTIMetaInfo <- function() {
	meta <- BrainMetaInfo(c(64,64,64), c(2,2,2), c(0,0,0), "DOUBLE")
	
	## as(meta, "NIfTIMetaInfo") requires a volume ....
	#nimeta <- as(meta, "NIfTIMetaInfo")
}
