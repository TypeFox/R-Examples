hellingerSupplementaryColsPreProcessing <- function(SUP.DATA,W=NULL,M=NULL){


	hell.preproc.again.all <- makeRowProfiles(SUP.DATA,hellinger=TRUE,weights=W,masses=M)
	t.hell.preproc <- t(hell.preproc.again.all$deviations)
	# if(is.null(W)){
		# W <- hell.preproc.again.all$weights
	# }else if(length(W)!=ncol(SUP.DATA)){
		# print('Length of W does not match column dim of SUP.DATA. Using default.')
		# W <- hell.preproc.again.all$weights
	# }
	# if(is.null(M)){
		# M <- hell.preproc.again.all$masses
	# }else if(length(M)!=nrow(SUP.DATA)){
		# print('Length of M does not match row dim of SUP.DATA. Using default.')
		# M <- hell.preproc.again.all$masses
	# }
	t.hell.preproc <- apply(t.hell.preproc,2,'*',hell.preproc.again.all$weights)
	t.hell.preproc <- t(apply(t.hell.preproc,1,'*',hell.preproc.again.all$masses)) ##this can be optional; masses in or derived from.	
	
}