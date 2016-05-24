######################################################################################
# Created 23Jan12: This file contains the functions needed to calculate H and F
######################################################################################

###########################################################################################
#
#							Calculate H
#
###########################################################################################
############ Simple and proportional case ###########################
calc_H = function(mat01, H00, Rmp, M, inBinMat){
	
	# Calculate ds using the matrices on p328 of the SINICA paper.
	Rmpvec = rep(NA, length = sum(2^(1:M)))
	Rmpvec[1:(sum(2^(1:M))/2)*2-1] = Rmp
	Rmpvec[1:(sum(2^(1:M))/2)*2] = 1-Rmp
	ds = H00*exp(mat01%*%log(Rmpvec))
	
	# Use inBinMat, which is a matrix that has a '1' for every bin the subject survived
	# through, and a ratio for the bin they failed in, and multiply that by the ds.
	# Take the row sum for each subject to get Ht.
	Ht = rowSums(inBinMat*matrix(ds, ncol = length(ds), nrow = nrow(inBinMat), byrow = TRUE))
	
	return(Ht)
}

###########################################################################################
#
#							Calculate F
#
###########################################################################################
################# Simple and proportional case ###########################
calc_F = function(mat01, Rmp, M, inBinMat){
	
	Rmpvec = rep(NA, length = sum(2^(1:M)))
	Rmpvec[1:(sum(2^(1:M))/2)*2-1] = Rmp
	Rmpvec[1:(sum(2^(1:M))/2)*2] = 1-Rmp
	
	# Splits holds the selected Rmp values for each bin. For example bin1 = R1,0*R2,0*R3,0 
	splits = exp(mat01%*%log(Rmpvec))
	
	F = rowSums(inBinMat*matrix(splits, ncol = length(splits), nrow = nrow(inBinMat), byrow = TRUE))
	
	# Divide the current F, which is currently only the numerator, but the sum of all the splits, 
	# which is the denominator of F.
	F = F/sum(splits)
	
	return(F)
}


