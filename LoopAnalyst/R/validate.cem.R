# validate.CEM validates a community effect matrix by testing that it is
# square, has only elements of NA, -1, 0 or 1, and has at least two
# parameters. It produces an error if the community effect matrix is 
# invalid by these tests, and returns nothing otherwise. It takes:
# CEM: a potential community matrix
validate.cem <- function(CEM) {

	# Is CEM a matrix? a square matrix?
 	if (!(is.matrix(CEM)) | !(identical( nrow(CEM), ncol(CEM) ) ) ) {
	 	stop("\nA Community Effect Matrix must be a square matrix with elements \nof values of only NA, 1, 0 and -1.")
	 	}

	# Is CEM big enough?
	if (nrow(CEM) == 1) {
	 	stop("\nA Community Effect Matrix must have two or more parameters.")
	 	}

	# Does CEM contain only values = 1, 0 or -1?
	for (i in 1:nrow(CEM)) {
		for (j in 1:ncol(CEM)) {
			if ( !is.na(CEM[i,j]) ) {
	 			if ( !( (CEM[i,j] == 1) | (CEM[i,j] == 0) | (CEM[i,j] == -1) ) ) {
	 				stop("\nA Community Effect Matrix must be a square matrix with elements \nof values of only NA, 1, 0 and -1.")
	 				}
	 			}
	 		}
	 	}

	# end validate.CEM
	}
