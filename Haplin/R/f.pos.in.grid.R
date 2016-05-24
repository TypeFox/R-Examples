f.pos.in.grid <- function(pos, A, comb){
#
# COMPUTES POSITION (ROW NUMBER) IN FULL GRID WITH length(A) COLUMNS, OR INVERSELY,
# COMPUTES GRID "COORDINATES" (COMBINATIONS OF COLUMN VALUES IN A ROW) FOR A GIVEN 
# POSITION (ROW NUMBER)
#
# pos = POSITION (numeric vector), A = vector WITH NUMBER OF ELEMENTS IN 
# EACH COLUMN, 
# comb = COMBINATION OF ELEMENTS IN A ROW (numeric matrix with length(A) columns)
#

	.ncol <- length(A)


	if(!missing(pos)){
	# FROM POSITION TO COMBINATION
		if(max(pos) > prod(A)) stop("Problem with position data")

		.comb <- matrix(NA, ncol = .ncol, nrow = length(pos))
		.postmp <- pos - 1
		
		for (i in seq(along = A)){
			.comb[,i] <- .postmp %% A[i]
			.postmp <- .postmp %/% A[i]
		}
		.comb <- .comb + 1
		return(.comb)
	}
#
#
	if(!missing(comb)){
	# FROM COMBINATION TO POSITION, CONVERTS SINGLE NUMERIC TO MATRIX IF NECESS.
		if(!is.matrix(comb)){
			if(is.numeric(comb) && (length(comb) == .ncol)){
				comb <- matrix(comb, ncol = .ncol)
			}else stop("Input is not a numeric matrix!")
		}
		if(any(apply(comb, 2, max) > A)) stop("Problem with position data")
		.weight <- c(1, cumprod(A[-.ncol]))
		.pos <- t(t(comb - 1) * .weight)
		.pos <- as.numeric(.pos %*% rep(1,.ncol)) + 1
		return(.pos)
	}



}
