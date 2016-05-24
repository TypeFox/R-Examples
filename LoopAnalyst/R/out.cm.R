# out.cm prints a community matrix or community effect matrix using the
# common format. It takes:
# M: a community matrix or community effect matrix

out.cm <- function(M) {
	
	N <- nrow(M)
	for (i in 1:N) {
		for (j in 1:N) {
			if (is.na(M[i,j])) {
				M[i,j] <- "?"
				}
			if (M[i,j] == 1) {
				M[i,j] <- "+"
				}
			if (M[i,j] == -1) {
				M[i,j] <- "-"
				}
			}
		}
	print(M,quote=FALSE)
	
	# end out.cm()
	}
