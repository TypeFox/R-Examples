# cem.corr() produces perturbation correlation tables from a community effect
# matrix. It takes:
# CEM: a community effect matrix


# check for version compatibility and notify user of version incompatibility
# and let them know i am ammenable to making back-compatible revisions.
cem.corr <- function(CEM) {


# out.cem.corr prints the correlation matrices using the format described in
# Puccia and Levins. It takes:
# M: a list of length N of correclation matrices of size N.

	out.cem.corr <- function(M) {

		N <- length(M)
		for (corr in 1:N) {
			for (i in 1:N) {
				for (j in 1:N) {
					if (is.na(M[[corr]][i,j])) {
						M[[corr]][i,j] <- "?"
						}
					if (M[[corr]][i,j] == 1) {
						M[[corr]][i,j] <- "+"
						}
					if (M[[corr]][i,j] == -1) {
						M[[corr]][i,j] <- "-"
						}
					if (!(j >= i)) {
						M[[corr]][i,j] <- " "
						}
					if ((i == j) & (M[[corr]][i,j] != "0")) {
						M[[corr]][i,j] <- 1
						}
					}
				}
			cat("\nInput to:", parameter.names[corr],"\n")
			print(M[[corr]],quote=FALSE)
			}
	
		# end out.cem.corr()
		}

	N <- nrow(CEM)
	parameter.names <- rownames(CEM)

# validate that the supplied matrix is a community effect matrix in several
# steps

	# validate that the matrix is square
	if (N != ncol(CEM)) {
		stop("\nsupplied matrix is not square; community effect matrix expected!\n")
		}

	# validate that the matrix contains elements that are one of 1, 0, -1 or NA
	for (i in 1:N) {
		for (j in 1:N) {
			a.ij <- CEM[i,j]
			if ( !( (1 == abs(a.ij)) | (0 == a.ij) | (is.na(a.ij)) ) ) {
				stop("\nSupplied matrix has at least one invalid element (i.e. \nnot 1, 0, -1 or NA); community effect matrix expected!\n")
				}
			}
		}
 

# Pith of cemcorr


	# create set of matrices	
	CorrelationMatrices <- rep(list(matrix(c(NA),N,N,dimnames=list(parameter.names,parameter.names))),N)
	for (corr in 1:N) {
		# populate upper diagonal
		for (i in 1:N) {
			for (j in 1:N) {
				if (j >= i) {
					CorrelationMatrices[[corr]][i,j] <- CEM[corr,i]*CEM[corr,j]
					}
				}
			}
		}


	return(out.cem.corr(CorrelationMatrices))
	}
