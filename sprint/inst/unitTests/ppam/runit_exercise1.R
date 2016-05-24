# SPRINT: Parallel computing with R on HECToR
# 01-Dec-11
# Exercise 1: pcor() and ppam()
#library("RUnit")
#library(cluster)
# Combined Test and Training Sets from the Golub Paper
#library(golubEsets)

#==============================================================================
#                        Load Data
#==============================================================================


# The data are from Golub et al. These are the combined training samples and 
# test samples. There are 47 patients with acute lymphoblastic leukemia (ALL) 
# and 25 patients with acute myeloid leukemia (AML). The samples were assayed 
# using Affymetrix Hgu6800 chips and data on the expression of 7129 genes 
# (Affymetrix probes) are available.

data(Golub_Merge)
data <- exprs(Golub_Merge)

#==============================================================================
#                        SPRINT 
#==============================================================================
#library(sprint)
#library(ff)

# test pcor and ppam together
test.pcor_and_ppam <- function(){
	
stime <- proc.time()["elapsed"]
# Calculate distance matrix using correlation function. The parallel version writes
# its output to a file that is loaded as an ff object in R and behaves (almost) 
# as if data was stored in RAM
tdata <- t(data)
distance_matrix <- cor(tdata)
pdistance_matrix <- pcor(tdata)

class(distance_matrix)
class(pdistance_matrix)
	
	invisible(checkEqualsNumeric(distance_matrix, pdistance_matrix[,]))


# Force memory clean-up
gc(reset=TRUE, verbose=FALSE)


# Find 6 medoids using the PAM algortithm
# You are passing a distance matrix on input, so set the diss parameter to TRUE
# The parallel version of the algorithm will automatically detect the format of
# the input data so no additional parameter is required
pam_result <- pam(distance_matrix, k=6, diss=TRUE)
ppam_result <- ppam(pdistance_matrix, k=6)
checkEquals(pam_result$medoids,
ppam_result$medoids, " Medoids should have same names.")

}
# do further analysis of data then exit MPI
# ...
#
# You can always save your data or objects on disk and continue to work on it
# in your interactive R window

# Shutdown the SPRINT library
# pterminate()

#quit(save="no")

