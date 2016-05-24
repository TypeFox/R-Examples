sessionInfo()
data(Golub_Merge)
data <- exprs(Golub_Merge)

#==============================================================================
#                        SPRINT 
#==============================================================================
#library(sprint)
#library(ff)

# test pcor and ppam together
test.ppam_loop_test <- function(){

stime <- proc.time()["elapsed"]
# Calculate distance matrix using correlation function. The parallel version writes
# its output to a file that is loaded as an ff object in R and behaves (almost) 
# as if data was stored in RAM
tdata <- t(data)
   sink(file = "/Users/egrant1/Documents/SPRINT/workspace/sprint/trunk/test_suite/ppam_loop_results", append = TRUE, type = c("output", "message"),
                 split = FALSE)

for(i in 1:10) {

        print("i is: ")
        print(paste(i))

        distance_matrix <- cor(tdata)
        pam_result <- pam(distance_matrix, k=6, diss=TRUE)

        pdistance_matrix <- pcor(tdata)

        invisible(checkEqualsNumeric(distance_matrix, pdistance_matrix[,]))


# Force memory clean-up
        gc(reset=TRUE, verbose=FALSE)


# Find 6 medoids using the PAM algortithm
# You are passing a distance matrix on input, so set the diss parameter to TRUE
# The parallel version of the algorithm will automatically detect the format of
# the input data so no additional parameter is required
        ppam_result <- ppam(pdistance_matrix, k=6)

	print(" medoids are ")
	print(paste(ppam_result$medoids))   
     checkEquals(pam_result$medoids,
                                ppam_result$medoids, " Medoids should have same names.")

        }
}

