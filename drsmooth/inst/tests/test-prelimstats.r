#' @name test_prelimstats
#' @title Test Preliminary Statistics

# test_that("Preliminary Stats Correct", {
#     
#     # test scripts do not seem to recognize loaded libraries:
#     library(car)
#     library(clinfun)
#     library(multcomp)
#     library(pgirmess)
#     library(DTK)
#     library(mgcv)
#     library(segmented)
#     
#     x <- dosefactor("dose", DRdata)
#     tests=c("outlier", "bartlett", "shapiro", "chisquare", "jonckheere")
#     
#     # Perform tests.
# 	prelim_stats <- matrix(nrow=((ncol(x))*length(tests)), ncol=2)
# 	for (i in 1:length(tests)) {
#         
# 		test_name <- tests[i]
# 		f <- match.fun(test_name)
# 		output_matrix <- f(x)
# 		start_row <- (i-1)*(ncol(x))+1
# 		for (j in 1:nrow(output_matrix)) {
# 			prelim_stats[start_row+j-1,] <- output_matrix[j,]
# 		}
# 	}
#     expect_that (as.numeric(prelim_stats[2,2]), equals(7.9377479841e-05))
#     expect_that (as.numeric(prelim_stats[3,2]), equals(0.0617825962))
#     expect_that (as.numeric(prelim_stats[4,2]), equals(0.0619620098))
#     expect_that (as.numeric(prelim_stats[7,2]), equals(4.1290876979e-07))
#     expect_that (as.numeric(prelim_stats[8,2]), equals(.07267942091))
#     expect_that (as.numeric(prelim_stats[9,2]), equals(.43505618970))
#     expect_that (as.numeric(prelim_stats[12,2]), equals(1.032261968e-05))
#     expect_that (as.numeric(prelim_stats[13,2]), equals(0.424770007))
#     expect_that (as.numeric(prelim_stats[14,2]), equals(0.592202950))
#     expect_that (as.numeric(prelim_stats[17,2]), equals(1.158544858e-10))
#     expect_that (as.numeric(prelim_stats[18,2]), equals(0.138242036))
#     expect_that (as.numeric(prelim_stats[19,2]), equals(0.607051425))
#     expect_that (as.numeric(prelim_stats[22,2]), equals(1.961126039e-08))
#     expect_that (as.numeric(prelim_stats[23,2]), equals(1.961126039e-08))
#     expect_that (as.numeric(prelim_stats[24,2]), equals(1.961126039e-08))
#     }
# )