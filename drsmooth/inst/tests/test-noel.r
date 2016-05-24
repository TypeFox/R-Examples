#' @name test_noel
#' @title Test NOEL Statistics

# test_that("NOEL Stats Correct", {
#     
#     # test scripts do not seem to recognize loaded libraries:
# #     library(car)
# #     library(clinfun)
# #     library(multcomp)
# #     library(pgirmess)
# #     library(DTK)
# #     library(mgcv)
# #     library(segmented)
# 
#     data <- dosefactor("dose", DRdata)
#     # tests = c("dunnetts1", "dunnetts2", "dunns1", "dunns2", "dunnettst3")
#     direction = "greater"
#     alpha = 0.05
#     
# #     dunnetts1_matrix <- dunnetts1("MF_Log", alpha, direction, data)
# #     dunnetts1_values <- dunnetts1_matrix[[2]]
# #     
# #     dunnetts2_matrix <- dunnetts2("MF_Log", alpha, direction, data)
# #     dunnetts2_values <- dunnetts2_matrix[[2]]
# #     
# #     dunns1_matrix <- dunns1("MF_Log", alpha, direction, data)
# #     dunns1_values <- dunns1_matrix[[2]]
# #     
# #     dunns2_matrix <- dunns2("MF_Log", alpha, direction, data)
# #     dunns2_values <- dunns2_matrix[[2]]
#         
#     dunnettst3_matrix <- dunnettst3("MF_Log", alpha, direction, data)
#     dunnettst3_values <- dunnettst3_matrix[[1]]
#     
# (targetcolumn, alternative, alpha, control, tot.obs, label, data) 
# 
#     # Dunnetts tests return different p-values every time.
#     # expect_that (as.numeric(dunnetts1_values[2,5]), equals(.9866647952))
#     # expect_that (as.numeric(dunnetts2_values[2,5]), equals(.9775326063))
#     # expect_that (dunns1_values[1,1], equals(8))
#     # expect_that (dunns2_values[9,1], equals(25.9))
#     expect_that (as.numeric(dunnettst3_values[7,3]), equals(-.1661))
#     }
# )