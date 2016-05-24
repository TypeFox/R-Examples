
#test_lazyScore <- function(){
#  pfm <- matrix(c(3, 5, 4, 2, 7, 0, 3, 4, 9, 1, 1, 3, 3, 6, 4, 1, 11, 
#                  0, 3, 0, 11, 0, 2, 1, 11, 0, 2, 1, 3, 3, 2, 6, 4, 1, 
#                  8, 1, 3, 4, 6, 1, 8, 5, 1, 0, 8, 1, 4, 1, 9, 0, 2, 3, 
#                  9, 5, 0, 0, 11, 0, 3, 0, 2, 7, 0, 5),
#                nrow = 4, dimnames = list(c("A","C","G","T"))
#                )
#  bg <- c(A=0.25, C=0.25, G=0.25, T=0.25)
#  pvalue <- 1e-5
#  type <- "PFM"
#  granularity <- 1e-5
#  #dyn.load("/Users/gtan/Repositories/Bitbucket/TFMPvalue/src/TFMPvalue.so")
#  score <- TFMLazyScore(pfm, pvalue, bg, type, granularity)
#  checkEqualsNumeric(8.77361, score, tolerance=0.01)
#
#}

