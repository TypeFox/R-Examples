# require(testthat)

context("heights_per_k.dendrogram works")


test_that("heights_per_k.dendrogram works",{
   
   dend <- as.dendrogram(hclust(dist(USArrests)))

   # '    dendextendRcpp:::heights_per_k.dendrogram(dend),
   # '    dendextend:::heights_per_k.dendrogram(dend)
   
#    require(dendextend)
   expect_identical(dendextendRcpp::dendextendRcpp_heights_per_k.dendrogram(dend),
                    old_heights_per_k.dendrogram(dend))
   

   
   
#    require(microbenchmark)
#    microbenchmark(
#       dendextendRcpp:::get_branches_heights(dend),
#       dendextend:::get_branches_heights(dend)
#    )
#    # Rcpp is 50 times faster...
#    
   
})


