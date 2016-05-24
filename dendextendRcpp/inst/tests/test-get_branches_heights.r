# require(testthat)

context("get_branches_heights works")


test_that("get_branches_heights works",{
   
   dend <- as.dendrogram(hclust(dist(USArrests)))

   # require(dendextend) # We Depend on this package, so no need to load it.
   # expect_identical(dendextendRcpp:::get_branches_heights(dend),
                    # dendextend:::get_branches_heights(dend))
   # The above will not work, since we mask this package inside the NAMESPACE...
   expect_identical(dendextendRcpp::dendextendRcpp_get_branches_heights(dend),
                    old_get_branches_heights(dend))

   
   
#    require(microbenchmark)
#    microbenchmark(
#       dendextendRcpp:::get_branches_heights(dend),
#       dendextend:::get_branches_heights(dend)
#    )
#    # Rcpp is 50 times faster...
#    
   
})


