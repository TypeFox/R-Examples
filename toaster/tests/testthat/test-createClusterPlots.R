context("create Cluster Plots")

test_that("Cluster plot functions throw errors", {
  
  expect_error(createCentroidPlot(km=structure(list(centers=matrix(c(1), nrow=1, byrow = TRUE)),
                                                    class = c("toakmeans", "kmeans")), 
                                  format="no-such-vis-format"),
               "'arg' should be one of \"line\", \"bar\", \"heatmap\", \"bar_dodge\"")
  
  expect_error(createCentroidPlot(),
               "Kmeans object must be specified.")
  
  expect_error(createCentroidPlot(km=structure(list(), class = c("toakmeans", "kmeans"))),
               "Kmeans object is missing cluster centers.")
  
  expect_error(createClusterPlot(),
               "Kmeans object must be specified.")
  
  expect_error(createClusterPlot(km=structure(list(), class = c("toakmeans", "kmeans"))),
               "Kmeans object is missing cluster aggregates.")
  
  expect_error(createClusterPairsPlot(),
               "Kmeans object must be specified.")
  
  expect_error(createClusterPairsPlot(km=structure(list(), class = c("toakmeans", "kmeans"))),
               "Kmeans object is missing sample data.")
  
  expect_error(createSilhouetteProfile(),
               "Kmeans object must be specified.")
  
  expect_error(createSilhouetteProfile(km=structure(list(), class = c("toakmeans", "kmeans"))),
               "Kmeans object is missing silhouette data.")
  
  expect_error(createSilhouetteProfile(km=structure(list(sil=list()), class = c("toakmeans", "kmeans"))),
               "Kmeans object is missing silhouette data.")
})