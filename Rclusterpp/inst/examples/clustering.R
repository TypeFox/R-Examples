#!/usr/bin/env Rscript

library(Rclusterpp)
suppressMessages(require(inline))

basic_clustering <- function(data) {
  fx <- cxxfunction( signature(data = "matrix"), '
    using namespace Rclusterpp;
    using namespace Eigen;
  
    // Convert from R data structure to Eigen without any data copying (the "Map" part)
    MapNumericMatrix data_e(as<MapNumericMatrix>(data));
    
    // Compute the distance matrix. Note this can be bettter vectorized. Also
    // note we create a strictly lower matrix.
    Eigen::NumericMatrix data_d = Eigen::NumericMatrix::Zero(data_e.rows(), data_e.rows());
    for (ssize_t i=0; i<(data_e.rows()-1); i++) {
      for (ssize_t j=i+1; j<data_e.rows(); j++) {
        data_d(j, i) = norm( data_e.row(i) - data_e.row(j) );  // Euclidean distance
      }
    }
    StrictlyLowerNumericMatrix data_t = data_d.triangularView<Eigen::StrictlyLower>(); 

    typedef NumericCluster::plain cluster_type;
    ClusterVector<cluster_type> clusters(data_t.rows());  // Create cluster vector

    // Perform average link clustering via recursive nearest neighbor method
    cluster_via_rnn( 
      average_link<cluster_type>(data_t, FromDistance), 
      init_clusters(data_t, clusters) 
    );

    // Return merge, height, etc. to R...
    return wrap( clusters );
  ',
  plugin = "Rclusterpp", verbose=FALSE)
  return(fx(data))
}

r <- basic_clustering(as.matrix(USArrests))
h <- hclust(dist(USArrests, method="euclidean"), method="average")

if (identical(r$merge, h$merge)) {
  message("Success! Clustering output is the same between Rclusterpp and stats::hclust ...")
} else {
  error("Failure! Clustering output doesn't match!")
}
