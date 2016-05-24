heatmap.samples <-
function(samples) {
  ## Provide a density map for the probability distribution of shift vectors
  ## Only works for two dimensions
  ## assuming there are two dimensions of input data
  ## produce a density map of the vectors in "samples"
  dens = kde2d(samples[,1], samples[,2]);
  filled.contour(dens);
}
