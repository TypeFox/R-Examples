"sim" <- 
function(x, data, coords, grid, method = "ik", ...) {
  # Generation of conditional simulation based on user specified method
  #
  #        x a multi.tpfit object
  #     data vector of data
  #   coords coordinates matrix
  #     grid simulation points
  #   method method to perform prediction and simulation c("ik", "ck", "path", "mcs")
  #      ... further option to pass to the function sim_*

  # Further arguments for Indicator Kriging (ik) and coKriging (ck)
  #      knn number of k-nearest neighbours
  # ordinary boolean (if TRUE ordinary Kriging is applied rather than simple Kriging)
  #       GA boolean (if TRUE genetic algorithm is applied rather than simulated annealing)
  #   optype character with the objective function to minimize after the simulation
  #   max.it maximum number of iteration for the optimization method

  # Further arguments for Fixed and Random Path methods (path)
  # radius radius to find neighbour points
  #  fixed boolean for random or fixed path algorithm

  # Further arguments for Multinomial Categorical Simulation (mcs)
  #    knn number of k-nearest neighbours (if NULL all data are neighbours)
  # radius radius to find neighbour points

  if (method == "ck") return(sim_ck(x, data, coords, grid, ...))
  if (method == "path") return(sim_path(x, data, coords, grid, ...))
  if (method == "mcs") return(sim_mcs(x, data, coords, grid, ...))
  if (method != "ik") warning("Simulation method not recognized. Indicator Kriging method (\"ik\") set by default.")
  return(sim_ik(x, data, coords, grid, ...))
}
