######################################################################################################################

# Function: SimParameters
# Argument: Multiple character strings.
# Description: This function is called by default.
#' @export
SimParameters = function(n.sims, seed, proc.load = 1) {

  # Error checks
  if (!is.numeric(n.sims)) stop("SimParameters: Number of simulation runs must be an integer.")
  if (length(n.sims) > 1) stop("SimParameters: Number of simulations runs: Only one value must be specified.")
  if (n.sims%%1 != 0) stop("SimParameters: Number of simulations runs must be an integer.")
  if (n.sims <= 0) stop("SimParameters: Number of simulations runs must be positive.")

  if (!is.numeric(seed)) stop("Seed must be an integer.")
  if (length(seed) > 1) stop("Seed: Only one value must be specified.")
  if (nchar(as.character(seed)) > 10) stop("Length of seed must be inferior to 10.")

  if (is.numeric(proc.load)){
    if (length(proc.load) > 1) stop("SimParameters: Processor load only one value must be specified.")
    if (proc.load %%1 != 0) stop("SimParameters: Processor load must be an integer.")
    if (proc.load <= 0) stop("SimParameters: Processor load must be positive.")
  }
  else if (is.character(proc.load)){
    if (!(proc.load %in% c("low", "med", "high", "full"))) stop("SimParameters: Processor load not valid")
  }

  sim.parameters = list(n.sims = n.sims,
                        seed = seed,
                        proc.load = proc.load)

  class(sim.parameters) = "SimParameters"
  return(sim.parameters)
  invisible(sim.parameters)
}