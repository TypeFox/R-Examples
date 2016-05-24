#######################################################################################################################

# Function: CreateDataStack.
# Argument: Data model and number of simulations.
# Description: Generates a data stack, which is a collection of individual data sets (one data set per simulation run).

CreateDataStack = function(data.model, n.sims, seed=NULL) {

  # Perform error checks for the data model and create an internal data structure
  data.structure = CreateDataStructure(data.model)

  # Check the seed if defined (the seed should be defined only when the user generate the data stack)
  if (!is.null(seed)){
    if (!is.numeric(seed))
      stop("Seed must be an integer.")
    if (length(seed) > 1)
      stop("Seed: Only one value must be specified.")
    if (nchar(as.character(seed)) > 10)
      stop("Length of seed must be inferior to 10.")
  }

  # Create short names for data model parameters
  outcome.dist = data.structure$outcome$outcome.dist
  outcome.type = data.structure$outcome$outcome.type
  outcome.dist.dim = data.structure$outcome$outcome.dist.dim
  data.sample.id = data.structure$id
  data.size = data.structure$sample.size.set
  data.event = data.structure$event.set
  rando.ratio = data.structure$rando.ratio
  data.design = data.structure$design.parameter.set
  data.outcome = data.structure$outcome.parameter.set

  # Number of outcome parameter sets, sample size sets and design parameter sets
  n.outcome.parameter.sets = length(data.structure$outcome.parameter.set)
  if (!is.null(data.structure$design.parameter.set)) {
    n.design.parameter.sets = length(data.structure$design.parameter.set)
  } else {
    n.design.parameter.sets = 1
  }

  # Determine if sample size or event were used
  # Determine which sample size set corresponds to the maximum of events or sample size for each data sample
  sample.size = any(!is.na(data.size))
  event = any(!is.na(data.event))
  if (sample.size) {
    n.sample.size.event.sets = dim(data.structure$sample.size.set)[1]
    max.sample.size = apply(data.size,2,max)
  } else if (event){
    n.sample.size.event.sets = dim(data.structure$event.set)[1]
    max.event = apply(data.event,2,max)
  }

  # Number of data samples specified in the data model
  n.data.samples = length(data.sample.id)

  # Create the data stack which is represented by a list of data sets (one data set for each simulation run)
  data.set = list()

  # Create a grid of the data scenario factors (outcome parameter, sample size and design parameter)
  data.scenario.grid = expand.grid(design.parameter.set = 1:n.design.parameter.sets,
                                   outcome.parameter.set = 1:n.outcome.parameter.sets,
                                   sample.size.set = 1:n.sample.size.event.sets)

  colnames(data.scenario.grid) = c("design.parameter", "outcome.parameter", "sample.size")

  # Number of data scenarios (number of unique combinations of the data scenario factors)
  n.data.scenarios = dim(data.scenario.grid)[1]

  # Create a grid of the ouctcome and design scenario factors (outcome parameter and design parameter)
  data.design.outcome.grid = expand.grid(design.parameter.set = 1:n.design.parameter.sets,
                                         outcome.parameter.set = 1:n.outcome.parameter.sets)

  colnames(data.design.outcome.grid) = c("design.parameter", "outcome.parameter")

  # Number of design and outcome scenarios (number of unique combinations of the design and outcome scenario factors)
  n.design.outcome.scenarios = dim(data.design.outcome.grid)[1]

  # Set the seed
  if (!is.null(seed)) set.seed(seed)

  # Loop over the simulations
  for (sim.index in 1:n.sims) {

    # If sample size is used (fixed number of sample size)
    if (sample.size) {

      design.outcome.variables = vector(n.design.outcome.scenarios, mode = "list")

      # Loop over the design and outcome grid
      for (design.outcome.index in 1:n.design.outcome.scenarios) {

        # Get the current design index and parameters
        current.design.index = data.design.outcome.grid[design.outcome.index, "design.parameter"]
        current.design.parameter = data.design[[current.design.index]]

        # Get the outcome index and parameters
        current.outcome.index = data.design.outcome.grid[design.outcome.index, "outcome.parameter"]
        current.outcome.parameter = data.outcome[[current.outcome.index]]

        # Initialized the data frame list
        df = vector(n.data.samples, mode = "list")

        # Loop over the data samples
        for (data.sample.index in 1:n.data.samples) {

          # Maximum sample size across the sample size sets for the current data sample
          current.max.sample.size = max.sample.size[data.sample.index]

          # Outcome parameter for the current data sample
          current.outcome = list(dist = outcome.dist, par = current.outcome.parameter[[data.sample.index]], type = outcome.type)

          # Get the current sample id
          current.sample.id = unlist(data.sample.id[[data.sample.index]])

          # Generate the data for the current design and outcome parameters
          df[[data.sample.index]] = GeneratePatients(current.design.parameter, current.outcome, current.sample.id, current.max.sample.size)

        } # Loop over the data samples

        design.outcome.variables[[design.outcome.index]] = list(design.parameter = current.design.index, outcome.parameter = current.outcome.index, sample = df)

      } # Loop over the design and outcome grid

      # Create the data scenario list (one element for each unique combination of the data scenario factors)
      data.scenario = list()

      # Loop over the data scenarios
      for (data.scenario.index in 1:n.data.scenarios) {

        design.index = data.scenario.grid[data.scenario.index, 1]
        outcome.index = data.scenario.grid[data.scenario.index, 2]
        sample.size.index = data.scenario.grid[data.scenario.index, 3]

        # Get the design.outcome variables corresponding to the current data scenario
        current.design.outcome.index = sapply(design.outcome.variables, function(x) x$design.parameter == design.index & x$outcome.parameter == outcome.index)
        current.design.outcome.variables = design.outcome.variables[current.design.outcome.index][[1]]$sample

        # Get the sample size
        current.sample.size = data.size[sample.size.index,]

        # Generate the data for the current data scenario
        data.scenario[[data.scenario.index]] = list(sample = CreateDataScenarioSampleSize(current.design.outcome.variables, current.sample.size))
      }

    } else if (event) {
      # If event is used (generate data until the number of event required for the first outcome is reached)

      design.outcome.variables = vector(n.design.outcome.scenarios, mode = "list")

      # Loop over the design and outcome grid
      for (design.outcome.index in 1:n.design.outcome.scenarios) {

        # Get the current design index and parameters
        current.design.index = data.design.outcome.grid[design.outcome.index, "design.parameter"]
        current.design.parameter = data.design[[current.design.index]]

        # Get the outcome index and parameters
        current.outcome.index = data.design.outcome.grid[design.outcome.index, "outcome.parameter"]
        current.outcome.parameter = data.outcome[[current.outcome.index]]

        # Initialized the data frame list
        df = vector(n.data.samples, mode = "list")

        # Initialized the temporary data frame list
        df.temp = vector(n.data.samples, mode = "list")

        # Set the Number of events
        n.observed.events = 0

        # Loop over the data samples to generate a first set of data corresponding to the maximum number of events required divided by randomization ratio
        for (data.sample.index in 1:n.data.samples) {

          # Outcome parameter for the current data sample
          current.outcome = list(dist = outcome.dist, par = current.outcome.parameter[[data.sample.index]], type = outcome.type)

          # Get the current sample id
          current.sample.id = unlist(data.sample.id[[data.sample.index]])

          # Generate the data for the current design and outcome parameters
          df.temp[[data.sample.index]] = GeneratePatients(current.design.parameter, current.outcome, current.sample.id, rando.ratio[data.sample.index] * ceiling(max.event / sum(rando.ratio)))

          # Merge the previous generated data with the temporary data
          if (!is.null(df[[data.sample.index]])) {
            data.temp = as.data.frame(mapply(rbind, lapply(df[[data.sample.index]], function(x) as.data.frame(x$data)), lapply(df.temp[[data.sample.index]], function(x) as.data.frame(x$data)), SIMPLIFY=FALSE))
            row.names(data.temp) = NULL
            df[[data.sample.index]] = lapply(df[[data.sample.index]], function(x) {return(list(id = x$id, outcome.type = x$outcome.type, data = as.matrix(data.temp)))})
          } else {
            df[[data.sample.index]] = df.temp[[data.sample.index]]
          }

        } # Loop over the data samples

        # Get the number of events observed accross all samples for the primary endpoint
        n.observed.events = sum(unlist(lapply(df, function(x) {return(!x[[1]]$data[,"patient.censor.indicator"])})))

        # Loop until the maximum number of events required is observed
        while(n.observed.events < max.event){

          # Loop over the data samples
          for (data.sample.index in 1:n.data.samples) {

            # Outcome parameter for the current data sample
            current.outcome = list(dist = outcome.dist, par = current.outcome.parameter[[data.sample.index]], type = outcome.type)

            # Get the current sample id
            current.sample.id = unlist(data.sample.id[[data.sample.index]])

            # Generate the data for the current design and outcome parameters
            df.temp[[data.sample.index]] = GeneratePatients(current.design.parameter, current.outcome, current.sample.id, rando.ratio[data.sample.index])

            # Merge the previous generated data with the temporary data
            if (!is.null(df[[data.sample.index]])) {
              data.temp = lapply(mapply(rbind, lapply(df[[data.sample.index]], function(x) as.data.frame(x$data)), lapply(df.temp[[data.sample.index]], function(x) as.data.frame(x$data)), SIMPLIFY=FALSE), function(x) as.matrix(x))
              df[[data.sample.index]] = mapply(function(x,y) {return(list(id=x$id, outcome.type = x$outcome.type, data = as.matrix(y, row.names = NULL)))}, x=df[[data.sample.index]], y=data.temp, SIMPLIFY=FALSE)
            } else {
              df[[data.sample.index]] = df.temp[[data.sample.index]]
            }

          } # Loop over the data samples

          # Get the number of events observed accross all samples for the primary endpoint
          n.observed.events = sum(unlist(lapply(df, function(x) {return(!x[[1]]$data[,"patient.censor.indicator"])})))

        } # Loop until the maximum number of events required is observed

        design.outcome.variables[[design.outcome.index]] = list(design.parameter = current.design.index, outcome.parameter = current.outcome.index, sample = df)

      } # Loop over the design and outcome grid

      # Create the data scenario list (one element for each unique combination of the data scenario factors)
      data.scenario = list()

      # Loop over the data scenarios
      for (data.scenario.index in 1:n.data.scenarios) {

        design.index = data.scenario.grid[data.scenario.index, 1]
        outcome.index = data.scenario.grid[data.scenario.index, 2]
        event.index = data.scenario.grid[data.scenario.index, 3]

        # Get the design.outcome variables corresponding to the current data scenario
        current.design.outcome.index = sapply(design.outcome.variables, function(x) x$design.parameter == design.index & x$outcome.parameter == outcome.index)
        current.design.outcome.variables = design.outcome.variables[current.design.outcome.index][[1]]$sample

        # Get the number of events
        current.events = data.event[event.index,]

        # Generate the data for the current data scenario
        data.scenario[[data.scenario.index]] = list(sample = CreateDataScenarioEvent(current.design.outcome.variables, current.events, rando.ratio))

      } # Loop over the data scenarios

    } # If event

    data.set[[sim.index]] = list(data.scenario = data.scenario)

  } # Loop over the simulations

  # Create the data stack
  data.stack = list(description = "data.stack",
                    data.set = data.set,
                    data.scenario.grid = data.scenario.grid,
                    data.structure = data.structure,
                    n.sims = n.sims,
                    seed = seed)

  class(data.stack) = "DataStack"
  return(data.stack)

}
# End of CreateDataStack
