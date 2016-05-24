MCp.RInSp = function(dataset, pop.diet = "sum", replicates = 999){
  #
  # The procedure replicates the two strategies of Monte Carlo resampling available in
  # the following commands: Emc, PSicalc, WTdMC, and WTcMC.
  #
  # Author: Nicola ZACCARELLI
  # E-mail: nicola.zaccarelli@gmail.com
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
  if (class(dataset) != "RInSp") stop("The input must be an object of class RInSp.")
  if (dataset$data.type %in% c("integer", "double") == FALSE) stop("Input data type must be 'integer' or 'double'.")
  if (pop.diet %in% c("sum", "average") == FALSE) stop("Wrong population diet option.")
  replicates = as.integer(replicates)
  if (replicates <=0) stop("Wrong value for replicates.")
  if (dataset$data.type == "integer") type.mc = 1 else type.mc= 2
  if (pop.diet == "sum") {
    diet = apply(dataset$resources, 2, sum) / sum(dataset$resources)} else {
    diet = apply(dataset$proportion, 2, mean)}
  ris = .Call("MCprocedure", dataset$resources, as.vector(type.mc), as.vector(diet), as.vector(replicates), PACKAGE="RInSp")
  ris = array(c(as.vector(dataset$resources), as.vector(ris)), c(dataset$num.individuals, dataset$num.prey, (replicates +1)))
 return(ris)}
