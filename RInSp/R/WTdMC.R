WTdMC = function(dataset, pop.diet = "sum", replicates=999, print.ris=TRUE){
  #
  # The program calculates the Total Niche Width (TNW), and breaks TNW down into its Between
  # Individual Component (BIC) and Within Individual Component (WIC).
  # Discrete data types are used.
  #
  # Author: Nicola ZACCARELLI, Giorgio MANCINELLI, Dan BOLNICK
  # E-mail: nicola.zaccarelli@gmail.com,
  #         giorgio.mancinelli@unisalento.it
  #         danbolnick@mail.texas.edu
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
  if (class(dataset) != "RInSp")
    stop("The input must be an object of class RInSp.")
  if(dataset$data.type != "integer") stop("Input data type must be integer.")
  if (pop.diet %in% c("sum", "average") == FALSE) stop("The specified population diet type is wrong.")
  if (pop.diet == "sum") diet.pop = 0 else diet.pop = 1
  replicates = as.integer(replicates)
  if (replicates <=1) stop("Wrong number of replicates.")
  if (print.ris == TRUE) cat("\n If your dataset is big, this can take time. Please be patient. \n")
# coerce vec to be double
  if (!is.double(dataset$resources)) dataset$resources = matrix(as.double(dataset$resources), dataset$num.individuals, dataset$num.prey)
  if(!is.integer(replicates)) replicates = abs(as.integer(replicates))
  Ris = .Call("WTdMC", dataset$resources, as.vector(diet.pop), as.vector(replicates), PACKAGE="RInSp")
  attributes(Ris)$dimnames[[2]] = c("Zero", "WIC", "BIC", "TNW", "WonT")
  cum.distr = ecdf(Ris[, 5])
  pvalue= cum.distr(Ris[1, 5])
# Calculate list of individuals with Shannon-Weaver value of zero
  checkZero = (dataset$proportion == 1)
  if (sum(checkZero) > 0) Zeros = dataset$ind.names[rowSums(checkZero)*c(1:length(dataset$ind.names))] else Zeros = Ris[1,1]
# Build list object for output  
  Ris2= list(WonT= Ris[1, 5], Zeros = Zeros, p.value= cum.distr(Ris[1, 5]), montecarlo= Ris, parameter = 5)
  class(Ris2) = "RInSp"
  if (print.ris == TRUE){
 cat("\n Using Roughgarden's 1979 equations, based on Shannon-Weaver diversity index: ")
 cat("\n Within-individual component          = ", Ris[1,2])
 cat("\n Between-individual component         = ", Ris[1,3])
 cat("\n Total Niche Width for the population = ", Ris[1,4])
 cat("\n The value of WIC/TNW is: ", Ris[1,5])
 cat("\n The p-value is: ", pvalue, "\n")
 if(Ris[1, 1] > 0) {
 if (as.integer(Ris[1,1]) == 1)  cat("\n Warning: ", as.integer(Ris[1,1]), " individual out of your population of ", as.integer(dataset$num.individuals))
    else cat("\n Warning: ", as.integer(Ris[1,1]), " individuals out of your population of ", as.integer(dataset$num.individuals))
 cat("\n have Shannon-Weaver scores equal to zero.  This may exaggerate")
 cat("\n the apparent degree of individual specialization, because these")
 cat("\n scores will drag the mean SWi score down and thus reduce WPC.")
 cat("\n This can be particularly misleading when an individual specializes")
 cat("\n entirely on one resource, which happens to be the most commonly")
 cat("\n consumed resource. This individual will have a Shannon-Weaver score")
 cat("\n of 0, even though its diet proportions are not unlike the proportions")
 cat("\n of the population as a whole.\n" )
 cat("\n" )
 if (as.integer(Ris[1,1]) == 1)  cat("\n The name of the individual is: \n") else cat("\n The names of the individuals are: \n")
 cat(Zeros, "\n")
}
  }
return(Ris2)
}
