WTcMC = function(dataset, replicates=999, weight="equal", print.ris=TRUE){
  #
  # The program calculates the Total Niche Width (TNW), and breaks TNW down into its
  # Between Individual Component (BIC) and Within Individual Component (WIC).
  # Continuous data are used.
  #
  # Author: Nicola ZACCARELLI, Giorgio MANCINELLI, Dan BOLNICK
  # E-mail: nicola.zaccarelli@gmail.com,
  #         giorgio.mancinelli@unisalento.it
  #         danbolnick@mail.texas.edu
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
  if (class(dataset) != "RInSp") stop("The input must be an object of class RInSp")
  if (dataset$data.type %in% c("integer", "double") != TRUE) stop("Input data type must be double or integer.")
  if (weight %in% c("equal","N_items") == FALSE) stop("`weight` must be either `equal` or `N_items`")
  if (weight == "equal") weight.opt = 1 else weight.opt = 2
  replicates = as.integer(replicates)
  if (replicates <=1) stop("Wrong number of replicates.")
  if (print.ris == TRUE) cat("\n If your dataset is big, this can take time. Please be patient. \n")
  # coerce vec to be double
  if (!is.double(dataset$resources)) dataset$resources = matrix(as.double(dataset$resources), dataset$num.individuals, dataset$num.prey)
  if(!is.integer(replicates)) replicates = as.integer(replicates)
  Ris = .Call("WTcMC", dataset$resources, as.vector(replicates), as.vector(weight.opt), PACKAGE="RInSp")
  attributes(Ris)$dimnames[[2]] = c("WIC", "BIC", "TNW", "WonT")
  cum.distr = ecdf(Ris[, 4])
  pvalue= cum.distr(Ris[1, 4])
  Ris2= list(WonT= Ris[1, 4], p.value= cum.distr(Ris[1, 4]), montecarlo= Ris, weight = weight, parameter = 4)
  class(Ris2) = "RInSp"
  if (print.ris==TRUE){
  cat("\n Results are based on Roughgarden's 1979 equations.")
  cat("\n The weighting is '")
  cat(weight)
  if (weight == "equal") { cat("' so each individual contributes")
                           cat("\n equally regardless of the diet items number.")}
      else { cat("' so each individual contributes")
             cat("\n proportionally to the diet items number.")}
  cat("\n Estimated parameters: ")
  cat("\n Within-individual component (WIC)          = ", Ris[1,1])
  cat("\n Between-individual component (BIC)         = ", Ris[1,2])
  cat("\n Total Niche Width for the population (TNW) = ", Ris[1,3])
  cat("\n The value of WIC/TNW is ", Ris[1, 4])
  cat("\n The p-value is ", pvalue)
  cat("\n The weighting option ", weight)
  cat("\n")
}
return(Ris2)
}
