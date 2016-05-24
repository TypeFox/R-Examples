pop.diet <-
function(dataset, prop = "sum"){
  #
  # The procedure calculates the population diet.
  #
  # Author: Nicola ZACCARELLI, Giorgio MANCINELLI, Dan BOLNICK
  # E-mail: nicola.zaccarelli@gmail.com,
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
if (class(dataset) != "RInSp") stop("The input must be an object of class RInSp.")
if (prop %in% c("sum", "average") == FALSE) stop("Wrong proportion option.")
if (prop == "sum") {
    diet = apply(dataset$resources, 2, sum) / sum(dataset$resources)} else {
    diet = apply(dataset$proportion, 2, mean)}
D = 1 / sum((diet)^2)
ris = list(popdiet = diet, popdtype = prop, richness= dataset$num.prey, D= D)
cat("\n The population diet is \n")
cat(diet)
cat("\n The population diet type is ", prop)
cat("\n The resource richness is ", dataset$num.prey)
cat("\n The resource D value is ", D)
cat("\n")
class(ris) = "RInSp"
return(ris)}

