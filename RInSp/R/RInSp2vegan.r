RInSp2vegan <-
function(dataset){
  #
  # The procedure converts an RIS class object into a dataframe to be used as input for the "Vegan" package its nestedness indices.
  #
  # Author: Nicola ZACCARELLI, Giorgio MANCINELLI, Dan BOLNICK
  # E-mail: nicola.zaccarelli@gmail.com,
  #         giorgio.mancinelli@unisalento.it
  #         danbolnick@mail.texas.edu
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
# some checking 
if (class(dataset) != "RInSp") stop("The input must be an object of class RInSp")
if (dataset$data.type != "proportions") 
{Ris = dataset$proportions
 Ris[Ris > 0] = 1
} else {
Ris = dataset$resources
 Ris[Ris > 0] = 1
}
Ris = as.data.frame(Ris)
return(Ris)
}
