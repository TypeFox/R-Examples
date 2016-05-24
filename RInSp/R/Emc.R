Emc <-
function(dataset, popd.type = "sum", index = "saramaki", replicates=999){
  #
  # The procedure will allow you to calculate likelihood measures of
  # niche breadth and overlap described in Petraitis (1979)
  #
  # Author: Nicola ZACCARELLI
  # E-mail: nicola.zaccarelli@gmail.com
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
 # for now only available for integer data
  if (class(dataset) != "RInSp") stop("The input must be an object of class RInSp.")
  if (dataset$data.type != "integer") stop("Input data type must be integer.")
  if (popd.type %in% c("sum", "average") == FALSE) stop("The specified population diet type is wrong.")
  if (popd.type == "sum")
    { dietpop = apply(dataset$resources, 2, sum) / sum(dataset$resources)} else {
      dietpop = apply(dataset$proportion, 2, mean)}
  if (index %in% c("saramaki", "barrat") == FALSE) stop("The specified distance index is wrong.")
  if (index == "saramaki") d.index = 1 else d.index = 2
  if(!is.double(replicates)) replicates = floor(abs(as.double(replicates)))
  if (replicates == 0) stop("The specified replicates value is wrong")
  if (replicates < 10) 
      { replicates = 10
        cat("\n Warning! Minimum number of replicates is 10. \n") }
  cat("\n Warning! resampling can take a while. Please be patient. \n")
totdieti = apply(dataset$resources, 1, sum)
pb = txtProgressBar(min=0, max=10,  char= "+", style = 3)
for (count in 1:10) {
     if (count < 10) StepCalc = floor(replicates / 10) else StepCalc = replicates - floor(replicates / 10) *9
     tmp = .Call("Emc", dataset$proportions, as.vector(d.index), as.vector(dietpop), as.vector(totdieti), as.vector(StepCalc), PACKAGE="RInSp")
     if (count ==1) Ris = tmp else Ris = rbind(Ris, tmp[-1, ])
     setTxtProgressBar(pb, count)
}
close(pb)
attributes(Ris)$dimnames[[2]] = c("Wmean", "E", "CW", "CwS")
cum.distr = ecdf(Ris[, 2])
pvalue= 1 - cum.distr(Ris[1, 2])
meannullE= mean(Ris[-1, 2])
Eadj= (Ris[1, 2] - meannullE) /(1 - meannullE)
CwS = Ris[1, 4]
Ris2= list(E = Ris[1, 2], meannullE= meannullE, Eadj= Eadj , p.value= pvalue, montecarlo= Ris, parameter = 2, pop.diet = popd.type, type.index = index)
class(Ris2) = "RInSp"
cat("\n Araujo's E index (Eobs): ", Ris[1,2])
cat("\n The mean Null E value is: ", meannullE)
cat("\n The E adjusted value is: ", Eadj)
cat("\n The p-value for P(Esim => Eobs) is: ", pvalue)
cat("\n Degree of clustering in the network (Cws) ", CwS)
# let's calculate the Monte Carlo resampling probability
cum.distr = ecdf(Ris[, 4])
pvalue= 1 - cum.distr(Ris[1, 4])
cat("\n The p-value for P(Cws.sim => Cws.obs) is: ", pvalue)
cat("\n Population diet type : ", popd.type)
cat("\n Cluster index : ", index)
cat("\n")
 return(Ris2)
}
