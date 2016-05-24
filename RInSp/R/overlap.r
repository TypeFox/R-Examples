overlap <- function(dataset){
  #
  # This option calculates the pairwise diet overlap between all individuals in the sample
  #
  # Author: Nicola ZACCARELLI, Giorgio MANCINELLI, Dan BOLNICK
  # E-mail: nicola.zaccarelli@gmail.com,
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
if (class(dataset) != "RInSp") stop("The input must be an object of class RInSp")
NRows = dataset$num.individuals
overlapmatrix = matrix(0, NRows, NRows)
for (i in 1:NRows) {
  for (k in 1:NRows) {
    overlapmatrix[i, k] = sum(pmin(dataset$proportions[i, ], dataset$proportions[k, ]))
  }
}
rownames(overlapmatrix) = rownames(dataset$resources)
colnames(overlapmatrix) = rownames(dataset$resources)
# calculating the mean overlap by excluding the diagonal elements
meanoverlap = (sum(overlapmatrix) - sum(diag(overlapmatrix))) / (NRows^2 - NRows)
# calculating the mean overlap for individual i excluding the diagonal elements
tmp = overlapmatrix
diag(tmp) = NA
meanindividualoverlap = matrix(rowMeans(tmp, na.rm = T), NRows, 1)
rownames(meanindividualoverlap) = rownames(dataset$resources)
Ris = list(meanoverlap = meanoverlap, meanindividualoverlap= meanindividualoverlap, meandissimilarity = 1- meanoverlap, overlapmatrix = overlapmatrix, parameter=0)
cat("\n The mean pairwise overlap is ", meanoverlap)
cat("\n The mean pairwise dissimilarity is ", 1 - meanoverlap, "\n")
return(Ris)}

