NODF = function(dataset, print.results= TRUE){
  #
  # NODF (acronym for "nestedness metric based on overlap and decreasing fill")
  # is calculated for a whole matrix following Almeida-Neto et al. (2008).
  #
  # Author: Nicola ZACCARELLI
  # E-mail: nicola.zaccarelli@gmail.com
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
if (class(dataset) != "RInSp") stop("The input must be an object of class RInSp.")
if (dataset$data.type %in% c("integer", "double") == FALSE) stop("Wrong data type. Must be 'integer' or 'double'")
rows = dataset$num.individuals
cols = dataset$num.prey
R = matrix(0, rows, cols)
R[dataset$resources > 0] = 1
MT = apply(R, 1, sum)
NpR = matrix(0, rows, rows)
for (i in 1:rows) {
   for (j in 1:rows) {
	  if (MT[i] > MT[j]) NpR[j, i] = sum(R[i, ] * R[j, ]) / MT[j]
   }
}
NpC = matrix(0, cols, cols)
MT = apply(R, 2, sum)
for (i in 1:cols) {
   for (j in 1:cols) {
	  if (MT[i] > MT[j]) NpC[j, i] = sum(R[ ,i] * R[ ,j]) / MT[j]
   }
}
NTotal = ( (sum(NpR) + sum(NpC)) * 200) /(rows*(rows-1) + cols*(cols-1))
Ris = list(NODF = NTotal, Nrows = sum(NpR)*200 / (rows *(rows -1)), Ncols = sum(NpC)*200 / (cols *(cols -1)), R = R, NpR = NpR, NpC = NpC)
if (print.results == TRUE) cat("\n The value of NODF is ", NTotal, "\n")
class(Ris) = "RInSp"
return(Ris)
}