like.Wi = function(dataset){
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
if (class(dataset) != "RInSp") stop("The input must be an object of class RSI")
Mat = t((apply(dataset$resources, 2, sum) / sum(dataset$resources)) / t(dataset$proportions))
Mat[Mat == "Inf"] = 1
Mat = Mat^dataset$resources
Ris = matrix(apply(Mat, 1, prod), dataset$num.individuals, 1)
# Calculate the chi-square degrees of freedom
tmp = dataset$resources
tmp[tmp > 0] = 1
dgf = apply(tmp, 1, sum)
pvalue = pchisq(-2*log(Ris), dgf, lower.tail=FALSE) 
tmp = cbind(Ris, pvalue, Ris^(1/apply(dataset$resources, 1, sum)))
rownames(tmp) = rownames(dataset$resources)
colnames(tmp) = c("Likelihood", "p-value", "Wi")
Ris = list(MeanWi = mean(tmp[, 3]), ResCat= dataset$num.prey, ind.vals= tmp)
cat("\nThe mean population Wi value is ", mean(tmp[, 3]), "\n")
cat("The number of resources is ", dataset$num.prey, "\n")
cat("Critical value of the chi-square distribution for", dataset$num.prey -1, "df and alpha of 0.05 is", qchisq(0.05, (dataset$num.prey -1), lower.tail=FALSE), "\n")
class(Ris) = "RInSp"
return(Ris)}
