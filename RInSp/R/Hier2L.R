Hier2L = function(dataset, factor=1, weight.type="N_items") {
  #
  # Variance partition and WIC/TNW
  #
  # Author: Nicola ZACCARELLI
  # E-mail: nicola.zaccarelli@gmail.com
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
  # some checking
 if (class(dataset) != "RInSp") stop("The input must be an object of class RIS")
 if (dataset$data.type != "double") stop("Input data type must be double.")
 if ((factor %in% c(1:dim(dataset$info)[2])) == FALSE) stop("Wrong factor column number.")
 f.values = as.factor(dataset$info[, factor])
 f.values = f.values[1:length(f.values), drop=T] #get rid unused levels
 f.levels = levels(f.values)
 rawdata = dataset$resources
 ris = matrix(0, 5, nlevels(f.values)+1)
 colnames(ris) = c("All", f.levels)
 rownames(ris) = c("WIC", "BIC", "TNC", "WIC/TNW", "BGC")
 # Calculate full set first
 tmpRIS = WTcMC(dataset, weight=weight.type, replicates=10, print.ris=FALSE)
 ris[1:4, 1]= t(tmpRIS$montecarlo[1, ])
 for (i in 1:nlevels(f.values)){
  tmp = subset(rawdata, f.values == f.levels[i])
  tmp = import.RInSp(tmp, print.messages=FALSE)
  tmpRIS = WTcMC(tmp, weight=weight.type, replicates=10, print.ris=FALSE)
  ris[1:4, i+1]= t(tmpRIS$montecarlo[1, ])
 }
 BGC = var(by(rowMeans(rawdata), f.values, mean)) * (nlevels(f.values) -1)/(nlevels(f.values))#var sample and not pop!
 ris[5, 1] = BGC
 cat("\n ***************************************")
 names = c("All", f.levels)
 for (i in 1:(nlevels(f.values)+1)) cat("\n WIC/TNW of Level: ", names[i], "= ",  ris[4, i])
 cat("\n BGC             :    ", ris[5, 1])
 cat("\n BGC/TNW (%)     :    ", ris[5, 1]*100/ris[3, 1] ,"\n")
 cat("*************************************** \n")
 return(ris)
}