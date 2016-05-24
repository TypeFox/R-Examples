"summary.treeshape" <-
function(object, ...){


tree=object
n=nrow(tree$merge)+1

cat("\n")
cat("Phylogenetic tree shape: ")
cat("object ",deparse(substitute(object))," of class 'treeshape'\n")	
cat("\n")
cat("Number of tips:", nrow(tree$merge)+1, "\n\n")
cat("Colless' shape statistic:", colless(tree), "\n")
cat("Expected value (Yule model): ",n*log(n)+(0.57721566-1-log(2))*n," ")
cat("Standard Deviation: ",sqrt(3-pi^2/6-log(2))*n,"\n")
cat("Expected value (PDA model): ", sqrt(pi)*n^(3/2)," ")
cat("Standard Deviation: ",sqrt(10/3-pi)*n^(3/2),"\n")
}

