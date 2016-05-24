print.epPCA.inference.battery <-
function (x,...) {
	
  res.epPCA.inference.battery <- x
  if (!inherits(res.epPCA.inference.battery, "epPCA.inference.battery")) stop ("no convenient data")
  cat("**Results for Principal Component Analysis Inference Battery**\n")
  cat ("Permutation was performed on ", ncol(res.epPCA.inference.battery$components$eigs.perm),
       "components and bootstrap performed on ", nrow(res.epPCA.inference.battery$fj.boots$tests$boot.ratios), " variables\n")
  cat("*The results are available in the following objects:\n\n")
  res <- array("", c(2, 2), list(1:2, c("name", "description")))
  
  res[1,] <- c("$components","p-values ($p.vals) and permutations ($eigs.perm) for each component.")
  res[2,] <- c("$fj.boots","A list with bootstrap data and tests.")
  
  print(res)

}
