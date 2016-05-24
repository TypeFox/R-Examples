`print.pedMap` <- function(x, nrow=5, ncol=10, ...){
   
  Xmap <- x$map[1:nrow,]
  Xfam <- x$fam[1:nrow,]
  Xgeno <- x$geno[1:nrow,1:ncol]
  
  cat("$map \n")
  print(Xmap,...)
  cat("...\n \n")
  cat("$fam \n")
  print(Xfam,...)
  cat("...\n \n")
  cat("$geno \n")
  print(Xgeno,...)
  cat("...\n \n")
} 
