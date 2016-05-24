#compute t value using matrixeQTL package

tMatFunction <- function(mytrait,mygt,fileName=""){
  myfileName <- paste0("tMat",fileName)
	#these names is useful to track the order of the output
	if( is.null(rownames(mygt) ) ) {
        rownames(mygt) <- paste( "mk", 1:nrow(mygt), sep="")
	}
  if( is.null(rownames(mytrait) ) ) {
        rownames(mytrait) <- paste( "pheno", 1:nrow(mytrait), sep="")
  }

	#library(MatrixEQTL)
  snps1       <- SlicedData$new( mygt)
  gene1       <- SlicedData$new( mytrait )
  #modelLINEAR <- 117348  #defined in MatrixEQTL package
  try(sink("message_matrixEQTL.txt"), silent=TRUE)
  me          <- Matrix_eQTL_main(snps = snps1,gene = gene1, output_file_name = myfileName , pvOutputThreshold = 1, useModel = modelLINEAR,
                                  errorCovariance = numeric(),verbose = FALSE, pvalue.hist = FALSE )
  try(sink(), silent=TRUE)
  unlink( myfileName )# remove the output file
  try(unlink("message_matrixEQTL.txt"), silent=TRUE)
  
  res.mat     <- me$all$eqtls
  
  my.t.mat    <- NULL
  for (i.mk in 1:nrow(mygt)){
    ind.i.mk  <- which(as.character(res.mat[,1])==(rownames(mygt)[i.mk])) #first colum: mk
    res.i.mk  <- res.mat[ind.i.mk,]
    t.i.mk    <- res.i.mk[ match(rownames(mytrait), as.character(res.i.mk[,2])),3] #2nd column: trait, 3rd colum: statistic with sign(t-statistic or sign*sqrt(f-statistic))
    my.t.mat  <- cbind(my.t.mat, t.i.mk)
  }
  my.t.mat           <- abs(my.t.mat)
  rownames(my.t.mat) <- rownames(mytrait)
  colnames(my.t.mat) <- rownames(mygt)
 
  return(my.t.mat)
}