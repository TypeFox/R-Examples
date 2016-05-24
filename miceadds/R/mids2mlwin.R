mids2mlwin <- function(imp, file.prefix, path = getwd(), sep = " ", dec = ".", silent = FALSE,
  X = NULL )
{
  m <- imp$m
  i <- imp$iteration
  n <- nrow(imp$data)
  varnames <- names(imp$imp)[!sapply(imp$imp, is.null)]
  impnames <- paste0(file.prefix, seq(1, m, 1), ".txt")
  f1 <- paste0( file.prefix , "__impvals.txt" ) 
  base::write(x = varnames, file = f1 , append = F, ncolumns = length(varnames), sep = " ")
  base::write(x = c(n, m),  file = f1 , append = T, ncolumns = 1)
  base::write(x = impnames, file = f1 , append = T, ncolumns = m, sep = "\t")
  for (k in 1:m){
    h1 <- mice::complete(imp, k)
	if ( ! is.null(X) ){
		h1 <- data.frame(  X , h1	)										
	  				   }	
    utils::write.table(x = h1 , file = impnames[k], dec = dec, sep = sep, 
                col.names = F, row.names = F, quote = F, eol = "\n", )
  }
  if (!silent){
    cat("Data Values written to", paste0(path, "/", file.prefix, 1, ".txt"), 
        "through", paste0(file.prefix, m, ".txt"), "\n")
    cat("Imputation master file written to", paste0(path, "/", f1 ), "\n")
  }
}
