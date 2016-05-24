#' Export to GenePop/LinkDos format
#' 
#' @param x profiles object
#' @param file (optional) filename
#' @param digits number of digits used to encode alleles (2 or 3)
#' @param description first line of GenePop file
#' @details Exports a profiles object as a text file in the GenePop/LinkDos format. See http://genepop.curtin.edu.au/help_input.html for more information.
#' @examples
#' data(freqsNLsgmplus)
#' 
#' set.seed(123)
#' 
#' # sample a small reference db
#' x <- sample.profiles(N = 1e3,freqs=freqsNLsgmplus)
#' 
#' write.genepop(x)
write.genepop <- function(x,file="",digits=2L,description= "DNAprofiles export"){
  x[is.na(x)] <- 0L
  x.char <- formatC(unclass(x),width=digits,flag = "0")
  
  df <- data.frame(row.names = paste(seq(nrow(x)),",",sep=""))
  
  for(l.i in seq(ncol(x)/2)){
    df[[l.i]] <- paste(x.char[,l.i*2-1], x.char[,l.i*2], sep="")    
  }
  
  if (file == "")  f <- stdout()  else f <- file(file, "w")
  
  cat(description, file = f)
  cat("\n", file = f)
  cat(substr(colnames(x),start=1,nchar(colnames(x))-2)[2*seq(ncol(x)/2)], sep=", ",file=f)
  cat("\nPOP \n", file = f)
  
  write.table(df, file = f, sep = " ", quote = FALSE, col.names = FALSE)
  if (file != "") close(f)
}