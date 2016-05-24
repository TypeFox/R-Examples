read.phy <- function(x){
	
	if (length(x) == 1)
		x <- scan(x, what = character(), sep = "\n", quiet 	= TRUE)
	semikolon <- grep(";", x)
	if (length(semikolon) > 0){
		x <- gsub(";", "", x)
		warning("\';\' eliminated from input file.")
	}
  ## number of taxa and characters
  ## -----------------------------
	n <- unlist(strsplit(x[1], split = "[[:space:]]+"))
  n <- as.integer(n[n != ""])
	names(n) <- c("ntax", "nchar")
  x <- x[-1]
  
	## taxon names
  ## -----------
	taxnames <- x[1:n["ntax"]]
	getTaxName <- function(x) 
    unlist(strsplit(x, split = "[[:blank:]]"))[1]
  taxnames <- sapply(taxnames, getTaxName, USE.NAMES = FALSE)
  
  ## concatenate characters
  ## ----------------------
  concatenate <- function(id, x, ntax, nblock){
    x <- x[seq(id, nblock, by = ntax)]
    x <- paste(x, collapse = "")
    x <- unlist(strsplit(x, split = "[[:blank:]]"))
    x <- tolower(paste(x[-1], collapse = ""))
    strsplit(x, "")
  }
  obj <- sapply(1:n["ntax"], concatenate, x = x, 
                ntax = n["ntax"], nblock = length(x))
	names(obj) <- taxnames
	obj <- as.DNAbin(obj)
	as.matrix(obj)
}
