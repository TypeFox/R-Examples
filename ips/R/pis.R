pis  <- function(x, what = "fraction", use.ambiguities = FALSE){

  if ( !inherits(x, "DNAbin") ) 
    stop("'x' is not of class 'DNAbin'")
  
  what <- match.arg(what, c("absolute", "fraction", "index"))
	
  if ( use.ambiguities ){
    warning("'use.ambiguities' is currently ignored ", 
            "and IUPAC ambiguity symbols are treated as missing data")
    use.ambiguities <- FALSE
  }
    
  pars.inf <- function(x){
		x <- table(x)
		x <- x[x > 1] # drop apomorphic chars
		n <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w", "y")
		if (length(x[!names(x) %in% n]) > 1)				
      x  <-  TRUE									
    else 
      x  <-  FALSE
	}
	x  <-  as.character(x)
	out  <-  apply(x, 2, pars.inf)
  if ( what %in% c("absolute", "fraction") ){
    out <- length(out[out])
    if ( what == "fraction" ){
      out <- round(out / ncol(x) * 100, digits = 2)
    } 
  } else {
    out <- which(out)
  }
	out
}