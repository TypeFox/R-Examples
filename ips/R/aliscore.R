aliscore <- function(x, gaps = "5state", w = 6, r, t, l, s, o, 
                     path = "/Applications/Aliscore_v.2.0"){
	
	rwd <- getwd()
	setwd(path)
	write.fas(x, "input.fas")
  
  ## parse options
  ## -------------
  N <- ifelse( gaps == "5state", "", "-N") # treatment of gaps
  w <- paste("-w", w) # window size
  r <- ifelse( missing(r), "", paste("-r", r ))
  if ( !missing(t) ) stop("option -t not yet implemented")
	if ( !missing(l) ) stop("option -l not yet implemented")
	if ( !missing(s) ) stop("option -s not yet implemented")
	o <- ifelse( missing(o), "", paste("-o", paste(o, collapse = ",") ))
  
	call <- paste("perl Aliscore.02.2.pl -i input.fas",
                N, w, r, o)
      
	system(call)
	id <- scan("input.fas_List_l_all.txt",
             sep = " ", quiet = TRUE)
	id <- as.numeric(id)
	nid <- length(id)
	if ( nid == 0 ) {
	  cat("\nALISCORE did not remove any characters.")
	} else {
	    x <- x[, -id]
	    cat("\nALISCORE removed", nid, "characters.")
	}
	setwd(rwd)
	x
}