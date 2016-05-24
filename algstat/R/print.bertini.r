#' Pretty Printing of Bertini Output
#'
#' Pretty printing of Bertini output.
#' 
#' @param x an object of class bertini
#' @param digits digits to round to
#' @param ... additional parameters
#' @usage \method{print}{bertini}(x, digits = 3, ...)
#' @return Invisible string of the printed object.
#' @export
#' @examples
#' \dontrun{
#' 
#' # see ?bertini
#' 
#' variety("x^2 + 1")
#' variety(c("x^2 + 1 + y","y"))
#'
#' }
#' 
print.bertini <- function(x, digits = 3, ...){
	
  ## argument checking and basic variable setting
  stopifnot(is.bertini(x))  
  
  ## get variables
  vars <- str_replace(x$main_data[2], "Variables:  ", "")
  vars <- strsplit(vars, " ")[[1]]
  tuple <- paste0("(", paste(vars, collapse = ","), ")")
  p <- length(vars)
  
  ## determine number of solutions and kinds
  nFSolns <- nrow(x$finite_solutions); if(is.null(nFSolns)) nFSolns <- 0
  nNsSolns <- nrow(x$nonsingular_solutions); if(is.null(nNsSolns)) nNsSolns <- 0
  nSSolns <- nrow(x$singular_solutions); if(is.null(nSSolns)) nSSolns <- 0
  nRSolns <- nrow(x$real_finite_solutions); if(is.null(nRSolns)) nRSolns <- 0
  
  ## round
  fSolns <- round(x$finite_solutions, digits = digits)
  nsSolns <- round(x$nonsingular_solutions, digits = digits)
  sSolns <- round(x$singular_solutions, digits = digits)
  rfSolns <- round(x$real_finite_solutions, digits = digits)        
  
  ## make solns data frame, add in type variable (real or complex)
  solns <- as.data.frame(fSolns)  
  complexSolnIndic <- apply(solns, 1, function(x) any(Im(x) != 0)) # T's and F's
  solns$type <- "real"
  solns$type[complexSolnIndic] <- "complex"  
  
  ## count up solutions (for multiplicities)  
  fSolnsString <- apply(fSolns, 1, paste, collapse = " ")
  fSolnsTab <- table(unname(fSolnsString))
  nDistinctSolns <- length(fSolnsTab)  
  
  
  ## add in regularity (singular or nonsingular)
  solns$regularity <- ""  
  
  # count up nonsingular solutions
  if(nNsSolns > 0){
    nsSolnsString <- apply(nsSolns, 1, paste, collapse = " ")
    nsSolnsTab <-table(unname(nsSolnsString))
    nDistinctNsSolns <- length(nsSolnsTab) 
  } else {
  	nsSolnsString <- character(0)
  	nDistinctNsSolns <- 0
  }  
  solns$regularity[fSolnsString %in% nsSolnsString] <- "nonsingular"
  
  # count up singular solutions
  if(nSSolns > 0){
    sSolnsString <- apply(sSolns, 1, paste, collapse = " ")
    sSolnsTab <-table(unname(sSolnsString)) 
    nDistinctSSolns <- length(sSolnsTab)  
  } else {
  	sSolnsString <- character(0)
  	nDistinctSSolns <- 0
  }
  solns$regularity[fSolnsString %in% sSolnsString] <- "singular"
  
  # count up real solutions
  if(nRSolns > 0){
    rfSolnsString <- apply(rfSolns, 1, paste, collapse = " ")
    rfSolnsTab <-table(unname(rfSolnsString))
    nDistinctRSolns <- length(rfSolnsTab)
  } else {
  	rfSolnsString <- character(0)
  	nDistinctRSolns <- 0
  }  
  
  
  ## message solutions found.
  if(nDistinctSolns == 1){
  	rc <- ifelse(nRSolns > 0, "real", "complex")
  	sing <- ifelse(nNsSolns > 0, "nonsinguar", "singular")  	
    cat(paste0("One ", rc, ", ", sing, " solution ", tuple, " found."))  	
  } else {
    cat(paste0(nDistinctSolns, " solutions ", tuple, " found.  (", 
      nDistinctRSolns, " real, ", nDistinctSolns - nDistinctRSolns, " complex; ", 
      nDistinctNsSolns, " nonsingular, ", nDistinctSSolns, " singular.)"))  	
  }  
  cat("\n")
  
  # add on multiplicities by joining a similar dataframe  
  solns$soln <- fSolnsString 
  mults <- as.data.frame(fSolnsTab)
  names(mults) <- c("soln", "mult")
  uniqueSolns <- unique(merge(solns, mults))
  uniqueSolns <- uniqueSolns[,-which(names(uniqueSolns) == "soln")]

  ## print out solutions
  printSolns <- uniqueSolns[,1:p, drop = FALSE]
  if(all(printSolns == round(printSolns))){
    formattedSolns <- apply(format(printSolns), 1, function(v){    	
      s <- paste0("(", paste(v, collapse = ","), ")")    
      s <- str_replace_all(s, "0\\+0i", "0")    
      s <- str_replace_all(s, "0\\-0i", "0")        
      s <- str_replace_all(s, "0\\+", "")
      s <- str_replace_all(s, "0\\-", "-")
      s <- str_replace_all(s, "\\+0i", "")
      s <- str_replace_all(s, "\\-0i", "")    
      s <- str_replace_all(s, "1i", "i") 
      s
    })    
  } else {
    formattedSolns <- apply(format(printSolns), 1, function(v){    	
      s <- paste0("(", paste(v, collapse = ","), ")")    
      s <- str_replace_all(s, "0.000\\+0.000i", "0")    
      s <- str_replace_all(s, "0.000\\-0.000i", "0")        
      s <- str_replace_all(s, "0.000\\+", "")
      s <- str_replace_all(s, "0.000\\-", "-")
      s <- str_replace_all(s, "\\+0.000i", "")
      s <- str_replace_all(s, "\\-0.000i", "")    
      s <- str_replace_all(s, "1.000i", "i") 
      s <- str_replace_all(s, "\\+0i", "")
      s <- str_replace_all(s, "\\-0i", "")      
      s
    })
  }


  
  for(k in 1:nrow(uniqueSolns)){
  	if(uniqueSolns$mult[k] == 1){
      regu <- ifelse(uniqueSolns[k,"regularity"] == "nonsingular", "(R)", "(S)")  		
      cat( paste("   ", formattedSolns[k], regu))
    } else {
      regu <- ifelse(uniqueSolns[k,"regularity"] == "nonsingular", "(R, ", "(S, ")
      cat( paste0("    ", formattedSolns[k], " ", regu, uniqueSolns$mult[k], ")"))
    }
    cat("\n")    
  }   
  
}



