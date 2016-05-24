#' Pretty printing of polyOptim (Bertini) output.
#'
#' Pretty printing of polyOptim (Bertini) output.
#' 
#' @param x an object of class polyOptim, bertini
#' @param lagrange show values of lagrange multipliers?
#' @param digits digits to round to
#' @param ... additional parameters
#' @usage \method{print}{polyOptim}(x, lagrange = FALSE, digits = 3, ...)
#' @return Invisible string of the printed object.
#' @export
#' @examples
#' 
#' # see ?polyOptim
#' 
#' 
print.polyOptim <- function(x, lagrange = FALSE, digits = 3, ...){
	
  ## argument checking and basic variable setting
  stopifnot(is.bertini(x))  
  
  ## get variables
  if(lagrange){
    vars <- c(x$variables$vars, x$variables$lams)
  } else {
    vars <- x$variables$vars
  }
  tuple <- paste0("(", paste(vars, collapse = ","), ")")
  p <- length(vars)
  if(is.null(x$variables$lams)){
    optimizationType <- "unconstrained"
  } else {
    optimizationType <- "constrained"
  }
  
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
  
  # cont up critical values
  nRCritVals <- nrow(x$real_optima)
  nGlobalMax <- sum(x$real_optima$optima == "global max")
  nGlobalMin <- sum(x$real_optima$optima == "global min")  
  
  ## message solutions found.
  if(nDistinctSolns == 1){
    cat(paste0("One real critical value ", tuple, " found.  "))
    #if(optimizationType == "constrained") cat(paste0(
    #  "(A global ", ifelse(nGlobalMax == 1, "maximum.)", "minimum.)")
    #))  	
  } else {
  	maxTerm <- ifelse(nGlobalMax == 1, "um", "a")
  	minTerm <- ifelse(nGlobalMin == 1, "um", "a")  	
    cat(paste0(nDistinctSolns, " critical values ", tuple, " found.  "))
    if(optimizationType == "constrained") cat(paste0(
      "(", 
        nGlobalMax, " global maxim", maxTerm, ", ", 
        nGlobalMin, " global minim", minTerm, ".",
      ")"
    ))  	
  }  
  cat("\n")
  



  ## print out solutions
  if(lagrange || is.null(x$variables$lams)){
  	printSolns <- x$real_optima
  } else {
    printSolns <- x$real_optima[,
      -which(names(x$real_optima) %in% x$variables$lams),
      drop = FALSE
    ]
  }
  
  printSolns[1:(p+1)] <- round(
    printSolns[1:(p+1)], digits = digits
  )
  
  if(all(printSolns[,1:p] == round(printSolns[,1:p]))){
    formattedSolns <- apply(
      format(printSolns[,1:p, drop = FALSE]), 1, function(v){    	
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
    formattedSolns <- apply(
      format(printSolns[,1:p, drop = FALSE]), 1, function(v){    	
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

  for(k in 1:length(formattedSolns)){
    cat(paste0("  ", formattedSolns[k], " -> ", format(printSolns$value)[k]))
    if(printSolns[k,p+2] == "global max" && nGlobalMax + nGlobalMin > 1)
      cat("  (global max)")
    if(printSolns[k,p+2] == "global min" && nGlobalMax + nGlobalMin > 1)
      cat("  (global min)")      
    cat("\n")
  }
  
}



