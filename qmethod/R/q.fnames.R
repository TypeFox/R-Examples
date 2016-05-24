q.fnames <- function(results, fnames) {
  # Error checks 
  if (class(results) != "QmethodRes") stop("The object provided is not of class 'QmethodRes'")
  comb <- array(sapply(fnames, function(x) substring(x,1,1)))
  nos <- 0:9
  if(sum(comb %in% nos) > 0) stop("The names should not begin with a number")
  if (length(fnames) != results$brief$nfactors) stop(paste0("The names provided (", length(names), ") does not match the number of factors in the results (", results$brief$nfactors, ")"))
  if (max(nchar(fnames)) > 50) stop("The names provided are longer than 50 characters.")
  
  # Change factor names for meaningful names
  q.objects <- c("loa", "flagged", "zsc", "zsc_n")
  for (i in q.objects) colnames(results[[i]]) <- fnames
  # Factor characteristics
  rownames(results[[7]]$characteristics) <- fnames
  dimnames(results[[7]]$cor_zsc) <- list(fnames, fnames)
  dimnames(results[[7]]$sd_dif)  <- list(fnames, fnames)
  return(results)
}