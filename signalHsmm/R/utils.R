#use cut-off to calculate new probability
rescale <- function(x)
  ifelse(x > 0.025, (x - 0.025)/0.975 * 0.5 + 0.5, x*20)

#later import function from biogram
degenerate <- function(seq, element_groups) {
  tmp_seq <- seq
  if (!all(unique(tmp_seq) %in% unlist(element_groups))) {
    warning("'seq' contains elements not present in 'element_groups'. Such elements will be replaced by NA.")
    tmp_seq[!(tmp_seq %in% unlist(element_groups))] <- NA
  }
  
  if(is.null(names(element_groups))) {
    warning("'element_groups' is unnamed. Assumed names of groups are their ordinal numbers.")
    names(element_groups) <- 1L:length(element_groups)
  }
  
  if(length(unique(names(element_groups))) != length(names(element_groups))) {
    stop("'element_groups' must have unique names.")
  }
  
  for (i in 1L:length(element_groups)) {
    tmp_seq[tmp_seq %in% element_groups[[i]]] <- names(element_groups)[i]
  }
  
  if(class(seq) == "matrix")
    dim(tmp_seq) <- dim(seq)
  
  tmp_seq
}
