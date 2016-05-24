network <- function(adjacency = NULL,  weight = "autoShreve", max.df = NULL){
#   by = NA,
#   by.var <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
#   if (by.var == ".") stop("by=. not allowed")
#   if(!is.na(by)) warning("'by' argument is not yet supported and will be ignored", immediate. = TRUE)
  ret <- list(adjacency = adjacency,  weight = weight, max.df = max.df, by = "NA")
  class(ret) <- "network.spec"
  ret
}
