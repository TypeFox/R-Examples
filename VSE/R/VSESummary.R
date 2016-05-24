#' VSESummary
#'
#' This function will compute the enrichment from a VSE matrix
#' @param data A matrix outputted by the function VariantSetEnrichment
#' @keywords VSE, enrichment
#' @examples
#' #Load pre-saved object "bca.vse" as an example VSE output
#' load(file.path(system.file("extdata", "vse_output.Rda", package="VSE")))
#' VSESummary(bca.vse)
#' @export
VSESummary <- function(data){
  p_values <- c()
  intersect_matrix <- data[[1]]
  matrix_scaled <- data[[3]]
  for (i in 1:nrow(matrix_scaled)){
    avs <- matrix_scaled[i,1]
    null <- matrix_scaled[i,c(2:ncol(matrix_scaled))]
    p <- 1 - pnorm(matrix_scaled[i,], mean(matrix_scaled[i,]), sd(matrix_scaled[i,]))
    p_values[i] <- p[1]
  }
  padjust_values <- p.adjust(p_values, method = "bonferroni")
  summary_df <- data.frame(region=row.names(matrix_scaled), avs=intersect_matrix[,1], enrichment=matrix_scaled[,1], p.value=p_values, p.adjusted=padjust_values)
  return(summary_df)
}
