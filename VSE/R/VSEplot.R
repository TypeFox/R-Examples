#' VSEplot
#'
#' This function will generate a figure for VSE analysis
#' @param data A list of matrices outputted by the function VariantSetEnrichment
#' @param padj Bonferroni adjusted p-value cutoff. Default: 0.01
#' @param ... Arguments from boxplot
#' @keywords VSE, boxplot
#' @examples
#' #Load pre-saved object "bca.vse" as an example VSE output
#' load(file.path(system.file("extdata", "vse_output.Rda", package="VSE")))
#' VSEplot(bca.vse, las=2,pch=20, cex=1, cex.main=0.6, padj=0.05)
#' @export
VSEplot <- function(data, padj=0.01, ...){
  matrix_norm <- data[[2]]
  matrix_scaled <- data[[3]]
  p_values <- c()
  for (i in 1:nrow(matrix_scaled)){
    avs <- matrix_scaled[i,1]
    null <- matrix_scaled[i,c(2:ncol(matrix_scaled))]
    p <- 1 - pnorm(matrix_scaled[i,], mean(matrix_scaled[i,]), sd(matrix_scaled[i,]))
    p_values[i] <- p[1]
  }
  padjust_values <- p.adjust(p_values, method = "bonferroni")
  boxplot(t(matrix_scaled), ylim=c( min(matrix_scaled),max(matrix_scaled) ), ylab="Enrichment Score", outline=FALSE, notch=TRUE, ...)
  for (i in 1:length(padjust_values)){
    points(i,matrix_scaled[i,1], pch=19, col=ifelse(padjust_values[i] <= padj, "red","black"))
  }
}
