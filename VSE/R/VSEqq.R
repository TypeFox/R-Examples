#' VSEqq
#'
#' This function will generate QQ plots of normal distribution
#' @param data A matrix outputted by the function VariantSetEnrichment
#' @param ... Arguments for qqnorm plot
#' @keywords VSE, boxplot
#' @examples
#' #Load pre-saved object "bca.vse" as an example VSE output
#' load(file.path(system.file("extdata", "vse_output.Rda", package="VSE")))
#' VSEqq(bca.vse)
#' @export
VSEqq <- function(data, ...){
  matrix_norm <- data[[2]]
  no_of_beds <- nrow(matrix_norm)
  beds <- row.names(matrix_norm)
  for (i in 1:no_of_beds){
    kst <- ks.test(matrix_norm[i,], "pnorm", mean=mean(matrix_norm[i,]), sd=sd(matrix_norm[i,]), exact = TRUE)
    tit <- paste0(beds[i], " (p-value = ", round(kst$p.value, 2), ")");
    qqnorm(matrix_norm[i,c(2:ncol(matrix_norm))], main = beds[i], sub=paste0(" (p-value = ", round(kst$p.value, 2), ")"), ...)
    qqline(matrix_norm[i,c(2:ncol(matrix_norm))], col=2)
  }
}
