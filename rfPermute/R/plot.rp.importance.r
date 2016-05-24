#' @title Plot Random Forest Importance Distributions.
#' @description Plot the Random Forest importance distributions,
#' with significant p-values as estimated in rfPermute.
#' 
#' @param x An object produced by a call to \code{\link{rp.importance}}.
#' @param alpha Critical alpha to identify "significant" predictors.
#' @param sig.only Plot only the significant (<= \code{alpha}) predictors?
#' @param type character vector listing which importance measures to plot.
#'   Can be class names or names of overall importance measures 
#'   (e.g., "MeanDecreaseAccuracy") in the \code{\link{rp.importance}} object.
#' @param n Plot \code{n} most important predictors.
#' @param main Main title for plot.
#' @param ... Optional arguments which will be ignored.
#' 
#' @details The function will generate a panel of plots, one for each 
#'   importance type.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom ggplot2 ggplot aes geom_bar labs coord_flip theme scale_fill_manual
#' @importFrom gridExtra grid.arrange
#' @export
#' 
plot.rp.importance <- function(x, alpha = 0.05, sig.only = FALSE, 
                               type = NULL, n = NULL, main = NULL, ...) { 
  cols <- if(is.null(type)) {
    lapply(seq(1, ncol(x), 2), function(i) c(i, i + 1))
  } else {
    type <- unique(gsub(".pval", "", type))
    not.found <- setdiff(type, colnames(x))
    if(length(not.found) > 0) {
      not.found <- paste(not.found, collapse = ", ")
      stop(paste("the following columns in 'type' can't be found in 'x':", not.found))
    }
    lapply(match(type, colnames(x)), function(i) c(i, i + 1))
  }
  
  imp.list <- lapply(cols, function(i) {
    imp.df <- as.data.frame(x[, i])
    colnames(imp.df) <- c("imp", "pval")
    imp.df$pred <- rownames(imp.df)
    imp.df$is.sig <- factor(imp.df$pval <= alpha)
    imp.df <- imp.df[order(imp.df$imp), ]
    rownames(imp.df) <- 1:nrow(imp.df)
    imp.df <- imp.df[order(imp.df$imp, decreasing = TRUE), ]
    if(sig.only) imp.df <- imp.df[as.logical(imp.df$is.sig), ]
    if(!is.null(n) & is.numeric(n)) imp.df <- imp.df[1:min(n, nrow(imp.df)), ]
    with(imp.df, 
         ggplot(imp.df, aes(reorder(pred, imp), imp)) + 
           geom_bar(aes(fill = is.sig), stat = "identity") +
           labs(title = colnames(x)[i], x = "", y = "Importance") + 
           coord_flip() + theme(legend.position = "none") +
           scale_fill_manual(values = c("FALSE" = "black", "TRUE" = "red"))
    )
  })
  suppressWarnings(do.call(gridExtra::grid.arrange, c(imp.list, top = main)))
}