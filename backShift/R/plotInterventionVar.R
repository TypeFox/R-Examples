#' Plots the estimated intervention variances.
#'
#' @param estIntVars A  (Gxp)-dimensional matrix with the estimated intervention 
#' variances returned by \code{backShift} (as \code{varianceEnv}). 
#' \code{G} is the number of unique environments, \code{p} is the number of variables.
#' @param trueIntVars A (Gxp)-dimensional matrix with the true intervention variances 
#'  if these are known (for simulations). By default this paramter is set to NULL.
#'
#'

plotInterventionVars <- function(estIntVars, trueIntVars = NULL){
  
  gg.df <- melt(as.data.frame(estIntVars))
  gg.df <- cbind(gg.df, rep(1:nrow(estIntVars), times = ncol(estIntVars)))
  
   if(!is.null(trueIntVars)){
    true.vars <- melt(trueIntVars)[,3]
    gg.df <- cbind(gg.df, true.vars)
    gg.df <- as.data.frame(gg.df)
    colnames(gg.df) <- c("variable", "est.variance", "env", "true.variance")
    gg.df <- melt(gg.df, c("variable","env"), variable.name = "Variance")
    
    p <- ggplot(gg.df, aes_string(x = 'env', y = 'value', aes_string(color = 'Variance')))
    p <- p + theme_classic(base_size = 8)
    p <- p + geom_line(aes_string(type = 'Variance', color = 'Variance'), size = 0.5)
    p <- p + scale_color_discrete(labels = c("Estimated int. variance", "True int. variance"))
  }else{
    gg.df <- as.data.frame(gg.df)
    colnames(gg.df) <- c("variable", "est.variance", "env")
    
    p <- ggplot(gg.df, aes_string(x = 'env', y = 'est.variance', aes_string(color = 'variable')))
    p <- p + theme_classic(base_size = 8)
    p <- p + geom_line(aes_string(group = 'variable', color = 'variable'), size = 0.5)
    p <- p + scale_color_discrete(name = "Variable")
  }
  
  p <- p + facet_wrap(~variable, scales = "free") 
  p <- p + scale_x_discrete(labels = 1:nrow(estIntVars))
  p <- p + theme(axis.text.x=element_text(angle=90, hjust=1))
  p <- p + xlab("Environment") + ylab("Intervention variance")
  p <- p + ggtitle("Estimated intervention variances (up to an offset)")
  p
}