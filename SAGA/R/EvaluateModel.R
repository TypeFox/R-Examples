EvaluateModel <-
function(data, model, cex.axis=1, cex.names=1, cex.main=1){
  result <- data[[1]][[model]]
  pars <- c("Mean", strsplit(names(data[[1]])[model], split = ", ")[[1]])
  par.est <- summary(result)$coefficients[, 1]
  se.est <- summary(result)$coefficients[, 2]
  max.val <- max(par.est + se.est)
  min.val <- min(par.est - se.est)
  if(max.val < 0){
    max.val <- 0
  }
  if(min.val > 0){
    min.val <- 0
  }
  mp <- barplot(par.est, names.arg = pars, 
                main = "Single Model Means and Cond. SE",
                ylim = c(min.val - 0.2, max.val + 0.2), cex.axis=cex.axis,
                cex.names=cex.names, cex.main=cex.main) 
  
  high.se <- par.est + se.est
  low.se <- par.est - se.est
  segments(mp, high.se, mp, low.se, lwd = 3)
  segments(mp - 0.1, high.se, mp + 0.1, high.se, lwd = 3)
  segments(mp - 0.1, low.se, mp + 0.1, low.se, lwd = 3)
  results <- matrix(, 2, length(pars))
  results[1, ] <- par.est
  results[2, ] <- se.est
  colnames(results) <- pars
  row.names(results) <- c("Estimate", "Cond. SE")
  model.p <- 1 - pchisq(result$deviance, result$df.residual, lower.tail = T)
  mult.test <- 0.05 / length(data)
  final.results <- list()
  final.results[[1]] <- results
  final.results[[2]] <- model.p
  final.results[[3]] <- mult.test
  names(final.results) <- c("Estiamtes", "P.val", "Bonf.Corr.Alpha")
  return(final.results)
}
