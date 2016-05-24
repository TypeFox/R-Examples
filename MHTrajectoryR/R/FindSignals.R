FindSignals <- function(results){
  output <- data.frame(ATC=as.character(results$best$drugs.names[ which(results$best$best.model == 1) ]),
                       beta=results$best$best.fit.coef[-1], 
                       stringsAsFactors = F, 
                       row.names = NULL)
  output <- output[order(output$beta, decreasing = T),]
  output <- output[which(output$beta>0),]
  return(output)
}

