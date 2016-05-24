summary.pampe <-
function(object, ... ){
  
  pampe.object <- object
  
  time.tr <- (length(pampe.object$model$fitted.values)+1):(length(pampe.object$counterfactual[,2]))
  controls <- pampe.object$controls
  avg.tr.effect <- mean(pampe.object$counterfactual[time.tr,1]-pampe.object$counterfactual[time.tr,2])
  model <- list(coefficients=summary(pampe.object$model)$coefficients,
                sigma=round(summary(pampe.object$model)$sigma,4),
                df=pampe.object$model$df.residual,
                r.squared=round(summary(pampe.object$model)$r.squared,3),
                adj.r.squared=round(summary(pampe.object$model)$r.squared,3),
                fstatistic=summary(pampe.object$model)$fstatistic) 
   res.table <- cbind(pampe.object$counterfactual, pampe.object$counterfactual[,1] - pampe.object$counterfactual[,2])
   colnames(res.table) <- c("Actual", "Counterfactual", "Estim. Tr. Effect")

  
 result <- list(time.tr=time.tr,controls=controls,avg.tr.effect=avg.tr.effect, model=model, res.table=res.table)
 class(result) <- "summary.pampe"
 return(result)
  
}


