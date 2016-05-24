robustness <-
function(pampe.object){

  data <- pampe.object$data
  time.pretr <- 1:length(pampe.object$model$fitted.values)
  time.tr <- (length(pampe.object$model$fitted.values)+1):(length(pampe.object$counterfactual[,2]))
  treated <- as.character(pampe.object$model$terms[[2]])
  controls <- as.character(pampe.object$controls)
  controls <- which(colnames(data) %in% controls)
 

  
##First create a matrix to save results
robust.pred <- matrix(NA,
                      ncol=length(pampe.object$controls),
                      nrow=length(c(time.pretr,time.tr)))
rownames(robust.pred) <- rownames(pampe.object$counterfactual)

##Iteratively remove each one of the ctrls in model
for (i in 1:length(pampe.object$controls)){
  #Each time remove one of the controls
  fmla <- paste(paste("`", treated, "`", sep=""), " ~ ", paste(paste("`", pampe.object$controls[-i], "`", sep=""), collapse= "+"))
  robust.temp <- lm(fmla, data=data[time.pretr,])
  #Save the results
  robust.pred[,i] <- predict(robust.temp, data)[c(time.pretr,time.tr)]
}
##Name rows and matrix of the results
robust.pred <- cbind(pampe.object$counterfactual[,1],pampe.object$counterfactual[,2], robust.pred)
colnames(robust.pred) <- c("Actual", "Predict w/ all", paste("w/o", pampe.object$controls))



class(robust.pred) <- "robustness"
return(robust.pred)

}
