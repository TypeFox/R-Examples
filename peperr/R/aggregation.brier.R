aggregation.brier <- function(full.data=NULL, response, x, model, cplx=NULL,  type=c("apparent", "noinf"), 
   fullsample.attr=NULL, ...){
   data <- as.data.frame(x)
   data$response <- response
   if(class(model)[1]=="penfit"){
      probs <- predict(model, data=data, penalized=x, ...)
   } else {
if(class(model)[1] == "glm") {
            probs <- predict(model, newdata = data, penalized = x, type = "response", ...)
        } else {
      probs <- predict(model, data=data, type="response", ...)
   }
}
   type <- match.arg(type)
   if (type=="apparent"){
      brier.score <- sum((probs-response)^2)/length(response)
   }
   if (type=="noinf"){
      brier.score <- mean((matrix(response, length(response), length(response), byrow=TRUE) - probs)^2)
   }
   brier.score
}