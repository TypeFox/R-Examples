### prediction type must be "logistic" or "linear"

predict.sofia <- function(object, newdata
  , prediction_type ### i removed default arguments, we can add them back but something confused me 
  , ...
) {

  sofia_facade <- new(RSofiaFacade)

  if(object$par$learner_type %in% c("logreg-pegasos") && prediction_type == "linear")
    warning(paste("linear prediction type used with learner type:", object$par$learner_type ))

  if(object$par$learner_type %in% 
    c("pegasos", "sgd-svm", "passive-aggressive", "margin-perceptron", "romma") &&
      prediction_type == "logistic"
    )
    warning(paste("logistic prediction type used with learner type:", object$par$learner_type))

  # Check if data needs processing
  # implies data was fit with sofia.fit (or sofia.character)
  if (is.null(object$formula)) {
    p <- sofia_facade$predict(object$weights, newdata, object$par$no_bias_term, prediction_type)  
  } else {
    parsed <- parse_formula(object$formula, newdata)
    p <- sofia_facade$predict(object$weights, parsed$data, object$par$no_bias_term, prediction_type)
  }
    
  return(p)

}
