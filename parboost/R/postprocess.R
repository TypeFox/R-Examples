##' Postprocess the submodels
##'
##' Either euqally weigh all submodels or perform (penalized) regression
##' on them to determine the weights.
##' @title Postprocess parboost ensemble components
##' @param response The response vector.
##' @param family Boosting family.
##' @param list_of_models List of submodels.
##' @param postprocessing Type of postprocessing.
##' @return Predictions of postprocessed ensemble and a list of mboost models
##' @author Ronert Obst
##' @keywords internal
postprocess <- function(response, family, list_of_models, postprocessing, cores_cv = detectCores()) {
    models <- lapply(list_of_models, function(x) x[[1]])
    predictions <- lapply(list_of_models, function(x) x[[2]])

    predictions <- t(laply(predictions, function(x) x))

    if (Sys.info()[1] != "Windows") registerDoParallel(cores = cores_cv)

    none <- function(predictions) {
        ### simple mean of the subsample models
        predict_parboost <- function(models, newdata, type) {
            model_predictions <- laply(models, predict, newdata, type = type)
            return(colMeans(model_predictions))
        }

        fitted <- rowMeans(predictions)
        return(list(predict = predict_parboost, fitted = fitted, ensemble_weights = rep(1/length(models), length(models))))
    }

    glmfit <- function(response, predictions) {
        ### fit glm to the submodel predictions
        if (family == "gaussian") glm_family <- gaussian()
        if (family == "binomial") glm_family <- binomial()
        if (family == "poisson") glm_family <- poisson()
        data <- data.frame(predictions, response)
        glm_model <- glm(formula = response ~ ., family = glm_family, data=data)
        predict_parboost <- function(models, newdata, type) {
            model_predictions <- as.data.frame(t(laply(models, predict, newdata)))
            names(model_predictions) <- paste("X", 1:length(names(model_predictions)), sep="")
            postprocessed_predictions <- predict(glm_model, model_predictions, type = type)
            return(postprocessed_predictions)
        }
        fitted <- predict(glm_model, type = "link")
        return(list(predict = predict_parboost, fitted = fitted, ensemble_weights = coef(glm_model)))
    }

    lasso <- function(response, predictions) {
        lasso_model <- cv.glmnet(x=predictions, y=response, family = family, alpha = 1, parallel = TRUE)
        predict_parboost <- function(models, newdata, type) {
            model_predictions <- t(laply(models, predict, newdata))
            postprocessed_predictions <- predict(lasso_model, model_predictions, type = type)
            return(postprocessed_predictions)
        }

        fitted <- predict(lasso_model, newx=predictions)
        return(list(predict = predict_parboost, fitted = fitted, ensemble_weights = as(coef(lasso_model), "matrix")))
    }

    ridge <- function(response, predictions) {
        ridge_model <- cv.glmnet(x=predictions, y=response, family = family, alpha = 0, parallel = TRUE)
        predict_parboost <- function(models, newdata, type) {
            model_predictions <- t(laply(models, predict, newdata))
            postprocessed_predictions <- predict(ridge_model, model_predictions, type = type)
            return(postprocessed_predictions)
        }

        fitted <- predict(ridge_model, newx=predictions)

        return(list(predict = predict_parboost, fitted = fitted, ensemble_weights = as(coef(ridge_model), "matrix")))
    }

    elasticnet <- function(response, predictions) {
        elasticnet_model <- cv.glmnet(x=predictions, y=response, family = family, alpha = 0.8, parallel = TRUE)
        predict_parboost <- function(models, newdata, type) {
            model_predictions <- t(laply(models, predict, newdata))
            postprocessed_predictions <- predict(elasticnet_model, model_predictions, type = type)
            return(postprocessed_predictions)
        }

        fitted <- predict(elasticnet_model, newx=predictions)

        return(list(predict = predict_parboost, fitted = fitted, ensemble_weights = as(coef(elasticnet_model), "matrix")))
    }

    postprocessed_model <- switch(postprocessing,
                                  "none" = none(predictions),
                                  "glm" = glmfit(response, predictions),
                                  "lasso" = lasso(response, predictions),
                                  "ridge" = ridge(response, predictions),
                                  "elasticnet" = elasticnet(response, predictions))

    return(list(postprocessed_model = postprocessed_model,
                models = models))
}
