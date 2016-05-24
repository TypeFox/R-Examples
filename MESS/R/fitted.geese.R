fitted.geese <- function(object, ...) {

    # Get the linear predictor
    glmcall <- object$call
    glmcall$id <- glmcall$jack <- glmcall$control <- glmcall$corstr <- glmcall$waves <- glmcall$zcor <- glmcall$std.err <- glmcall$scale.fix <- glmcall$scale.value<- glmcall$z <- NULL
    glmcall[[1]] <- as.name("model.frame")

    mf <- eval(glmcall, parent.frame())

    X <- model.matrix(formula(object), data=mf)
    N <- nrow(X)


    offset <- model.offset(mf)
    if (is.null(offset)) 
        offset <- rep(0, N)

  value <- object
    
    
    value$offset <- offset
    if (is.null(value$offset)) 
        value$linear.predictors <- X %*% object$beta
    else value$linear.predictors <- value$offset + X %*% 
        object$beta


    fitted.values <- family(value)$linkinv(value$linear.predictors)

    fitted.values

}
