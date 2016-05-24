initial <- function(y, X, offset = NULL, type = "Poi", alpha = 1) {

    if (type == "Poi"){
        (GLM <- glm(y ~ -1 + X, offset = offset, family = poisson, x = TRUE))
        GLMInt <- glm(y ~ X, offset = offset, family = poisson, x = TRUE)
    }


    if (type == "NegBin"){
      offset.nb<-if(is.null(offset)) rep(0,length(y)) else offset
        (GLM <- glm.nb(y ~ -1 + X + offset(offset.nb),
                       init.theta = alpha, x = TRUE))
        GLMInt <- glm.nb(y ~ X + offset(offset.nb),
                         init.theta = alpha, x = TRUE)
    }

    if (type == "Bin"){
        (GLM <- glm(y ~ -1 + X, offset = offset,
                    family = binomial(link = logit),
                    na.action = na.omit, x = TRUE))
        GLMInt <- glm(y ~ X, offset = offset,
                      family = binomial(link = logit),
                      na.action = na.omit, x = TRUE)
    }

    list(beta = GLM$coefficients, y = y, X = X, alpha = GLM$theta, type = type,
         null.deviance = GLMInt$null.deviance, df.null = GLMInt$df.null)

}
