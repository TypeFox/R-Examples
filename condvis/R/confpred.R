confpred <- 
function (model, newdata)
{
    if (identical(class(model), "lm")){
        pred <- predict(model, newdata, interval = "confidence", type = 
            "response")
        upr <- pred[, "upr"]
        lwr <- pred[, "lwr"]
        return(cbind(lwr, upr))
    }
    if (identical(class(model), c("glm", "lm"))){
        pred <- predict(model, newdata, type = "link", se.fit = TRUE)
        upr <- pred$fit + (2 * pred$se.fit)
        lwr <- pred$fit - (2 * pred$se.fit)
        upr <- model$family$linkinv(upr)
        lwr <- model$family$linkinv(lwr)
        return(cbind(lwr, upr))
    }
    if (identical(class(model), c("gam", "glm", "lm")) && "mgcv.conv" %in% 
        names(model)){
        pred <- predict(model, newdata, type = "link", se.fit = TRUE)
        upr <- pred$fit + (2 * pred$se.fit)
        lwr <- pred$fit - (2 * pred$se.fit)
        upr <- model$family$linkinv(upr)
        lwr <- model$family$linkinv(lwr)
        return(cbind(lwr, upr))
    }
    NULL
} 