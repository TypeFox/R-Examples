R2 <- function(model) {
    if (!any(class(model)=="lm")) stop("model must be class ``lm'' or ``glm''")
    yhat <- fitted(model)
    ehat <- residuals(model)
    y <- yhat + ehat
    cor(yhat,y)^2
}  
