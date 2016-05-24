gamEst <- function(X1, x1, y1, v1){
    tmp <- lsfit(X1, y1, wt = 1/v1, intercept = FALSE)
    beta1 <- lsfit(log(x1), log(tmp$residuals^2))$coef
    g1 <- beta1["X"]
    g1
}
