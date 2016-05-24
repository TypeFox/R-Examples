dropModel <- function(model, drop) {
    model2 <- update(model,
                         dropFormula(model, drop))
    plot(resid(model), resid(model2),
         xlab = "Original Residuals",
         ylab = paste("Residuals after dropping", drop))
    abline(0,1)
    model2
}
