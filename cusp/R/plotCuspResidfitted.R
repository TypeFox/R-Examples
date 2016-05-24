`plotCuspResidfitted` <-
function (object, caption = "Residual vs Fitted", xlab = paste("Fitted (", 
    colnames(fitted(object))[1], " convention)", sep = ""), ylab = "Residual", 
    ...) 
{
    plot(fitted(object), resid(object), xlab = xlab, ylab = ylab, 
        ..., main = caption)
}

