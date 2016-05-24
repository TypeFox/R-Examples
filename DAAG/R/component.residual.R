"component.residual" <-
function (lm.obj, which = 1, xlab = "Component", ylab = "C+R") 
{
res <- residuals(lm.obj)
data <- model.matrix(lm.obj)
if (var(data[,1]) == 0) {data <- data[, -1]
lm.obj$coef <- lm.obj$coef[-1]
}
bx <- lm.obj$coef[which]
plot(data[,which], bx*data[,which]+res, xlab = xlab, ylab = ylab)
panel.smooth(data[,which], bx*data[,which]+res) 
}
