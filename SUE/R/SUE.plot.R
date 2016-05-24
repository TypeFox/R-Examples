SUE.plot <-
function(fit){
plot(SUE.fitted.values(fit),SUE.residuals(fit),xlab="fitted.values",ylab="residuals")
}
