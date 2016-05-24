Plot.coef.glmnet <- function(cvfit, betas){

op <- par(no.readonly = TRUE)
par(mar=c(4,4,2.5,1), mgp=c(2.5,1,0), mfrow=c(2,2))
fit <- cvfit$glmnet.fit
bet <- fit$beta[match(betas, rownames(fit$beta)),]

plot(fit, xvar="lambda", col="gray")
glmnet::plotCoef(bet, lambda = fit$lambda, df = fit$df, dev = fit$dev.ratio, xvar = "lambda", add=TRUE, col="red")
abline(v=log(cvfit$lambda.min), lty=3)
abline(v=log(cvfit$lambda.1se), lty=3)

glmnet::plotCoef(bet, lambda = fit$lambda, df = fit$df, dev = fit$dev.ratio, xvar = "lambda", add=FALSE, col="red")
abline(v=log(cvfit$lambda.min), lty=3)
abline(v=log(cvfit$lambda.1se), lty=3)

norm <- apply(abs(fit$beta), 2, sum)
plot(fit, xvar="norm", col="gray")
glmnet::plotCoef(bet, xvar = "norm", add=TRUE, col="red",
                  norm = norm, lambda = fit$lambda, df = fit$df, dev = fit$dev.ratio)
abline(v=norm[match(cvfit$lambda.min, cvfit$lambda)], lty=3)
abline(v=norm[match(cvfit$lambda.1se, cvfit$lambda)], lty=3)

plot(fit, xvar="dev", col="gray")
glmnet::plotCoef(bet, lambda = fit$lambda, df = fit$df, dev = fit$dev.ratio, xvar = "dev", add=TRUE, col="red")
abline(v=fit$dev.ratio[match(cvfit$lambda.min, cvfit$lambda)], lty=3)
abline(v=fit$dev.ratio[match(cvfit$lambda.1se, cvfit$lambda)], lty=3)

par(op)
}