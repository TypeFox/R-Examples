.refit <-
function (object, fitting.method = "quadratic",
jackknife.estimation = "quadratic", asymptotic = TRUE,
allowed.fitting = c("quad", "line", "nonl", "logl", "log2"),
allowed.jackknife = c("quad", "line", "nonl", "logl", FALSE), ...)
{
fitting.method <- substr(fitting.method, 1, 4)
if (object$fitting.method == fitting.method)
stop("Model is already fitted with the specified fitting method", call. = FALSE)
if (!any(fitting.method == allowed.fitting)) {
warning("Fitting method not implemented. Using: quadratic", call. = FALSE)
fitting.method <- "quad"
}
if (jackknife.estimation != FALSE)
jackknife.estimation <- substr(jackknife.estimation, 1, 4)
if (!any(jackknife.estimation == allowed.jackknife)) {
warning("Fitting method (jackknife) not implemented. Using: quadratic",
call. = FALSE)
jackknife.estimation <- "quad"
}
if (!any(names(object) == "variance.jackknife") && jackknife.estimation != FALSE) {
warning("Jackknife variance estimation is not possible, due to the lack of it in the supplied model. Will be ignored.", call. = FALSE)
jackknife.estimation <- FALSE
}
if (!any(names(object) == "variance.asymptotic") && asymptotic) {
warning("Asymptotic variance estimation is not possible, due to the lack of it in the supplied model. Will be ignored.", call. = FALSE)
asymptotic <- FALSE
}
cl <- class(object)
if (any(names(object) == "variance.asymptotic") && asymptotic == FALSE) {
# removing unwanted parts of the object
object <- object[setdiff(names(object), c("PSI", "c11", "a11", "sigma",
"sigma.gamma", "g", "s", "variance.asymptotic"))]
}
if (any(names(object) == "variance.jackknife") && jackknife.estimation == FALSE) {
# removing unwanted parts of the object
       object <- object[setdiff(names(object), c("extrapolation.variance",
"variance.jackknife", "variance.jackknife.lambda"))]
}
class(object) <- cl
estimates <- object$SIMEX.estimates[-1, -1]
lambda <- object$lambda
ncoef <- length(coef(object))
ndes <- dim(object$model$model)[1]
p.names <- names(coef(object))
SIMEX.estimate <- vector(mode = "numeric", length = ncoef)
switch(fitting.method,
"quad" = extrapolation <- lm(estimates ~ lambda + I(lambda^2)),
"line" = extrapolation <- lm(estimates ~ lambda),
"logl" = extrapolation <- lm(I(log(t(t(estimates) +
(abs(apply(estimates, 2, min)) + 1) *
(apply(estimates, 2, min) <= 0)))) ~ lambda),
"log2" = extrapolation <- fit.logl(lambda, p.names, estimates),
"nonl" = extrapolation <- fit.nls(lambda, p.names, estimates)
)
# security if nls does not converge
if (any(class(extrapolation) == "lm") && fitting.method == "log2")
fitting.method <- "logl"
# predicting the SIMEX estimate
switch(fitting.method,
"quad" = SIMEX.estimate <- predict(extrapolation,
newdata = data.frame(lambda = -1)),
"line" = SIMEX.estimate <- predict(extrapolation,
newdata = data.frame(lambda = -1)),
"nonl" = for (i in 1:length(p.names))
SIMEX.estimate[i] <- predict(extrapolation[[p.names[i]]],
newdata = data.frame(lambda = -1)),
"log2" = for (i in 1:length(p.names))
SIMEX.estimate[i] <- predict(extrapolation[[p.names[i]]],
newdata = data.frame(lambda = -1)) -
((abs(apply(estimates, 2, min)) + 1) *
 (apply(estimates, 2, min) <= 0))[i],
"logl" = SIMEX.estimate <- exp(predict(extrapolation,
newdata = data.frame(lambda = -1))) -
(abs(apply(estimates, 2, min)) + 1) *
(apply(estimates, 2, min) <= 0)
)
# jackknife estimation
if (jackknife.estimation != FALSE) {
variance.jackknife <- object$variance.jackknife.lambda[-1, -1]
switch(jackknife.estimation,
"quad" = extrapolation.variance <-
lm(variance.jackknife ~ lambda + I(lambda^2)),
"line" = extrapolation.variance <- lm(variance.jackknife ~ lambda),
"logl" = extrapolation.variance <- lm(I(log(t(t(variance.jackknife) +
(abs(apply(variance.jackknife, 2, min)) + 1) *
(apply(variance.jackknife, 2, min) <= 0)))) ~ lambda),
"nonl" = extrapolation.variance <-
fit.nls(lambda, 1:NCOL(variance.jackknife), variance.jackknife)
)
# variance.jackknife <- rbind(predict(extrapolation.variance, newdata = data.frame(lambda = -1)), variance.jackknife)
variance.jackknife2 <- vector("numeric",ncoef^2)
switch(jackknife.estimation,
"nonl" = for (i in 1:NCOL(variance.jackknife))
variance.jackknife2[i] <- predict(extrapolation.variance[[i]],
newdata = data.frame(lambda = -1)),
"quad" =  variance.jackknife2 <- predict(extrapolation.variance,
newdata = data.frame(lambda = -1)),
"line" = variance.jackknife2 <- predict(extrapolation.variance,
newdata = data.frame(lambda = -1)),
"logl" = variance.jackknife2 <- exp(predict(extrapolation.variance,
newdata = data.frame(lambda = -1))) -
(abs(apply(variance.jackknife, 2, min)) + 1) *
(apply(variance.jackknife, 2, min) <= 0)
)
variance.jackknife <- rbind(variance.jackknife2, variance.jackknife)
variance.jackknife.lambda <- cbind(c(-1, lambda), variance.jackknife)
variance.jackknife <- matrix(variance.jackknife[1, ], nrow = ncoef,
ncol = ncoef, byrow = TRUE)
dimnames(variance.jackknife) <- list(p.names, p.names)
object$variance.jackknife.lambda <- variance.jackknife.lambda
object$variance.jackknife <- variance.jackknife
object$extrapolation.variance <- extrapolation.variance
}
if (asymptotic) {
sigma <- object$sigma
s <- construct.s(ncoef, lambda, fitting.method, extrapolation)
d.inv <- solve(s %*% t(s))
sigma.gamma <- d.inv %*% s %*% sigma %*% t(s) %*% d.inv
g <- list()
switch(fitting.method,
"quad" = g <- c(1, -1, 1),
"line" = g <- c(1, -1),
"logl" = for (i in 1:ncoef) g[[i]] <-
c(exp(coef(extrapolation)[1, i] - coef(extrapolation)[2, i]),
-exp(coef(extrapolation)[1, i] - coef(extrapolation)[2, i])),
"log2" = for (i in 1:ncoef) g[[i]] <-
c(exp(coef(extrapolation[[i]])[1] - coef(extrapolation[[i]])[2]),
-exp(coef(extrapolation[[i]])[1] - coef(extrapolation[[i]])[2])),
"nonl" = for (i in 1:ncoef) g[[i]] <-
c(-1, -(coef(extrapolation[[i]])[3] - 1)^-1,
coef(extrapolation[[i]])[2] / (coef(extrapolation[[i]])[3] - 1)^2)
)
g <- diag.block(g, ncoef)
variance.asymptotic <- (t(g) %*% sigma.gamma %*% g) / ndes
dimnames(variance.asymptotic) <- list(p.names, p.names)
object$sigma.gamma <- sigma.gamma
object$g <- g
object$s <- s
object$variance.asymptotic <- variance.asymptotic
}
object$call$fitting.method <- fitting.method
object$call$jackknife.estimation <- jackknife.estimation
object$call$asymptotic <- asymptotic
object$SIMEX.estimates[1, ] <- c(-1, SIMEX.estimate)
object$coefficients <- as.vector(SIMEX.estimate)
names(object$coefficients) <- p.names
fitted.values <- predict(object,
newdata = object$model$model[, -1, drop = FALSE], type = "response")
object$fitted.values <- fitted.values
if (is.factor(object$model$model[, 1]))
object$residuals <-
as.numeric(levels(object$model$model[, 1]))[object$model$model[, 1]] -
fitted.values else object$model$model[, 1] - fitted.values
object$extrapolation <- extrapolation
return(object)
}

