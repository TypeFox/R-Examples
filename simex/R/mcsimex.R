mcsimex <-
function(
model, # The naive model
SIMEXvariable, # string/vector of strings with names of the SIMEXvariable
mc.matrix, # the misclassification Matrix must be a list
lambda = c(0.5, 1, 1.5, 2), # the values for Lambda
B = 100, # number of simulations
fitting.method = "quadratic", # fitting method for the extrapolation step
jackknife.estimation = "quadratic", # logical: do jackknife variance estimation
asymptotic = TRUE) # fitting method
{
fitting.method <- substr(fitting.method, 1, 4)
if (!any(fitting.method == c("quad", "line", "nonl", "logl", "log2"))) {
warning("Fitting Method not implemented. Using: quadratic", call. = FALSE)
fitting.method <- "quad"
}
if (jackknife.estimation != FALSE)
jackknife.estimation <- substr(jackknife.estimation, 1, 4)
if (!any(jackknife.estimation == c("quad", "line", "nonl", "logl", FALSE))) {
warning("Fitting Method (jackknife) not implemented. Using: quadratic",
call. = FALSE)
jackknife.estimation <- "quad"
}
if (any(lambda <= 0)) {
warning("lambda should not contain 0 or negative values. 0 or negative values will be ignored",
call. = FALSE)
lambda <- lambda[lambda >= 0]
}
if (!any(names(model) == "x") && asymptotic)
stop("The option x must be enabled in the naive model for asymptotic variance estimation",
call. = FALSE)
if (is.matrix(mc.matrix)) {
mc.matrix <- list(mc.matrix)
names(mc.matrix) <- SIMEXvariable
}
if (is.list(mc.matrix)) {
if (!all(check.mc.matrix(mc.matrix)))
stop("mc.matrix may contain negative values for exponents smaller than 1")
if (length(mc.matrix) != length(SIMEXvariable))
stop("mc.matrix and SIMEXvariable do not match")
}
if (any(!sapply(as.data.frame(model$model), is.factor)[SIMEXvariable]))
stop("SIMEXvariable must be a factor")
cl <- match.call()
ncoef <- length(model$coefficients)
ndes <- length(model$y)
nlambda <- length(lambda)
p.names <- names(coef(model))
factors <- lapply(SIMEXvariable, levels)
estimates <- matrix(data = NA, length(lambda) + 1, length(model$coefficients))
theta <- matrix(data = NA, B, ncoef)
colnames(theta)<- p.names
theta.all <- vector(mode = "list", nlambda)
if (jackknife.estimation != FALSE) {
var.exp <- list()
var.exp[[1]] <- extract.covmat(model)
}
if (asymptotic) {
psi <- matrix(rep(0, ndes * ncoef), ncol = ncoef, nrow = ndes)
psi <- residuals(model, type = "response") * model$x
PSI <- psi
am <- list()
xi <- model$x
a <- list()
dh <- model$family$mu.eta(model$linear.predictors)
for (k in 1:ndes) a[[k]] <- dh[k] * xi[k, ] %*% t(xi[k, ])
a.mat <- matrix(unlist(a), nrow = length(a), byrow =TRUE)
ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]), byrow = FALSE)
am[[1]] <- -ab / ndes
}
estimates[1, ] <- model$coefficients
for (i in 1:length(lambda)) {
if (jackknife.estimation != FALSE)
variance.est <- matrix(0, ncol = ncoef, nrow = ncoef)
if (asymptotic) {
psi <- matrix(0, ncol = ncoef, nrow = ndes)
a <- list()
for (k in 1:ndes) a[[k]] <- matrix(0, nrow = ncoef, ncol = ncoef)
}
for (j in 1:B) {
SIMEXdata <- data.frame(model$model)
# doing the misclassification
SIMEXv<- data.frame(SIMEXdata[, SIMEXvariable])
colnames(SIMEXv) <- SIMEXvariable
if (is.character(mc.matrix)) {
SIMEXdata[, SIMEXvariable] <-
eval(call(mc.matrix, SIMEXdata, lambda[i]))
} else {
SIMEXdata[, SIMEXvariable] <- misclass(SIMEXv, mc.matrix, lambda[i])
}
# updating the model and calculating the estimates
model.SIMEX <- update(model, data = data.frame(SIMEXdata))
theta[j, ] <- model.SIMEX$coefficients
if (jackknife.estimation != FALSE) {
variance.est <- variance.est + extract.covmat(model.SIMEX)
}
if (asymptotic) {
xi <- model.SIMEX$x
psi <- psi + (residuals(model.SIMEX, type = "response") * xi)
dh <- model$family$mu.eta(model.SIMEX$linear.predictors)
for(k in 1:ndes) a[[k]] <- a[[k]] - dh[k] * xi[k, ] %*% t(xi[k, ])
}
}
# taking the mean of the estimate -> SIMEX estimate
estimates[i + 1, ] <- colMeans(theta)
theta.all[[i]] <- theta
if (jackknife.estimation != FALSE) {
variance.est <- variance.est / B
s2 <- cov(theta)
var.exp[[i + 1]] <- variance.est - s2
}
if (asymptotic) {
xiB <- psi / B
PSI <- cbind(PSI, xiB)
a.mat <- matrix(unlist(a), nrow = length(a), byrow =TRUE)
ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]), byrow = FALSE)
am[[i + 1]] <- ab / (B * ndes)
}
}
lambda <- c(0, lambda)
colnames(estimates) <- p.names
# fitting the extrapolation function
switch(fitting.method,
"quad" = extrapolation <- lm(estimates ~ lambda + I(lambda^2)),
"line" = extrapolation <- lm(estimates ~ lambda),
"logl" = extrapolation <- lm(I(log(t(t(estimates) +
(abs(apply(estimates, 2, min)) + 1) *
(apply(estimates, 2, min) <= 0)))) ~ lambda),
"log2" = extrapolation <- fit.logl(lambda, p.names, estimates),
"nonl" = extrapolation <- fit.nls(lambda, p.names, estimates)
)
if (any(class(extrapolation) == "lm") && fitting.method == "log2")
fitting.method <- "logl"
# predicting the SIMEX estimate
SIMEX.estimate <- vector(mode = "numeric", length = ncoef)
switch(fitting.method,
"quad" = SIMEX.estimate <- predict(extrapolation,
newdata = data.frame(lambda = -1)),
"line" = SIMEX.estimate <- predict(extrapolation,
newdata = data.frame(lambda = -1)),
"nonl" = for(i in 1:length(p.names))
SIMEX.estimate[i] <- predict(extrapolation[[p.names[i]]],
newdata = data.frame(lambda = -1)),
"logl" = SIMEX.estimate <- exp(predict(extrapolation,
newdata = data.frame(lambda = -1))) -
(abs(apply(estimates, 2, min)) + 1) * (apply(estimates, 2, min) <= 0),
"log2" = for(i in 1:length(p.names))
SIMEX.estimate[i] <- predict(extrapolation[[p.names[i]]],
newdata = data.frame(lambda = -1)) -
((abs(apply(estimates, 2, min)) + 1) * (apply(estimates, 2, min) <= 0))[i]
)
if (jackknife.estimation != FALSE) {
variance.jackknife <- matrix(unlist(var.exp), ncol = ncoef^2 , byrow = TRUE)
switch(jackknife.estimation,
"quad" = extrapolation.variance <-
lm(variance.jackknife ~ lambda + I(lambda^2)),
"line" = extrapolation.variance <-
lm(variance.jackknife ~ lambda),
"logl" = extrapolation.variance <-
lm(I(log(t(t(variance.jackknife) +
(abs(apply(variance.jackknife, 2, min)) + 1) *
(apply(variance.jackknife, 2, min) <= 0)))) ~ lambda),
"nonl" = extrapolation.variance <-
fit.nls(lambda, 1:NCOL(variance.jackknife), variance.jackknife)
)
variance.jackknife2 <- vector("numeric",ncoef^2)
switch(jackknife.estimation,
"nonl" = for(i in 1:NCOL(variance.jackknife))
variance.jackknife2[i] <- predict(extrapolation.variance[[i]],
newdata = data.frame(lambda = -1)),
"quad" = variance.jackknife2 <- predict(extrapolation.variance,
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
}
if (asymptotic) {
c11 <- cov(PSI)
a11 <- diag.block(am)
a11.inv <- solve(a11)
sigma <- a11.inv %*% c11 %*% t(a11.inv)
s <- construct.s(ncoef, lambda, fitting.method, extrapolation)
d.inv <- solve(s %*% t(s))
sigma.gamma <- d.inv %*% s %*% sigma %*% t(s) %*% d.inv
g <- list()
switch(fitting.method,
"quad" = g <- c(1, -1, 1),
"line" = g <- c(1, -1),
"logl" = for(i in 1:ncoef)
g[[i]] <-
c(exp(coef(extrapolation)[1, i] - coef(extrapolation)[2, i]),
-exp(coef(extrapolation)[1, i] - coef(extrapolation)[2,i])),
"log2" = for(i in 1:ncoef)
g[[i]] <-
c(exp(coef(extrapolation[[i]])[1] - coef(extrapolation[[i]])[2]),
-exp(coef(extrapolation[[i]])[1] - coef(extrapolation[[i]])[2])),
"nonl" = for(i in 1:ncoef)
g[[i]] <-
c(-1, -(coef(extrapolation[[i]])[3] - 1)^-1,
coef(extrapolation[[i]])[2] / (coef(extrapolation[[i]])[3] - 1)^2)
)
g <- diag.block(g, ncoef)
variance.asymptotic <- (t(g) %*% sigma.gamma %*% g) / ndes
dimnames(variance.asymptotic) <- list(p.names, p.names)
}
# creating class "mcsimex"
theta <- matrix(unlist(theta.all), nrow = B)
theta.all <- list()
for (i in 1:ncoef)
theta.all[[p.names[i]]] <-
data.frame(theta[, seq(i, ncoef * nlambda, by = ncoef)])
z <- cbind(lambda, estimates)
z <- rbind(c(-1, SIMEX.estimate), z) # returning the estimated values
colnames(z) <- c("lambda", names(model$coefficients))
erg <- list(
coefficients = z[1, -1],  # SIMEX corrected Coefficients
SIMEX.estimates = z,  # all thetas as a matrix
lambda = lambda,   # vector for the values for lambda
model = model,  # the naive model
mc.matrix = mc.matrix,  # vector of values of measurement.error
B = B,  # number of Simulations
extrapolation = extrapolation,  # model of the extrapolation
fitting.method = fitting.method,# which fitting method was used
SIMEXvariable = SIMEXvariable, # string containing the colnames SIMEXvariables
call = cl,
theta = theta.all
)
class(erg) <- ("mcsimex")
fitted.values <- predict(erg, newdata = model$model[, -1, drop = FALSE],
type = "response")
erg$fitted.values <- fitted.values
if (is.factor(model$model[, 1]))
erg$residuals <- as.numeric(levels(model$model[, 1]))[model$model[, 1]] -
fitted.values else erg$residuals <- model$model[, 1] - fitted.values

if (jackknife.estimation != FALSE) {
erg$extrapolation.variance <- extrapolation.variance
erg$variance.jackknife <- variance.jackknife
erg$variance.jackknife.lambda <- variance.jackknife.lambda
}
if (asymptotic) {
erg$PSI <- PSI
erg$c11 <- c11
erg$a11 <- a11
erg$sigma <- sigma
erg$sigma.gamma <- sigma.gamma
erg$g <- g
erg$s <- s
erg$variance.asymptotic <- variance.asymptotic
}
return(erg)
}

