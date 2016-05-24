library(grplasso)

data(splice)

tol <- 5 * 10^-8
control <- grpl.control(tol = tol, trace = 0)

##############################################
##                                          ##
## See whether correct solution is obtained ##
##                                          ##
##############################################

fit <- grplasso(y ~ Pos.1 * Pos.2, data = splice, lambda = 25, control = control,
                center = FALSE, contrast = list(Pos.1 = "contr.sum",
                                  Pos.2 = "contr.sum"))

sol <- structure(c(-0.13977233, 0.0226494585308358, 0.0897927902302861, 
-0.0311859296853295, 0.519851153317134, -0.25546398114419, -0.209672224889056, 
0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(16L, 1L), .Dimnames = list(
    c("(Intercept)", "Pos.11", "Pos.12", "Pos.13", "Pos.21", 
    "Pos.22", "Pos.23", "Pos.11:Pos.21", "Pos.12:Pos.21", "Pos.13:Pos.21", 
    "Pos.11:Pos.22", "Pos.12:Pos.22", "Pos.13:Pos.22", "Pos.11:Pos.23", 
    "Pos.12:Pos.23", "Pos.13:Pos.23"), "25"))

stopifnot(all.equal(coef(fit), sol, tol = 10^-6))

##################################################################
##                                                              ##
## Check whether different contrasts lead to the same solutions ##
## I.e. check whether the (back-) transformations work right    ##
##                                                              ##
##################################################################

contr.A <- list(Pos.1 = "contr.sum",
                Pos.2 = "contr.sum",
                Pos.3 = "contr.sum")
contr.B <- list(Pos.1 = "contr.helmert",
                Pos.2 = "contr.helmert",
                Pos.3 = "contr.helmert")

fit.A <- grplasso(y ~ Pos.1 * Pos.2 * Pos.3, nonpen = ~ 1, data = splice,
                  standardize = TRUE, lambda = 61:1,
                  contrasts = contr.A, control = control)
fit.B <- grplasso(y ~ Pos.1 * Pos.2 * Pos.3, nonpen = ~ 1, data = splice,
                  standardize = TRUE, lambda = 61:1,
                  contrasts = contr.B, control = control)
#par(mfrow = c(1, 2))
#plot(fit.A, log = "x")
#plot(fit.B, log = "x")

if(max(abs(fit.A$fn.val - fit.B$fn.val) / fit.A$fn.val) > tol)
  stop("Inconsistent result when changing the encoding scheme (fn.val)")

pred.A <- predict(fit.A, newdata = splice, type = "response")
pred.B <- predict(fit.B, newdata = splice, type = "response")
m      <- abs(pred.A - pred.B)

if(max(m) > sqrt(tol))
  stop("Inconsistent result when changing the encoding scheme (prediction)")

#range(pred.A - pred.B)
#range((pred.A - pred.B) / (1 + pred.A))

###############################################
##                                           ##
## Check whether offset is working correctly ##
##                                           ##
###############################################

## Fit an ordinary model, with unpenalized intercept
lambda.max <- lambdamax(y ~ ., data = splice,
                        model = LogReg(), standardize = TRUE)

fit1 <- grplasso(y ~ ., data = splice, lambda = 0.5 * lambda.max,
                 standardize = TRUE, control = control)

## Plugging in the intercept as an offset should lead to an intercept
## which is close to zero
shift <- 4
intercept <- rep(shift, nrow(splice))
fit2      <- grplasso(y ~ . + offset(intercept), lambda = 0.5 * lambda.max,
                      data = splice, standardize = TRUE, control = control)

d.coef    <- coef(fit2) - coef(fit1)
d.coef[1] <- d.coef[1] + shift

if(max(abs(d.coef) / (1 + abs(coef(fit1)))) > 5 * 10^-4)
  stop("Inconsistent result when using offset (d.coef)")

pred.new <- predict(fit2, newdata = splice) ## didn't work in earlier version
pred.in  <- predict(fit2)

stopifnot(all.equal(pred.new, pred.in))

#################################################
##                                             ##
## Check whether max.iter is working correctly ##
##                                             ##
#################################################

fit.maxiter <- grplasso(y ~ Pos.1 * Pos.2, data = splice, lambda = c(1, 0.1),
                        control = grpl.control(max.iter = 2, trace = 2,
                          inner.loops = 0),
                        contrast = list(Pos.1 = "contr.sum",
                          Pos.2 = "contr.sum"))

stopifnot(all(!fit.maxiter$converged))
