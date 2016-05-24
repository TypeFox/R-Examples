pkgname <- "simex"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('simex')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("build.mc.matrix")
### * build.mc.matrix

flush(stderr()); flush(stdout())

### Name: mc.matrix
### Title: Build and check misclassification matrices from empirical
###   estimations
### Aliases: build.mc.matrix check.mc.matrix
### Keywords: regression

### ** Examples

Pi <- matrix(data = c(0.989, 0.01, 0.001, 0.17, 0.829, 0.001, 0.001, 0.18, 0.819),
    nrow = 3, byrow = FALSE)
check.mc.matrix(list(Pi))
check.mc.matrix(list(build.mc.matrix(Pi)))
build.mc.matrix(Pi)

Pi3 <- matrix(c(0.8, 0.2, 0, 0, 0, 0.8, 0.1, 0.1, 0, 0.1, 0.8, 0.1, 0, 0, 0.3, 0.7),
    nrow = 4)
check.mc.matrix(list(Pi3))
build.mc.matrix(Pi3)
check.mc.matrix(list(build.mc.matrix(Pi3)))

P1 <- matrix(c(1, 0, 0, 1), nrow = 2)
P2 <- matrix(c(0.8, 0.15, 0, 0.2, 0.7, 0.2, 0, 0.15, 0.8), nrow = 3, byrow = TRUE)
P3 <- matrix(c(0.4, 0.6, 0.6, 0.4), nrow = 2)
mc.matrix <- list(P1, P2, P3)
check.mc.matrix(mc.matrix) # TRUE FALSE FALSE



cleanEx()
nameEx("diag.block")
### * diag.block

flush(stderr()); flush(stdout())

### Name: diag.block
### Title: Constructs a block diagonal matrix
### Aliases: diag.block
### Keywords: datagen

### ** Examples

a <- matrix(rep(1, 4), nrow = 2)
b <- matrix(rep(2, 6), nrow = 2)
e <- c(3, 3, 3, 3)
f <- t(e)
d <- list(a, b, e, f)
diag.block(d)
diag.block(a, 3)



cleanEx()
nameEx("mcsimex")
### * mcsimex

flush(stderr()); flush(stdout())

### Name: mcsimex
### Title: Misclassification in models using MCSIMEX
### Aliases: mcsimex print.mcsimex summary.mcsimex print.summary.mcsimex
###   plot.mcsimex predict.mcsimex refit.mcsimex
### Keywords: models

### ** Examples

x <- rnorm(200, 0, 1.142)
z <- rnorm(200, 0, 2)
y <- factor(rbinom(200, 1, (1 / (1 + exp(-1 * (-2 + 1.5 * x -0.5 * z))))))
Pi <- matrix(data = c(0.9, 0.1, 0.3, 0.7), nrow = 2, byrow = FALSE)
dimnames(Pi) <- list(levels(y), levels(y))
ystar <- misclass(data.frame(y), list(y = Pi), k = 1)[, 1]
naive.model <- glm(ystar ~ x + z, family = binomial, x = TRUE, y = TRUE)
true.model  <- glm(y ~ x + z, family = binomial)
simex.model <- mcsimex(naive.model, mc.matrix = Pi, SIMEXvariable = "ystar")

op <- par(mfrow = c(2, 3))
invisible(lapply(simex.model$theta, boxplot, notch = TRUE, outline = FALSE,
    names = c(0.5, 1, 1.5, 2)))
plot(simex.model)

simex.model2 <- refit(simex.model, "line")
plot(simex.model2)
par(op)

# example for a function which can be supplied to the function mcsimex()
# "ystar" is the variable which is to be misclassified
# using the example above
## Not run: 
##D my.misclass <- function (datas, k) {
##D 	ystar <- datas$"ystar"
##D 	p1 <- matrix(data = c(0.75, 0.25, 0.25, 0.75), nrow = 2, byrow = FALSE)
##D 	colnames(p1) <- levels(ystar)
##D 	rownames(p1) <- levels(ystar)
##D 	p0 <- matrix(data = c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = FALSE)
##D 	colnames(p0) <- levels(ystar)
##D 	rownames(p0) <- levels(ystar)
##D 	ystar[datas$x < 0] <-
##D 	    misclass(data.frame(ystar = ystar[datas$x < 0]), list(ystar = p1), k = k)[, 1]
##D 	ystar[datas$x > 0] <-
##D 	    misclass(data.frame(ystar = ystar[datas$x > 0]), list(ystar = p0), k = k)[, 1]
##D 	ystar <- factor(ystar)
##D 	return(data.frame(ystar))
##D }
##D 
##D simex.model.differential <- mcsimex(naive.model, mc.matrix = "my.misclass", SIMEXvariable = "ystar")
## End(Not run)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("misclass")
### * misclass

flush(stderr()); flush(stdout())

### Name: misclass
### Title: Generates misclassified data
### Aliases: misclass
### Keywords: datagen

### ** Examples

x1 <- factor(rbinom(100, 1, 0.5))
x2 <- factor(rbinom(100, 2, 0.5))

p1 <- matrix(c(1, 0, 0, 1), nrow = 2)
p2 <- matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8), nrow = 3)

colnames(p1) <- levels(x1)
colnames(p2) <- levels(x2)

x <- data.frame(x1 = x1, x2 = x2)
mc.matrix <- list(x1 = p1, x2 = p2)

x.mc <- misclass(data.org = x, mc.matrix = mc.matrix, k = 1)

identical(x[, 1], x.mc[, 1]) # TRUE
identical(x[, 2], x.mc[, 2]) # FALSE



cleanEx()
nameEx("simex-package")
### * simex-package

flush(stderr()); flush(stdout())

### Name: simex-package
### Title: Error or misclassification correction in models using (MC)SIMEX
### Aliases: simex-package
### Keywords: package models regression

### ** Examples

# See example(simex) and example(mcsimex)
## Seed
set.seed(49494)

## simulating the measurement error standard deviations
sd_me1 <- 0.3
sd_me2 <- 0.4
temp <- runif(100, min = 0, max = 0.6)
sd_me_het1 <- sort(temp)
temp2 <- rnorm(100, sd = 0.1)
sd_me_het2 <- abs(sd_me_het1 + temp2)

## simulating the independent variables x (real and with measurement error):
x_real1 <- rnorm(100)
x_real2 <- rpois(100, lambda = 2)
x_real3 <- -4*x_real1 + runif(100, min = -2, max = 2)  # correlated to x_real

x_measured1 <- x_real1 + sd_me1 * rnorm(100)
x_measured2 <- x_real2 + sd_me2 * rnorm(100)
x_het1 <- x_real1 + sd_me_het1 * rnorm(100)
x_het2 <- x_real3 + sd_me_het2 * rnorm(100)

## calculating dependent variable y:
y1  <- x_real1 + rnorm(100, sd = 0.05)
y2 <- x_real1 + 2*x_real2 + rnorm(100, sd = 0.08)
y3 <- x_real1 + 2*x_real3 + rnorm(100, sd = 0.08)


### one variable with homoscedastic measurement error
(model_real <- lm(y1  ~ x_real1))

(model_naiv <- lm(y1  ~ x_measured1, x = TRUE))

(model_simex <- simex(model_naiv, SIMEXvariable = "x_measured1", measurement.error = sd_me1))
plot(model_simex)


### two variables with homoscedastic measurement errors
(model_real2 <- lm(y2 ~ x_real1 + x_real2))

(model_naiv2 <- lm(y2 ~ x_measured1 + x_measured2, x = TRUE))

(model_simex2 <- simex(model_naiv2, SIMEXvariable = c("x_measured1", "x_measured2"), 
                       measurement.error = cbind(sd_me1, sd_me2)))

plot(model_simex2)


### one variable with increasing heteroscedastic measurement error
model_real 

(mod_naiv1 <- lm(y1  ~ x_het1, x = TRUE))

(mod_simex1 <- simex(mod_naiv1, SIMEXvariable = "x_het1", measurement.error = sd_me_het1, asymptotic = FALSE))

plot(mod_simex1)

## Not run: 
##D ### two correlated variables with heteroscedastic measurement errors
##D (model_real3 <- lm(y3 ~ x_real1 + x_real3))
##D 
##D (mod_naiv2 <- lm(y3 ~ x_het1 + x_het2, x = TRUE))
##D 
##D (mod_simex2 <- simex(mod_naiv2, SIMEXvariable = c("x_het1", "x_het2"), 
##D                      measurement.error = cbind(sd_me_het1, sd_me_het2), asymptotic = FALSE))
##D plot(mod_simex2)
##D 
##D 
##D ### two variables, one with homoscedastic, one with heteroscedastic measurement error
##D model_real2
##D 
##D (mod_naiv3 <- lm(y2 ~ x_measured1 + x_het2, x = TRUE))
##D 
##D (mod_simex3 <- simex(mod_naiv3, SIMEXvariable = c("x_measured1", "x_het2"), 
##D                      measurement.error = cbind(sd_me1, sd_me_het2), asymptotic = FALSE))
##D 
##D 
##D ### glm: two variables, one with homoscedastic, one with heteroscedastic measurement error
##D t <- x_real1 + 2*x_real2 
##D g <- 1 / (1 + exp(-t))
##D u <- runif(100)
##D ybin <- as.numeric(u < g)
##D 
##D 
##D (logit_real <- glm(ybin ~ x_real1 + x_real2, family = binomial))
##D 
##D (logit_naiv <- glm(ybin ~ x_measured1 + x_het2, x = TRUE, family = binomial))
##D 
##D (logit_simex <- simex(logit_naiv, SIMEXvariable = c("x_measured1", "x_het2"), 
##D                       measurement.error = cbind(sd_me1, sd_me_het2), asymptotic = FALSE))
##D summary(logit_simex)
##D print(logit_simex)
##D plot(logit_simex)
## End(Not run)



cleanEx()
nameEx("simex")
### * simex

flush(stderr()); flush(stdout())

### Name: simex
### Title: Measurement error in models using SIMEX
### Aliases: simex print.simex summary.simex print.summary.simex plot.simex
###   predict.simex refit refit.simex
### Keywords: models

### ** Examples

## Seed
set.seed(49494)

## simulating the measurement error standard deviations
sd_me <- 0.3
sd_me2 <- 0.4
temp <- runif(100, min = 0, max = 0.6)
sd_me_het1 <- sort(temp)
temp2 <- rnorm(100, sd = 0.1)
sd_me_het2 <- abs(sd_me_het1 + temp2)

## simulating the independent variables x (real and with measurement error):
x_real <- rnorm(100)
x_real2 <- rpois(100, lambda = 2)
x_real3 <- -4*x_real + runif(100, min = -10, max = 10)  # correlated to x_real

x_measured <- x_real + sd_me * rnorm(100)
x_measured2 <- x_real2 + sd_me2 * rnorm(100)
x_het1 <- x_real + sd_me_het1 * rnorm(100)
x_het2 <- x_real3 + sd_me_het2 * rnorm(100)

## calculating dependent variable y:
y <- x_real + rnorm(100, sd = 0.05)
y2 <- x_real + 2*x_real2 + rnorm(100, sd = 0.08)
y3 <- x_real + 2*x_real3 + rnorm(100, sd = 0.08)


### one variable with homoscedastic measurement error
(model_real <- lm(y ~ x_real))

(model_naiv <- lm(y ~ x_measured, x = TRUE))

(model_simex <- simex(model_naiv, SIMEXvariable = "x_measured", measurement.error = sd_me))
plot(model_simex)


### two variables with homoscedastic measurement errors
(model_real2 <- lm(y2 ~ x_real + x_real2))

(model_naiv2 <- lm(y2 ~ x_measured + x_measured2, x = TRUE))

(model_simex2 <- simex(model_naiv2, SIMEXvariable = c("x_measured", "x_measured2"), 
                           measurement.error = cbind(sd_me, sd_me2)))

plot(model_simex2)


### one variable with increasing heteroscedastic measurement error
model_real 

(mod_naiv1 <- lm(y ~ x_het1, x = TRUE))

(mod_simex1 <- simex(mod_naiv1, SIMEXvariable = "x_het1", measurement.error = sd_me_het1, asymptotic = FALSE))

plot(mod_simex1)

## Not run: 
##D ### two correlated variables with heteroscedastic measurement errors
##D (model_real3 <- lm(y3 ~ x_real + x_real3))
##D 
##D (mod_naiv2 <- lm(y3 ~ x_het1 + x_het2, x = TRUE))
##D 
##D (mod_simex2 <- simex(mod_naiv2, SIMEXvariable = c("x_het1", "x_het2"), 
##D                          measurement.error = cbind(sd_me_het1, sd_me_het2), asymptotic = FALSE))
##D plot(mod_simex2)
##D 
##D 
##D ### two variables, one with homoscedastic, one with heteroscedastic measurement error
##D model_real2
##D 
##D (mod_naiv3 <- lm(y2 ~ x_measured + x_het2, x = TRUE))
##D 
##D (mod_simex3 <- simex(mod_naiv3, SIMEXvariable = c("x_measured", "x_het2"), 
##D                          measurement.error = cbind(sd_me, sd_me_het2), asymptotic = FALSE))
##D 
##D 
##D ### glm: two variables, one with homoscedastic, one with heteroscedastic measurement error
##D t <- x_real + 2*x_real2 + rnorm(100, sd = 0.01)
##D g <- 1 / (1 + exp(t))
##D u <- runif(100)
##D ybin <- as.numeric(u < g)
##D 
##D 
##D (logit_real <- glm(ybin ~ x_real + x_real2, family = binomial))
##D 
##D (logit_naiv <- glm(ybin ~ x_measured + x_het2, x = TRUE, family = binomial))
##D 
##D (logit_simex <- simex(logit_naiv, SIMEXvariable = c("x_measured", "x_het2"), 
##D                           measurement.error = cbind(sd_me, sd_me_het2), asymptotic = FALSE))
##D summary(logit_simex)
##D print(logit_simex)
##D plot(logit_simex)
## End(Not run)



### * <FOOTER>
###
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
