library(dglm)

n <- 400

TestFunc <- function(mean.formula, var.formula, model.df) {
  
  my.glm <- glm(formula = mean.formula,
                data = model.df)
  
  x <- predict(my.glm, se.fit = TRUE)
  
  my.dglm <- dglm(formula = mean.formula,
                  dformula = var.formula,
                  data = model.df)
  
  y <- predict(my.dglm, se.fit = TRUE)
  
  z <- predict(my.dglm$dispersion.fit, se.fit = TRUE)
  
  return(list(glm = x, dglm.mean = y, dglm.var = z))
}


a <- runif(n)
b <- runif(n)
c <- runif(n)
d <- runif(n)

betas <- c(runif(n = 3, min = 10, max = 20), runif(n = 3, min = 1, max = 2))

mu <- betas[1] + a*betas[2] + b*betas[3]
sd <- sqrt(exp(betas[4] + c*betas[5] + d*betas[6]))
y <- mu + rnorm(n = n, sd = sd)

my.df <- data.frame(y, a, b, c, d)

mean.formula <- as.formula('y ~ a + b')
var.formula <- as.formula('~ c + d')

l <- TestFunc(mean.formula = mean.formula, var.formula = var.formula, model.df = my.df)

par(mfrow = c(1, 2))
plot(mu, l$dglm.mean$fit); abline(0, 1)
legend(x = 'topleft', legend = paste(c('intercept  ', '           a', '          b'), round(betas[1:3], 1)))
plot(sd, exp(l$dglm.var$fit/l$dglm.var$residual.scale)); abline(0, 1)
legend(x = 'topleft', legend = paste(c('intercept  ', '           c', '           d'), round(betas[4:6], 2)))

my.dglm <- dglm(formula = mean.formula, dformula = var.formula, data = my.df)
