

irls_glogit <- function(formula, data,
                        tol = 0.000001, offset = 0, m = 1 ) {

## Set up the model components
   mf <- model.frame(formula, data)
   y <- model.response(mf, "numeric")
   X <- model.matrix(formula, data = data)

## Check for missing values; stop if any.
   if (any(is.na(cbind(y, X)))) stop("Some data are missing.")

## Check for nonsensical values; stop if any.
   if (any(y < 0 | y > m)) stop("Some data are absurd.")

## Initialize mu, eta, the deviance, etc.
   mu <- rep(mean(y), length(y))
   eta <- log(mu/(m - mu))
   dev <- 2 * sum(y * log(pmax(y,1)/mu) +
              (m - y)* log(pmax(y,1)/(m - mu)) )
   deltad <- 1
   i <- 1

## Loop through the IRLS algorithm
   while (abs(deltad) > tol ) {
     w <- mu*(m-mu)/m
     z <- eta + (y-mu)/w - offset
     mod <- lm(z ~ X-1, weights=w)
     eta <- mod$fit + offset
     mu <- m/(1+exp(-eta))
     dev.old <- dev
     dev <- 2 * sum(y*log(pmax(y,1)/mu) +
                    (m - y)* log(pmax(y,1)/(m - mu)))
     deltad <- dev - dev.old
     i <- i + 1
   }

## Build some post-estimation statistics
   df.residual <- summary(mod)$df[2]
   pearson.chi2 <- sum((y - mu)^2 /
                       (mu * (1 - mu/m))) / df.residual

## Return a brief result
   return(list(coef = coef(mod),
               se = sqrt(diag(summary(mod)$cov.unscaled)),
               pearson.chi2 = pearson.chi2))
 }

