`metacor.DSL` <-
function(r, n, labels, alpha = 0.05, plot = TRUE, xlim = c(-1, 1), transform = TRUE){
z <- 0.5 * log( (1 + r) / (1 - r))
z.var <- 1 / (n - 3)
z.lower<-z-qnorm(alpha/2)*sqrt(z.var)
z.upper<-z+qnorm(alpha/2)*sqrt(z.var)
fixed.z.mean <- sum( (n - 3) * z) / sum(n - 3)
Q <- sum( (n - 3) * (z - fixed.z.mean)^2)
k <- length(r)
a <- sum(n - 3) - (sum((n - 3)^2) / sum(n - 3))
tau <- (Q - (k - 1)) / a
if (tau < 0) tau <- 0
w <- ((1 / (n - 3)) + tau)^ - 1
z.mean <- sum(w * z) / sum(w)
z.mean.se <- sqrt(1 / sum(w))
g <- z.mean / z.mean.se
z.mean.lower <- z.mean - qnorm(alpha / 2) * z.mean.se
z.mean.upper <- z.mean + qnorm(alpha / 2) * z.mean.se
z.to.r <- function(z) (exp(2 * z) - 1) / (exp(2 * z) + 1)
r.lower <- z.to.r(z.lower)
r.upper <- z.to.r(z.upper)
r.se <- ((r.lower - r) / (-qnorm(alpha / 2)))
r.mean <- z.to.r(z.mean)
r.mean.lower <- z.to.r(z.mean.lower)
r.mean.upper <- z.to.r(z.mean.upper)
r.mean.se <- ((r.mean.lower - r.mean) / (-qnorm(alpha/2)))
p <- pnorm(abs(g), lower.tail = F)
if (plot){
  if (transform) metaplot(r, r.se, labels = labels, xlab = quote("Correlation coefficient"~italic(r)), ylab = "" , summn = r.mean, sumse = r.mean.se, sumnn = r.mean.se^-2, xlim = xlim)
  else metaplot(z, sqrt(z.var), labels = labels, xlab = quote("Fisher"~italic(z)), ylab = "", summn = z.mean, sumse = z.mean.se, sumnn = z.mean.se^-2, xlim = xlim)

  }
res<-list()
res$z <- z
res$z.var <- z.var
res$z.lower <- z.lower
res$r.lower <- r.lower
res$z.upper <- z.upper
res$r.upper <- r.upper

res$z.mean <- z.mean
res$r.mean <- r.mean
res$z.mean.se <- z.mean.se
res$z.mean.lower <- z.mean.lower
res$r.mean.lower <- r.mean.lower
res$z.mean.upper <- z.mean.upper
res$r.mean.upper <- r.mean.upper
res$p <- p
invisible(res)


}

