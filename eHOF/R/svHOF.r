"svHOF" <- function (
		x, y, M, model, mod, famname = binomial, ...) {
#x <- scale01(x)
if(model == 'I') a = log((1 - mean(y/M))/mean(y/M))
if(model == 'II') {
        m1 <- glm(cbind(y, M - y) ~ x, family=famname)
        k <- coef(m1)
        names(k) <- NULL
        a <- -k[1]
        b <- -k[2]
      }
if(model == 'III') { 
      a <- mod$par[[1]]
      b <- mod$par[[2]]
      c <- 0
      }

if(model %in% c('IV','VI','VII')) {
        m2 <- glm(cbind(y, M - y) ~ x + I(x^2), family=famname)
        k <- coef(m2)
        if (k[3] > 0) {
          m2 <- glm(cbind(y, M-y) ~ x + offset(-x^2), family=famname)
          k <- c(coef(m2), -1)
        }
        names(k) <- NULL
        u <- -k[2]/2/k[3]
        h <- plogis(k[1] - k[2]^2/4/k[3])
        h <- min(0.98, h)
        u <- min(0.98, u)
        u <- max(0.02, u)
        r <- 1/h * (-2 * h + 2 * sqrt(h))/2
        r <- log(r)
        b <- 5.07 - 0.227 * k[3]
        a <- -b * u + r
        c <- -a
      }

if(model %in% c('V')) {
     a <- mod$par[[1]]
     b <- mod$par[[2]]
     c <- mod$par[[3]]
     d <- mod$par[[2]]
    }
  
  out <- switch(model,
	    I = list(a = a),
	   II = list(a = a, b = b),
	  III = list(a = a, b = b, c = c),
	   IV = list(a = a, b = b, c = c),
	    V = list(a = a, b = b, c = c, d = b),
	   VI = list(a = a, b = b, c = c, d = diff(range(x))/2),
    VII = list(a = a, b = b, c = c, d = diff(range(x))/2, e = diff(range(x))/2)
        )
  unlist(out)
}
