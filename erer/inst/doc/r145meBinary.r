# A. Load data; estimate a binary logit model; marginal effect by maBina()
library(erer); data(daIns)
ra <- glm(formula = Y ~ Injury + HuntYrs + Nonres + Lspman + Lnong + 
                    Gender + Age + Race + Marital + Edu + Inc + TownPop,
          family = binomial(link = "logit"), data = daIns, x = TRUE)
(ca <- data.frame(summary(ra)$coefficients))
(magin <- maBina(w = ra))

# B. Calculate marginal effects (ME) and standard errors
# B1. Input data
w <- ra; x.mean <- TRUE; rev.dum <- TRUE

# B2. Calcualate ME
x <- as.matrix(w$x)
x.bar <- as.matrix(colMeans(x))
b.est <- as.matrix(coef(w)); K <- nrow(b.est)
xb <- t(x.bar) %*% b.est
link <- w$family$link
pfun <- switch(link, probit = pnorm, logit = plogis)
dfun <- switch(link, probit = dnorm, logit = dlogis)  
if (x.mean) {f.xb <- dfun(xb)} else {f.xb <- mean(dfun(x %*% b.est))}
me <- f.xb * coef(w)

# B3. Standard errors for ME
if (link == "probit") {s <- -xb} else {s <- 1 - 2 * pfun(xb)}
dr <- c(f.xb) * (diag(1, K, K) + c(s) * (b.est %*% t(x.bar))) 
va <- dr %*% vcov(w) %*% t(dr) 
se <- sqrt(diag(va))

# B4. Revise ME and error for dummy independent variables
if (rev.dum) {
  for (i in 1:ncol(x)) {
    if (identical(sort(unique(x[, i])), c(0, 1))) {
      x.d1 <- x.bar; x.d1[i, 1] <- 1
      x.d0 <- x.bar; x.d0[i, 1] <- 0
      me[i] <- pfun(t(x.d1) %*% b.est) - 
               pfun(t(x.d0) %*% b.est)
      dr2 <- dfun(t(x.d1) %*% b.est) %*%  t(x.d1) -  
             dfun(t(x.d0) %*% b.est) %*%  t(x.d0)           
      se[i] <- sqrt(c(dr2 %*% vcov(w) %*% t(dr2)))
    }
  }
}

# C. Output
out <- data.frame(effect = me, error = se)
out$t.value <- out$effect / out$error
out$p.value <- 2 * (1- pt(abs(out[, 3]), nrow(x) - K))
out <- round(out, digits = 3); out