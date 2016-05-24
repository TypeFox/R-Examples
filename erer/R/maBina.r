maBina <- function(w, x.mean = TRUE, rev.dum = TRUE, digits = 3, 
  subset.name = NULL, subset.value)
{
  # 1. Check inputs
  if (!inherits(w, "glm")) {stop("Need an object from 'glm()'.\n")}
  link <- w$family$link
  if (link != "probit" & link != "logit") {
    stop("Need a binary probit or logit model.\n")}
  if (is.null(dim(w$x))) stop("Please specify 'x = TRUE' in glm().\n")

  # 2. Calcualate marginal effects (ME)        
  x <- as.matrix(w$x)
  if(!is.null(subset.name)) x <- x[x[, subset.name] == subset.value, ] 
  x.bar <- as.matrix(colMeans(x))

  b.est <- as.matrix(coef(w)); K <- nrow(b.est)
  xb <- t(x.bar) %*% b.est
  pfun <- switch(link, probit = pnorm, logit = plogis)
  dfun <- switch(link, probit = dnorm, logit = dlogis)  
  if (x.mean) {f.xb <- dfun(xb)} else {f.xb <- mean(dfun(x %*% b.est))}
  me <- f.xb * coef(w)
  
  # 3. Standard errors for ME
  if (link == "probit") {s <- -xb} else {s <- 1 - 2 * pfun(xb)}
  dr <- c(f.xb) * (diag(1, K, K) + c(s) * (b.est %*% t(x.bar))) 
  va <- dr %*% vcov(w) %*% t(dr) 
  se <- sqrt(diag(va))  

  # 4. Revise ME and error for dummy variable
  if (rev.dum) {
    for (i in 1:ncol(x)) {
      if (identical(sort(unique(x[,i])), c(0, 1))) {
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
  
  # 5. Output
  out <- data.frame(effect = me, error = se)
  out$t.value <- out$effect / out$error
  out$p.value <- 2 * (1- pt(abs(out[, 3]), nrow(x) - K))
  out <- round(out, digits = digits)

  result <- listn(link, f.xb, w, x, out)
  class(result) <- "maBina"
  return(result)
}

print.maBina <- function(x, ...) {print(x$out)}