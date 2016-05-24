comb.samples <- 
function (svy.A, svy.B, svy.C = NULL, y.lab, z.lab, form.x, estimation = NULL, 
          micro = FALSE, ...) 
{
  data.A <- svy.A$variables
  if(is.factor(data.A[, y.lab])){
    y.lev <- levels(data.A[, y.lab])
    levels(data.A[, y.lab]) <- 1:nlevels(data.A[, y.lab])
    svy.A$variables <- data.A
  }
  n.A <- nrow(data.A)
  w.A <- weights(svy.A)

  data.B <- svy.B$variables
  if(is.factor(data.B[, z.lab])){
    z.lev <- levels(data.B[, z.lab])
    levels(data.B[, z.lab]) <- 1:nlevels(data.B[, z.lab])
    svy.B$variables <- data.B
  }
  n.B <- nrow(data.B)
  w.B <- weights(svy.B)

  X.A <- model.matrix(form.x, data = data.A)
  X.B <- model.matrix(form.x, data = data.B)
  form.y <- as.formula(paste("~", y.lab, "-1", collapse = ""))
  form.z <- as.formula(paste("~", z.lab, "-1", collapse = ""))
  Y.A <- model.matrix(form.y, data = data.A)
  Z.B <- model.matrix(form.z, data = data.B)

# estimates B_yx
  QR.A <- qr(X.A * sqrt(w.A))
  beta.yx.A <- qr.coef(QR.A, Y.A * sqrt(w.A))
  beta.yx.A[is.na(beta.yx.A)] <- 0

# estimates B_zx
  QR.B <- qr(X.B * sqrt(w.B))
  beta.zx.B <- qr.coef(QR.B, Z.B * sqrt(w.B))
  beta.zx.B[is.na(beta.zx.B)] <- 0
#
# derives X'X  
  XX.wA <- t(X.A) %*% (X.A * w.A)
  XX.wB <- t(X.B) %*% (X.B * w.B)
  gamma.p <- n.A/(n.A + n.B)
  XX.pool <- gamma.p * XX.wA + (1 - gamma.p) * XX.wB
#
# CIA estimate  
  YZ.CIA <- t(beta.yx.A) %*% XX.pool %*% beta.zx.B
  if(is.factor(data.A[, y.lab])) rownames(YZ.CIA) <- y.lev
  if(is.factor(data.B[, z.lab])) colnames(YZ.CIA) <- z.lev   
#  dimnames(YZ.CIA) <- list(y.lev, z.lev)
  out <- list(yz.CIA = YZ.CIA, call = match.call())
  if (!is.null(svy.C)) {
    data.C <- svy.C$variables
    if(is.factor(data.C[, y.lab])){
      y.lev.C <- levels(data.C[, y.lab])
      if (all.equal(y.lev, y.lev.C)) 
        levels(data.C[, y.lab]) <- 1:nlevels(data.C[, y.lab])
      else stop("The levels of y.lab in svy.A and in svy.C do not match")
    }    
    if(is.factor(data.C[, z.lab])) {
      z.lev.C <- levels(data.C[, z.lab])
      if (all.equal(z.lev, z.lev.C)) 
        levels(data.C[, z.lab]) <- 1:nlevels(data.C[, z.lab])
      else stop("The levels of z.lab in svy.B and in svy.C do not match") 
    }
    svy.C$variables <- data.C
    n.C <- nrow(data.C)
    w.C <- weights(svy.C)
    if (estimation == "ITWS" || estimation == "i2ws" || estimation == 
        "incomplete") {
      tot.y.A <- colSums(Y.A * w.A)
      tot.z.B <- colSums(Z.B * w.B)
      tot.yz <- c(tot.y.A, tot.z.B[-1])
      form.yz <- as.formula(paste("~", paste(y.lab, z.lab, 
                                             sep = "+"), "- 1", sep = ""))
      cal.C <- calibrate(design = svy.C, formula = form.yz, 
                         population = tot.yz, ...)
    }
    if (estimation == "STWS" || estimation == "s2ws" || estimation == 
        "synthetic") {
      X.C <- model.matrix(form.x, data = data.C)
      Y.C <- model.matrix(form.y, data = data.C)
      Z.C <- model.matrix(form.z, data = data.C)
      resY.C <- Y.C - (X.C %*% beta.yx.A)
      resZ.C <- Z.C - (X.C %*% beta.zx.B)
      c.y <- ncol(Y.C)
      c.z <- ncol(Z.C)
      new.YZ <- matrix(NA, nrow = n.C, ncol = (c.y * c.z))
      for (i in 1:n.C) {
        m1 <- cbind(Y.C[i, ]) %*% rbind(Z.C[i, ])
        m2 <- cbind(resY.C[i, ]) %*% rbind(resZ.C[i, ])
        new.YZ[i, ] <- c(m1) - c(m2)
      }
      lab1 <- rep(colnames(Y.C), c.z)
      lab2 <- rep(colnames(Z.C), each = c.y)
      lab <- paste(lab1, lab2, sep = "_")
      colnames(new.YZ) <- lab
      orig.vars <- colnames(svy.C$variables)
      svy.C$variables <- data.frame(svy.C$variables, new.YZ)
      vec.tot <- c(YZ.CIA)
      names(vec.tot) <- lab
      form.yz <- as.formula(paste("~", paste(lab, collapse = "+"), 
                                  "- 1", sep = ""))
      cal.C <- calibrate(design = svy.C, formula = form.yz, 
                         population = vec.tot, ...)
      cal.C$variables <- cal.C$variables[, orig.vars]
    }
    ww.C <- weights(cal.C)
    f.yz <- paste("ww.C", paste(y.lab, z.lab, sep = "+"), 
                  sep = "~")
    YZ.noCIA <- xtabs(as.formula(f.yz), data = data.C)
    if(is.factor(data.A[, y.lab])) rownames(YZ.noCIA) <- y.lev
    if(is.factor(data.B[, z.lab])) colnames(YZ.noCIA) <- z.lev   
          
#    dimnames(YZ.noCIA) <- list(y.lev, z.lev)
    out <- list(yz.CIA = YZ.CIA, cal.C = cal.C, yz.est = YZ.noCIA, 
                call = match.call())
  }
  if (micro) {

# prediction in A
    pred.Y.A <- X.A %*% beta.yx.A
    res.Y.A <- Y.A - pred.Y.A
    pred.Z.A <- X.A %*% beta.zx.B

#predictions in B
    pred.Z.B <- X.B %*% beta.zx.B
    res.Z.B <- Z.B - pred.Z.B
    pred.Y.B <- X.B %*% beta.yx.A

    if (!is.null(svy.C) & (estimation == "STWS" || estimation == 
                           "s2ws" || estimation == "synthetic")) {
      pred.Y.C <- X.C %*% beta.yx.A
      res.Y.C <- Y.C - pred.Y.C
      pred.Z.C <- X.C %*% beta.zx.B
      res.Z.C <- Z.C - pred.Z.C
      alfa2.1 <- t(res.Y.A) %*% (res.Y.A * w.A)
      alfa2.2 <- t(res.Y.C) %*% (res.Z.C * ww.C)
      qr.alfa <- qr(alfa2.1)
      alfa2 <- qr.coef(qr.alfa, alfa2.2)
      alfa2[is.na(alfa2)] <- 0
      cat(alfa2, fill = T)
      pred.Z.A <- pred.Z.A + res.Y.A %*% alfa2
      beta2.1 <- t(res.Z.B) %*% (res.Z.B * w.B)
      beta2.2 <- t(res.Z.C) %*% (res.Y.C * ww.C)
      qr.beta <- qr(beta2.1)
      beta2 <- qr.coef(qr.beta, beta2.2)
      beta2[is.na(beta2)] <- 0
      cat(beta2, fill = T)
      pred.Y.B <- pred.Y.B + res.Z.B %*% beta2
    }
    pred <- list(Y.A=pred.Y.A, Z.A = pred.Z.A, 
                 Z.B=pred.Z.B, Y.B = pred.Y.B)
    out <- c(out, pred)
  }
  out
}