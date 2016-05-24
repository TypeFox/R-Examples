summary.ciTarFit <- function(object, digits = 3, ...)
{
  w <- object
  mn <- ifelse(w$model == "tar", 
          ifelse(w$thresh == 0, w$model, "c.tar"),
              ifelse(w$thresh == 0, w$model, "c.mtar"))                 
 
  dia <- data.frame(matrix(NA, nrow = 14, ncol = 2))
  colnames(dia) <- c("item", mn)
  dia[1, ] <- c("lag",       w$lag)
  dia[2, ] <- c("thresh",    w$thresh)
  dia[3, ] <- c("total obs", NROW(w$y))
  dia[4, ] <- c("coint obs", nrow(w$data.CI))
  dia[5, ] <- c("sse", deviance(w$CI))
  dia[6, ] <- c("aic", AIC(w$CI, k = 2))
  dia[7, ] <- c("bic", AIC(w$CI, k=log(nrow(w$data.CI))))
  dia[8, ] <- c("H1: no CI",    w$f.phi[, 5][2])
  dia[9, ] <- c("H2: no APT",   w$f.apt[, 5][2] )
  dia[10, ] <- c("H1: p.value", w$f.phi[, 6][2])
  dia[11, ] <- c("H2: p.value", w$f.apt[, 6][2])
  dia[12, ] <- c("LB test(4)", Box.test(resid(w$CI), lag = 4, 
    type = "Ljung")$p.value)
  dia[13, ] <- c("LB test(8)", Box.test(resid(w$CI), lag = 8, 
    type = "Ljung")$p.value)
  dia[14, ] <- c("LB test(12)",Box.test(resid(w$CI), lag = 12,
    type = "Ljung")$p.value)
  dia[, 2] <- round(as.numeric(dia[, 2]), digits)
  
  c.LR <- round(data.frame(coef(summary(w$LR))), digits = digits)
  c.CI <- round(data.frame(coef(summary(w$CI))), digits = digits)    
  out  <- rbind(cbind(model=mn, reg="LR", c.LR), 
    cbind(model = mn, reg = "CI", c.CI))
  colnames(out) <- c("model", "reg", "estimate", "st.error", 
    "t.value", "p.value")
  out$sign <- ifelse(out$p.value <= 0.01, "***", 
    ifelse(out$p.value > 0.01 & out$p.value <= 0.05,   "**",
    ifelse(out$p.value > 0.05 & out$p.value <= 0.10, "*",
    ifelse(out$p.value > 0.10 & out$p.value <= 0.15, ".", " ")))) 
  out <- data.frame(variable = rownames(out), out)
  rownames(out) <- 1:nrow(out)
  out <- out[, c(2, 3, 1, 4:8)]
  result <- listn(dia, out)
  return(result)
} 