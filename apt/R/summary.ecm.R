summary.ecm <- function(object, digits = 3, ...) {
  x <- object
  DepVar.x <- c(paste("diff.", x$name.x, ".t_0  |", sep = ""), 
    rep("|", times = nrow(coef(summary(x$ecm.x))) - 1)) 
  DepVar.y <- c(paste("diff.", x$name.y, ".t_0  -", sep = ""),  
    rep("-", times = nrow(coef(summary(x$ecm.y))) - 1))
              
  coef.x <- cbind(DepVar = DepVar.x, 
    IndVar = rownames(coef(summary(x$ecm.x))), 
    round(data.frame(coef(summary(x$ecm.x))), digits = digits))
  coef.y <- cbind(DepVar = DepVar.y, 
    IndVar = rownames(coef(summary(x$ecm.y))), 
    round(data.frame(coef(summary(x$ecm.y))), digits = digits))
  rownames(coef.x) <- 1:nrow(coef.x)
  rownames(coef.y) <- 1:nrow(coef.y)
  
  coeff  <- rbind(coef.x, coef.y)
  colnames(coeff) <- c("DepVar  ","IndVar  ", "estimate","error",
    "t.value", "p.value")
  coeff$signif <- ifelse(coeff$p.value <= 0.01, "***", 
    ifelse(coeff$p.value > 0.01 & coeff$p.value <= 0.05, "**",
    ifelse(coeff$p.value > 0.05 & coeff$p.value <= 0.10, "*", 
    ifelse(coeff$p.value > 0.10 & coeff$p.value <= 0.15, ".", " ")))) 
  return(coeff)
} 