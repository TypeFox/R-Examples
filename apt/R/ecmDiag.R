ecmDiag <- function(m, digits = 2)
{
  if (!inherits(m, "ecm")) {
    stop("\n Provide an 'ecm' objectt from ecmSymFit or ecmAsyFit.\n")}

  dia.name <- c("R-squared", "Adj-R2", "F-stat", "Stat DW", "p-value DW", 
      "AIC", "BIC", "LB(4)", "LB(8)", "LB(12)") 
  dia <- data.frame(matrix(NA, nrow = length(dia.name), ncol = 2))
  colnames(dia) <- c(m$name.x, m$name.y)
  
  dia[1, ] <- c(summary(m$ecm.x)$r.squared,     summary(m$ecm.y)$r.squared)
  dia[2, ] <- c(summary(m$ecm.x)$adj.r.squared, summary(m$ecm.y)$adj.r.squared)
  dia[3, ] <- c(summary(m$ecm.x)$fstatistic[1], summary(m$ecm.y)$fstatistic[1])
  dia[4, ] <- c(durbinWatsonTest(m$ecm.x)$dw,   durbinWatsonTest(m$ecm.y)$dw)   
  dia[5, ] <- c(durbinWatsonTest(m$ecm.x)$p,    durbinWatsonTest(m$ecm.y)$p)
      
  dia[6, ] <- c(AIC(m$ecm.x, k = 2), AIC(m$ecm.y, k = 2))
  dia[7, ] <- c(AIC(m$ecm.x, k = log(nrow(m$IndVar))), 
                AIC(m$ecm.y, k = log(nrow(m$IndVar))))
  dia[8, ] <- c(Box.test(residuals(m$ecm.x), lag=4,  type="Ljung")$p.value,
                Box.test(residuals(m$ecm.y), lag=4,  type="Ljung")$p.value)
  dia[9, ] <- c(Box.test(residuals(m$ecm.x), lag=8,  type="Ljung")$p.value,
                Box.test(residuals(m$ecm.y), lag=8,  type="Ljung")$p.value)
  dia[10, ]<- c(Box.test(residuals(m$ecm.x), lag=12, type="Ljung")$p.value,
                Box.test(residuals(m$ecm.y), lag=12, type="Ljung")$p.value) 
  dia <- round(dia, digits = digits)
  dia$item <- dia.name   
  dia <- dia[, c(3, 1, 2)]          
  return(dia)
} 