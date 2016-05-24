linearRegressionTable <- function(model.fit) {
  
  tab <- data.frame(summary(model.fit)$coefficients)
  colnames(tab) <- colnames(summary(model.fit)$coefficients)
  tab[ ,1:3] <- round(tab[ ,1:3], 2)
  tab[ ,4] <- ifelse(tab[ ,4] < .001, "< 0.001", round(tab[ ,4], 3))
  
  design.matrix <- model.fit$model
  numeric.columns <- design.matrix[,unlist(lapply(design.matrix,is.numeric))]
  scaled.numeric.columns <- scale(numeric.columns)
  design.matrix[,unlist(lapply(design.matrix,is.numeric))] <- scaled.numeric.columns 
  model.fit.b <- stats::update(model.fit, data = design.matrix)
  full.table <- data.frame(B = tab[,1], SE_B = tab[,2], BETA = round(model.fit.b$coefficients, 2), t = tab[,3], p = tab[,4])
  full.table
}
