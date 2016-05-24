TestFunc <- function(mean.formula, var.formula, model.df) {
  
  my.glm <- glm(formula = mean.formula,
                data = model.df)
  
  x <- predict(my.glm, se.fit = TRUE)
  
  my.dglm <- dglm(formula = mean.formula,
                  dformula = var.formula,
                  data = model.df)
  
  y <- predict(my.dglm, se.fit = TRUE)
  
  z <- predict(my.dglm$dispersion.fit, se.fit = TRUE)
  
  return(list(glm = x, dglm.mean = y, dglm.var = z))
}

