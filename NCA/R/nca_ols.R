nca_ols <-
function (loop.data, ...) {
  x         <- loop.data$x
  y         <- loop.data$y
  ols       <- lm(y~x)
  intercept <- unname(coef(ols)["(Intercept)"])
  slope     <- unname(coef(ols)["x"])
  ceiling   <- p_ceiling(loop.data, slope, intercept)
  effect    <- ceiling / loop.data$scope
  ineffs    <- p_ineffs(loop.data, intercept, slope)
  above     <- p_above(loop.data, slope, intercept)
  
  return(list(line=ols,
              ceiling=ceiling, slope=slope, effect=effect,
              intercept=intercept, above=above, ineffs=ineffs,
              bottleneck=NULL))
}