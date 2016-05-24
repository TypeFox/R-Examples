nca_cols <-
function (loop.data, mpy, cutoff, bottleneck.x) {
  x         <- loop.data$x
  y         <- loop.data$y  
  cols      <- lm(y + max (residuals(lm(y~x)))~x)
  
  intercept   <- unname(coef(cols)["(Intercept)"])
  slope       <- unname( coef(cols)["x"])
  ceiling     <- p_ceiling(loop.data, slope, intercept)
  effect      <- ceiling / loop.data$scope
  ineffs      <- p_ineffs(loop.data, intercept, slope)
  bottleneck  <- p_bottleneck(loop.data, mpy, slope, intercept, cutoff, bottleneck.x)
  above       <- if(slope > 0) 0 else NA

  return(list(line=cols,
              ceiling=ceiling, slope=slope, effect=effect,
              intercept=intercept, above=above, ineffs=ineffs,
              bottleneck=bottleneck))
}