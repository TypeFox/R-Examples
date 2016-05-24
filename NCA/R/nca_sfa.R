nca_sfa <-
function (loop.data, mpy, cutoff, bottleneck.x) {
  x           <- loop.data$x
  y           <- loop.data$y

  sfa         <- NULL
  options(warn=-1)
  try(sfa <- sfa(y~x, form="production"), silent=TRUE)
  options(warn=0)
  if (is.null(sfa)) {
    message()
    message("Ignoring SFA due to errors !")
    message()

    return(list(line=NULL, ceiling=NA, slope=NA, effect=NA,
              intercept=NA, above=NA, ineffs=list(x=NA, y=NA, abs=NA, rel=NA),
              bottleneck=NA))
  }

  intercept   <- unname(coef(sfa)["Intercept"])
  slope       <- unname(coef(sfa)["x"])
  ceiling     <- p_ceiling(loop.data, slope, intercept)
  effect      <- ceiling / loop.data$scope
  ineffs      <- p_ineffs(loop.data, intercept, slope)
  above       <- p_above(loop.data, slope, intercept)
  bottleneck  <- p_bottleneck(loop.data, mpy, slope, intercept, cutoff, bottleneck.x)

  sfa$coef <- sfa$coef[c(-3, -4)]

  return(list(line=sfa,
              ceiling=ceiling, slope=slope, effect=effect,
              intercept=intercept, above=above, ineffs=ineffs,
              bottleneck=bottleneck))
}