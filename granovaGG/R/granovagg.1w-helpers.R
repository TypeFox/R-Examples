IsFSignificant <- function(model.summary) {
  f.critical  <- qf(
    p   = 0.95,
    df1 = model.summary$fstatistic["numdf"],
    df2 = model.summary$fstatistic["dendf"]
  )
  return(
    as.logical(
      model.summary$fstatistic["value"] > f.critical
    )
  )
}

## Plotting Functions
InitializeGgplot_1w <- function() {
  return(ggplot())
}

GrandMeanLine <- function(owp) {
  return(
    geom_hline(
      color      = brewer.pal(n = 8, name = "Set1")[3],
      alpha      = I(1/2),
      size       = I(0.25),
      yintercept = owp$stats$grand.mean
    )
  )
}

GrandMeanPoint <- function(owp) {
  return(
    geom_point(
      aes_string(
        x = "0",
        y = "mean(score)",
        color = 'factor(paste("Grand Mean"))'
      ),
      size = 2.5,
      data = owp$data
    )
  )
}

JitteredScoresByGroupContrast <- function(owp, jj) {
  only.jitter.in.x.direction <- position_jitter(
    height = 0,
    width = GetDegreeOfJitter_1w(owp, jj)
  )

  return(
    geom_point(
      aes_string(
        x = "contrast",
        y = "score"
      ),
      alpha    = 1,
      size     = 2,
      data     = owp$data,
      position = only.jitter.in.x.direction
    )
  )
}

GetDegreeOfJitter_1w <- function(owp, jj) {
  result <- owp$params$horizontal.percent

  if (!is.null(jj)) {
    result <- (jj / 200) * owp$params$contrast.range.distance
  }
  else if (IsSmallestContrastDifferenceSmallerThanOnePercentOfDataResolution(owp)) {
      result <- GetSmallestDistanceBetweenAdjacentContrasts(owp$summary$contrast)
  }

  return(result)
}

GetSmallestDistanceBetweenAdjacentContrasts <- function(contrasts) {
  ordered.contrasts             <- sort(contrasts)
  adjacent.contrast.differences <- 1:(length(contrasts) - 1)

  for (i in adjacent.contrast.differences) {
    contrast.difference              <- abs(ordered.contrasts[i + 1] - ordered.contrasts[i])
    adjacent.contrast.differences[i] <- contrast.difference
  }

  return(min(adjacent.contrast.differences))
}

IsSmallestContrastDifferenceSmallerThanOnePercentOfDataResolution <- function(owp) {
  return(
    abs(GetSmallestDistanceBetweenAdjacentContrasts(owp$summary$contrast)) < owp$params$horizontal.percent
  )
}






ScaleX_1w <- function(owp) {
  return(
    scale_x_continuous(
      breaks = (owp$params$aggregate.x.breaks),
      labels = signif(owp$params$aggregate.x.breaks, digits = 2),
      limits = owp$params$x.range,
      expand = c(0.00, 0)
    )
  )
}

ScaleY_1w <- function(owp) {
  return(
    scale_y_continuous(
      breaks = (owp$params$aggregate.y.breaks),
      labels = signif(owp$params$aggregate.y.breaks, digits = 2),
      limits = owp$params$y.range,
      expand = c(0.00, 0)
    )
  )
}

