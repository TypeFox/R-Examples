fortify.zoo <- function(model, data, melt = FALSE, ...)
{
  ## dimensions
  n <- NROW(model)
  k <- NCOL(model)

  ## series labels
  lab <- colnames(model)
  if(is.null(lab)) lab <- rep.int(deparse(substitute(model)), k)
  lab <- make.unique(lab)
  
  ## either long format (melt = TRUE) or wide format (melt = FALSE)
  if(melt) {
    df <- if(k == 1L) {
      data.frame(index(model), factor(rep.int(1, n), labels = lab), coredata(model))
    } else {
      data.frame(index(model)[rep.int(1:n, k)],
        factor(rep(1:k, each = n), levels = 1:k, labels = lab),
	as.vector(coredata(model)))
    }
    names(df) <- c("Index", "Series", "Value")
  } else {
    df <- cbind(data.frame(index(model)), coredata(model))
    names(df) <- c("Index", lab)  
  }
  
  return(df)
}


autoplot.zoo <- function(object, geom = "line", facets, ...)
{
  ## convert to data.frame (and assure correct label
  ## processing by fortify.zoo)
  lab <- deparse(substitute(object))
  if(NCOL(object) == 1L) {
    if(is.null(dim(object))) dim(object) <- c(NROW(object), 1L)
    if(is.null(colnames(object))) colnames(object) <- lab
  }
  if(is.null(colnames(object))) colnames(object) <- paste(lab, 1:NCOL(object), sep = ".")
  df <- fortify.zoo(object, melt = TRUE)

  ## default for facets
  single <- nlevels(df$Series) == 1L
  if(missing(facets)) {
    auto <- TRUE
    facets <- if(single) NULL else Series ~ .
  } else {
    auto <- FALSE
  }  

  ## "fake" variables for nonstandard evaluation
  Index <- Value <- Series <- NULL
  
  ## call qplot
  gg <- if(single | (!is.null(facets) & auto)) {
    ggplot2::qplot(Index, Value, data = df, geom = geom, facets = facets, ...) + ggplot2::ylab(if(single) levels(df$Series) else "") + ggplot2::xlab("Index")
  } else {
    ggplot2::qplot(Index, Value, data = df, group = Series, geom = geom, facets = facets, colour = Series, ...) + ggplot2::ylab("") + ggplot2::xlab("Index")
  }
  return(gg)
}

facet_free <- function (facets = Series ~ ., margins = FALSE, scales = "free_y", ...) {
  ggplot2::facet_grid(facets, margins = margins, scales = scales, ...)
}

yearmon_trans <- function(format = "%b %Y", n = 5) {
  breaks. <- function(x) as.yearmon(scales::pretty_breaks(n)(x))
  format. <- function(x) format(x, format = format)
  scales::trans_new("yearmon", transform = as.numeric, inverse = as.yearmon,
    breaks = breaks., format = format.)
}

scale_x_yearmon <- function(..., format = "%b %Y", n = 5) {
  ggplot2::scale_x_continuous(..., trans = yearmon_trans(format, n))
}
scale_y_yearmon <- function(..., format = "%b %Y", n = 5) {
  ggplot2::scale_y_continuous(..., trans = yearmon_trans(format, n))
}

yearqtr_trans <- function(format = "%Y-%q", n = 5) {
  breaks. <- function(x) as.yearqtr(scales::pretty_breaks(n)(x))
  format. <- function(x) zoo::format.yearqtr(x, format = format)
  scales::trans_new("yearqtr", transform = as.numeric, inverse = as.yearqtr,
    breaks = breaks., format = format.)
}

scale_x_yearqtr <- function(..., format = "%Y-%q", n = 5) {
  ggplot2::scale_x_continuous(..., trans = yearqtr_trans(format, n))
}
scale_y_yearqtr <- function(..., format = "%Y-%q", n = 5) {
  ggplot2::scale_y_continuous(..., trans = yearqtr_trans(format, n))
}

