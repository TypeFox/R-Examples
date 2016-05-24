setGeneric(
  name = "tsMake",
  def = function(object, ...)
    standardGeneric("tsMake")
)

setMethod(
  f = "tsMake",
  signature = "WqData",
  definition = function(object, focus, layer, type = c("ts.mon", "zoo"),
    qprob = NULL)
  {

    # Validate args
    d <- data.frame(object)
    if ( missing(focus) || length(focus) > 1 )
      stop("'focus' must be the name of a single site or variable.")
    if (match(focus, d$site, nomatch = 0) > 0) {
      d <- d[d$site == focus, ]
      if (nrow(d) == 0) 
        stop("No data for this site.")
    } else {
      if (match(focus, d$variable, nomatch = 0) > 0) {
        d <- d[d$variable == focus, ]
        if (nrow(d) == 0) 
          stop("No data for this variable.")
      } else {
        stop("'focus' does not match any sites or variables")
      }
    }
    type <- match.arg(type)

    # Assemble all depths
    depths <- NULL
    if (missing(layer))
      layer <- list(c(-Inf, Inf))
    if (identical(layer, 'max.depths')) {
      ans <- aggregate(depth ~ time + site + variable, data = d, max, na.rm = TRUE)
      d <- merge(d, ans)
      d$depth <- depths <- max(d$depth, na.rm=TRUE) + 1
    } else {
      if (!is(layer, "list"))
        layer <- list(layer)
      for (el in layer) {
        if ( !is(el, "numeric") || length(el) > 2 )
          stop("layer is not specified correctly")
        if (length(el) > 1) {
          depths1 <- unique(d[d$depth >= el[1] & d$depth <= el[2], "depth"])
          depths <- c(depths, depths1)
        } else {
          depths <- c(depths, el)
        }
      }
    }
    d <- d[d$depth %in% depths, ]
    if (nrow(d) == 0) 
      stop("No data for this layer.")

    # Define aggregation function
    if (is.null(qprob)) {
      f = mean
    } else {
      f = function(x, ...) quantile(x, probs = qprob, ...)
    }

    # Reshape data
    if (match(focus, d$site, nomatch = 0) > 0) {
      c1 <- dcast(d, time ~ variable, fun.aggregate = f, na.rm = TRUE)
    } else {
      c1 <- dcast(d, time ~ site, fun.aggregate = f, na.rm = TRUE)
    } 

    # Create zoo or ts object
    z1 <- zoo(c1[, -1], c1[, 1])
    if (type == 'ts.mon') {
      z1 <- aggregate(z1, as.yearmon, f, na.rm = TRUE)
      z1names <- colnames(z1)
      z1 <- as.ts(z1)
      colnames(z1) <- z1names
    }
    z1
  }
)
