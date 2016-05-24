
setClass(
  Class = 'WqData',
  contains = 'data.frame',
  validity = function(object) {
    if (!identical(object@names[1:5], c("time", "site", "depth",
    	"variable", "value")))
      stop("columns are not all named correctly")
    if (!all(
        is(object$time, "DateTime"),
        is(object$site, "factor"),
        is(object$depth, "numeric"),
        is(object$variable, "factor"),
        is(object$value, "numeric")))
      stop("columns are not all of correct class")
  }
)

setMethod(
  f = `[`,
  signature = "WqData",
  definition = function(x, i, j="MISSING", drop="MISSING") {
    if (missing(i))
      r <- TRUE
    else {
      if (is(i, "numeric"))
        r <- i
      else {
        e <- substitute(i)
        r <- eval(e, x, parent.frame())
        if (!is.logical(r))
          stop("'i' must be logical or numeric")
        r <- r & !is.na(r)
      }
    }
    x.stored <- x
    df <- `[.data.frame`(x, r, j=1:5, drop=FALSE)
    for (slot.name in names(getSlots("data.frame"))) {
      slot(x.stored, slot.name) <- slot(df, slot.name)
    }
    x.stored
  }
)

setMethod(
  f = "summary",
  signature = "WqData",
  definition = function(object, ...) {
    trange <- range(as.Date(format(object$time)), na.rm = TRUE)
    cat("date range: ", paste(trange[1], "to", trange[2]), "\n\n")
    nums <- table(object$site, object$variable)
    quarts <- tapply(object$value, object$variable, summary)
    quarts1 <- matrix(unlist(quarts), byrow = TRUE,  ncol = 6)
    colnames(quarts1) <- names(quarts[[1]])
    rownames(quarts1) <- names(quarts)[seq_len(nrow(quarts1))]
    sumry <- list(observations = nums, quartiles = quarts1)
    sumry
  }
)

setMethod(
  f = "plot",
  signature = "WqData",
  definition =  function(x, y = "missing", vars, num.col = NULL) {
    if (missing(vars))
      vars <- unique(x$variable)
      num.plots <- max(10, length(vars))
      vars <- vars[1:num.plots]
    d <- data.frame(x)
    d <- d[d$variable %in% vars, ]
    ggplot(d, aes_string(x = "site", y = "value", z = "variable")) +
      geom_boxplot(outlier.colour = 'blue', outlier.shape = 1) +
      facet_wrap(~ variable, scales = "free_y", ncol = num.col)
  }
)
