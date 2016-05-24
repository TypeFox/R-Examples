#' Calculate distance scores on data in preparation for composite scoring
#'
#' @param d The data
#' @param g A grouping variable
#' @param thresholds Thresholds to use when calculating distances
#'   (e.g., median, clinical thresholds, etc.).  If groups are used,
#'   must be thresholds for each group (e.g., to allow separate thresholds for
#'   females and males).
#' @param higherisbetter A logical vector for each biomarker whether higher scores
#'   are better or not.
#' @param winsorize Whether to winsorize the data or not.  Defaults to \code{FALSE}.
#'   If not \code{FALSE}, the percentile to winsorize at.  For example, .01 would be
#'   the .01 and the 1 - .01 percentiles.
#' @param better Logical indicating whether \dQuote{better} values than the threshold
#'   are allowed. Defaults to \code{TRUE}.
#' @param na.rm A logical whether missing values should be ommitted. Defaults to
#'   \code{TRUE}.
#' @param saveall A logical whether to save all intermediary datasets and graphs.
#'   Defaults to \code{FALSE}.
#' @return A list of results.
#' @export
#' @family composite
#' @examples
#' # this example creates distances for the built in mtcars data
#' # see ?mtcars for more details
#' # The distances are calculated from the "best" in the dataset
#' # defined by these thresholds
#' thresholds <- with(mtcars, c(
#'   mpg = max(mpg),
#'   hp = max(hp),
#'   wt = min(wt),
#'   qsec = min(qsec)))
#'
#' # higher mpg and hp are better,
#' # whereas lower wt and qsec are better
#' dres <- distanceScores(mtcars[, c("mpg", "hp", "wt", "qsec")],
#'   thresholds = list(thresholds),
#'   higherisbetter = c(TRUE, TRUE, FALSE, FALSE),
#'   saveall = TRUE)
#'
#' # see a density plot of the distance scores
#' dres$Density
#'
#' # cleanup
#' rm(thresholds, dres)
distanceScores <- function(d, g, thresholds, higherisbetter, winsorize = FALSE, better = TRUE, na.rm = TRUE, saveall = FALSE) {
  fcall <- match.call()

  # data and input checks
  stopifnot(all(apply(d, 2, is.numeric)))
  k <- ncol(d)

  if (missing(higherisbetter)) {
      higherisbetter <- rep(0, k)
  }

  if (missing(thresholds)) {
      thresholds <- list(rep(0, k))
  }

  if (missing(g)) {
      g <- rep("1", nrow(d))
      stopifnot(identical(length(thresholds), 1L))
      if (is.null(names(thresholds))) {
        names(thresholds) <- "1"
      }
  }
  ng <- length(unique(g))

  stopifnot(identical(length(higherisbetter), k))
  stopifnot(all(sapply(thresholds, length) == k))
  stopifnot(identical(length(g), nrow(d)))

  # handle missing values
  if (na.rm) {
      okindex <- which(rowSums(is.na(cbind(d, g))) == 0)
      d <- d[okindex, ]
      g <- g[okindex]
      if (!identical(length(unique(g)), ng)) {
          warning("After removing missing values, levels of grouping variable no longer equal")
      }
  }

  if (any(unlist(thresholds) != 0) & saveall) d.original <- d


  # create the threshold matrix
  thresholdmatrix <- t(sapply(g, function(i) thresholds[[as.character(i)]]))
  rownames(thresholdmatrix) <- NULL
  colnames(thresholdmatrix) <- colnames(d)

  if (winsorize) {
      d <- winsorizor(d, percentile = winsorize, na.rm = TRUE)
      winsorized <- attr(d, "winsorized")
      if (saveall) d.winsorize <- d
  } else {
      winsorized <- NULL
  }


  d <- as.data.frame(sapply(1:k, function(i) {
    if (higherisbetter[i] == 0) {
      d[, i] - thresholdmatrix[, i]
    } else if (higherisbetter[i] == 1) {
      thresholdmatrix[, i] - d[, i]
    }
  }))

  colnames(d) <- colnames(thresholdmatrix)


  if (saveall) d.distances <- d

  # if "better" than threshold are not allowed
  # than truncate at zero, otherwise leave as is
  if (!better) {
      d <- as.data.frame(apply(d, 2, pmax, 0))
  }

  if (saveall) {
    data <- list(
      Distance = if(!better) d.distances else NULL,
      Winsorized = if(winsorize) d.winsorize else NULL,
      Raw = if(exists("d.original")) d.original else NULL)

    plots <- list(
      Distance = if(!better) ldensity(cbind(d.distances, Group = g), melt = TRUE, g = "Group") else NULL,
      Winsorized = if(winsorize) ldensity(cbind(d.winsorize, Group = g), melt = TRUE, g = "Group") else NULL,
      Raw = if(exists("d.original")) ldensity(cbind(d.original, Group = g), melt = TRUE, g = "Group") else NULL)
  } else {
    plots <- data <- list(Distance = NULL, Windsorized = NULL, Raw = NULL)
  }

  out <- list(
      Distances = d,
      Density = ldensity(cbind(d, Group = g), melt = TRUE, g = "Group"),
      Groups = g,
      SavedPlots = plots, SavedData = data,
      thresholds = thresholds, higherisbetter = higherisbetter,
      winsorize = winsorize, winsorized = winsorized,
      better = better, na.rm = na.rm,
      call = fcall)
  class(out) <- c("distancescores", "list")

  return(out)
}



#' Prepare data to have a composite calculated
#'
#' @param object An object ready for use
#' @param covmat The covariance matrix to use.  If missing,
#'   austomatically calculated from the data.
#' @param standardize A logical value whether to standardize the data or not.
#'   Defaults to \code{TRUE}.
#' @return A list of results.
#' @export
#' @family composite
#' @examples
#' # this example creates distances for the built in mtcars data
#' # see ?mtcars for more details
#' # The distances are calculated from the "best" in the dataset
#' # defined by these thresholds
#' thresholds <- with(mtcars, c(
#'   mpg = max(mpg),
#'   hp = max(hp),
#'   wt = min(wt),
#'   qsec = min(qsec)))
#'
#' # higher mpg and hp are better,
#' # whereas lower wt and qsec are better
#' dres <- distanceScores(mtcars[, c("mpg", "hp", "wt", "qsec")],
#'   thresholds = list(thresholds),
#'   higherisbetter = c(TRUE, TRUE, FALSE, FALSE),
#'   saveall = TRUE)
#'
#' # see a density plot of the distance scores
#' dres$Density
#'
#' # now prepare to create the composite
#' # covariance matrix will be calculated from the data
#' # and data will be standardized to unit variance by default
#' cprep <- prepareComposite(dres)
#' # cleanup
#' rm(thresholds, dres, cprep)
prepareComposite <- function(object, covmat, standardize = TRUE) {
    if (!inherits(object, "distancescores")) {
        warning(paste("Object is not of type 'distancescores'.",
                      "prepareComposite() may not work correctly."))
    }

    k <- ncol(object$Distances)
    if (missing(covmat)) {
        covmat <- cov(object$Distances)
    }

    stopifnot(identical(k, ncol(covmat)))

    sigma <- sqrt(diag(covmat))

    if (standardize) {
        data <- sweep(object$Distances, 2, sigma, "/")
    } else {
        data <- object$Distances
    }

    out <- c(list(
        data = data,
        covmat = covmat,
        sigma = sigma,
        standardize = standardize,
        k = k), object)

    class(out) <- c("compositedata", "list")

    return(out)
}


#' Score Data Using the Mahalanobis Distance
#'
#' Create a composite using the Mahalanobis Distance
#'
#' @param object An object of class \code{compositedata} ready for use
#' @param ncomponents the number of components to use from the
#'   principal component analysis. If missing, defaults to the
#'   number of columns in the data.
#' @return A list of results.
#' @export
#' @family composite
#' @examples
#' # this example creates distances for the built in mtcars data
#' # see ?mtcars for more details
#' # The distances are calculated from the "best" in the dataset
#' # defined by these thresholds
#' thresholds <- with(mtcars, c(
#'   mpg = max(mpg),
#'   hp = max(hp),
#'   wt = min(wt),
#'   qsec = min(qsec)))
#'
#' # higher mpg and hp are better,
#' # whereas lower wt and qsec are better
#' dres <- distanceScores(mtcars[, c("mpg", "hp", "wt", "qsec")],
#'   thresholds = list(thresholds),
#'   higherisbetter = c(TRUE, TRUE, FALSE, FALSE),
#'   saveall = TRUE)
#'
#' # see a density plot of the distance scores
#' dres$Density
#'
#' # now prepare to create the composite
#' # covariance matrix will be calculated from the data
#' # and data will be standardized to unit variance by default
#' cprep <- prepareComposite(dres)
#'
#' # now we can create the composite based on mahalanobis distances
#' # from our defined thresholds
#' mcomp <- mahalanobisComposite(cprep)
#'
#' # view a histogram of the composite scores
#' mcomp$ScoreHistogram
#'
#' # summarize the composite scores
#' summary(mcomp$Scores)
#'
#' # check the screeplot and loadings
#' mcomp$Screeplot
#' mcomp$LoadingGraph
#' # examine the loadings as a table
#' mcomp$LoadingTable
#'
#' # one component is adequate to explain these data
#' # to be safe can pick first two and re-run model
#'
#' # use only first two components
#' mcomp2 <- mahalanobisComposite(cprep, ncomponents = 2)
#'
#' # view a histogram of the updated composite scores
#' mcomp2$ScoreHistogram
#'
#' # summarize the composite scores
#' summary(mcomp2$Scores)
#'
#' # compare using all versus two components
#' plot(mcomp$Scores, mcomp2$Scores)
#'
#' # cleanup
#' rm(thresholds, dres, cprep, mcomp, mcomp2)
mahalanobisComposite <- function(object, ncomponents) {
  if (!inherits(object, "compositedata")) {
      warning(paste("Object is not of type 'compositedata'.",
                    "mahalanobisComposite() may not work correctly."))
  }

  if (missing(ncomponents)) {
    ncomponents <- object$k
  }

  stopifnot(ncomponents <= object$k)

  pca <- princomp(scores = FALSE, cor = object$standardize, covmat = object$covmat)

  screeplot <- ggplot(data.frame(Component = 1:object$k, Eigenvalue = pca$sdev^2),
    aes(Component, Eigenvalue)) +
    geom_line() + geom_point() +
    geom_hline(aes(yintercept = 1)) +
    scale_x_continuous(breaks = 1:object$k) +
    theme_classic()

  Lbase <- as.matrix(pca$loadings[])

  ltab <- Lbase
  colnames(ltab) <- paste0("C", 1:object$k)
  ltab[] <- format(round(ltab, 2), nsmall=2, digits=2)

  L <- as.data.frame(pca$loadings[, 1:ncomponents])
  colnames(L) <- paste0("Comp", 1:ncomponents)
  L$Variable <- factor(rownames(L), levels = colnames(object$data)[1:ncomponents])
  L <- melt(L, id.vars = "Variable")

  loadingsplot <- ggplot(L, aes(Variable, abs(value), fill = ifelse(value > 0, "positive", "negative"))) +
    geom_bar(stat = "identity") +
    scale_fill_manual("Direction", values = c("positive" = "black", "negative" = "grey80")) +
    facet_grid(variable ~ .) +
    ylab("Absolute Loading") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust=.5))

  c.scores <- as.matrix(object$data) %*% Lbase %*% diag(1/pca$sdev)
  c.scores <- as.data.frame(c.scores)
  colnames(c.scores) <- paste0("C", 1:object$k)

  finalScores <- sqrt(rowSums(c.scores[, 1:ncomponents]^2))
  # alternate way
  # finalScores <- sqrt(rowSums((data %*% solve(cov2cor(covmat))) * data))

  out <- c(list(
      Scores = finalScores,
      ScoreHistogram = ldensity(data.frame(Scores = finalScores), x = "Scores", hist = TRUE),
      Screeplot = screeplot,
      LoadingGraph = loadingsplot,
      LoadingTable = ltab), object)
  class(out) <- c("mahalanobiscomposite", "list")

  return(out)
}

# clear R CMD CHECK notes
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Component", "Eigenvalue", "Variable", "value"))

#' Score Data Using a simple sum
#'
#' Create a composite using summation
#'
#' @param object An object ready for use
#' @param transform A character string indicating the type of transformation to use.
#'   One of \dQuote{square}, \dQuote{abs}, or \dQuote{none}, which either sums the raw data,
#'   sums the squared data and then takes the square root, or sums the absolute values of the
#'   data.
#' @param type A character string indicating the type of aggregation to use.
#'   One of \dQuote{sum} or \dQuote{mean}.
#' @param systems An optional list where each element is a character vector of the
#'   variable names within a particular system.  If given, scores are first averaged
#'   within a system, before being aggregated across systems.
#' @return A list of results.
#' @export
#' @family composite
#' @examples
#' # this example creates distances for the built in mtcars data
#' # see ?mtcars for more details
#' # The distances are calculated from the "best" in the dataset
#' # defined by these thresholds
#' thresholds <- with(mtcars, c(
#'   mpg = max(mpg),
#'   hp = max(hp),
#'   wt = min(wt),
#'   qsec = min(qsec)))
#'
#' # higher mpg and hp are better,
#' # whereas lower wt and qsec are better
#' dres <- distanceScores(mtcars[, c("mpg", "hp", "wt", "qsec")],
#'   thresholds = list(thresholds),
#'   higherisbetter = c(TRUE, TRUE, FALSE, FALSE),
#'   saveall = TRUE)
#'
#' # see a density plot of the distance scores
#' dres$Density
#'
#' # now prepare to create the composite
#' # covariance matrix will be calculated from the data
#' # and data will be standardized to unit variance by default
#' cprep <- prepareComposite(dres)
#'
#' # now we can create the composite based on summing the (standardized)
#' # distances from our defined thresholds
#' # by default, distances are squared, then summed, and then square rooted
#' # to be back on the original scale
#' scomp <- sumComposite(cprep, "square", "sum")
#'
#' # view a histogram and summary of the composite scores
#' scomp$ScoreHistogram
#' summary(scomp$Scores)
#'
#' # calculate average (mean) instead of sum
#' scomp2 <- sumComposite(cprep, "square", "mean")
#'
#' # view a histogram and summary of the composite scores
#' scomp2$ScoreHistogram
#' summary(scomp2$Scores)
#'
#' # scores are still the same
#' plot(scomp$Scores, scomp2$Scores)
#'
#' # first average scores within a system, then sum
#' # note that within a system, scores are always averaged,
#' # never summed.
#' scomp3 <- sumComposite(cprep, "square", "sum",
#'   systems = list(
#'     environment = c("mpg"),
#'     performance = c("hp", "qsec", "wt")))
#'
#' # view a histogram and summary of the composite scores
#' scomp3$ScoreHistogram
#' summary(scomp3$Scores)
#'
#' # compare all three scores
#' # because of the different number of indicators within each system
#' # there is a re-weighting for S3
#' plot(data.frame(S1 = scomp$Scores, S2 = scomp2$Scores, S3 = scomp3$Scores))
#'
#' # cleanup
#' rm(thresholds, dres, cprep, scomp, scomp2, scomp3)
sumComposite <- function(object, transform = c("square", "abs", "none"), type = c("sum", "mean"),
                         systems) {
  if (!inherits(object, "compositedata")) {
      warning(paste("Object is not of type 'compositedata'.",
                    "sumComposite() may not work correctly."))
  }

  transform <- match.arg(transform)
  aggregator <- switch(type,
                       sum = rowSums,
                       mean = rowMeans)

  if (!missing(systems)) {
      x <- sapply(systems, function(v) {
        switch(transform,
          square = rowMeans(object$data[, v, drop = FALSE]^2),
          abs = rowMeans(abs(object$data[, v, drop = FALSE])),
          none = rowMeans(object$data[, v, drop = FALSE]))
       })
      finalScores <- switch(transform,
         square = sqrt(aggregator(x)),
         abs = aggregator(x),
         none = aggregator(x))
  } else {
      finalScores <- switch(transform,
        square = sqrt(aggregator(object$data^2)),
        abs = aggregator(abs(object$data)),
        none = aggregator(object$data))
      systems <- NULL
  }

  out <- c(list(
      Scores = finalScores,
      ScoreHistogram = ldensity(data.frame(Scores = finalScores), x = "Scores", hist = TRUE),
      transform = transform,
      type = type,
      systems = systems), object)
  class(out) <- c("sumcomposite", "list")

  return(out)
}


#' Score Data Using a Factor Model
#'
#' Create a composite using a Factor Model
#'
#' @param object An object ready for use
#' @param type A character string indicating the type of factor model to use
#' @param factors A named list where names are the factor names and each
#'   element is a character string of the indicator names.
#' @return A list of results.
#' @import lavaan
#' @export
#' @family composite
#' @examples
#' # this example creates distances for the built in mtcars data
#' # see ?mtcars for more details
#' # The distances are calculated from the "best" in the dataset
#' # defined by these thresholds
#' thresholds <- with(mtcars, c(
#'   mpg = max(mpg),
#'   hp = max(hp),
#'   wt = min(wt),
#'   qsec = min(qsec)))
#'
#' # higher mpg and hp are better,
#' # whereas lower wt and qsec are better
#' dres <- distanceScores(mtcars[, c("mpg", "hp", "wt", "qsec")],
#'   thresholds = list(thresholds),
#'   higherisbetter = c(TRUE, TRUE, FALSE, FALSE),
#'   saveall = TRUE)
#'
#' # see a density plot of the distance scores
#' dres$Density
#'
#' # now prepare to create the composite
#' # covariance matrix will be calculated from the data
#' # and data will be standardized to unit variance by default
#' cprep <- prepareComposite(dres)
#'
#' # now we can create the composite based on summing the (standardized)
#' # distances from our defined thresholds
#' # by default, distances are squared, then summed, and then square rooted
#' # to be back on the original scale
#' fcomp <- factorComposite(cprep, type = "onefactor")
#'
#' # view a histogram of the composite scores
#' fcomp$ScoreHistogram
#'
#' # summarize the composite scores
#' summary(fcomp$Scores$Comp)
#'
#' # cleanup
#' rm(thresholds, dres, cprep, fcomp)
factorComposite <- function(object, type = c("onefactor", "secondorderfactor", "bifactor"), factors = "") {
    type <- match.arg(type)
    if (!inherits(object, "compositedata")) {
        warning(paste("Object is not of type 'compositedata'.",
                      "factorComposite() may not work correctly."))
    }

    vars <- colnames(object$data)
    unused <- setdiff(vars, unlist(factors))

    onefactor <- paste0("Comp =~ ", paste(vars, collapse = " + "))

    m.factors <- paste(sapply(1:length(factors), function(i) {
        sprintf("%s =~ %s", names(factors)[i],
                ifelse(length(factors[[i]]) < 3,
                       paste(paste0("b", i, "1", " * ", factors[[i]]), collapse = " + "),
                       paste(paste0("b", i, 1:length(factors[[i]]), " * ", factors[[i]]), collapse = " + ")))
    }), collapse = "\n")

    if (length(factors) > 1) {
        m.factors.covariances <- sapply(1:(length(factors) - 1), function(i) {
            sapply((i + 1):(length(factors)), function(j) {
                sprintf("%s ~~ 0 * %s", names(factors)[i], names(factors)[j])
            })
        })
    } else {
        m.factors.covariances <- ""
    }

    m.overall <- paste0("Comp =~ ", paste(c(names(factors), unused), collapse = " + "))

    secondorderfactor <- paste(m.factors, m.overall, sep = "\n")

    bifactor <- paste(
        m.factors,
        onefactor,
        paste(paste("Comp ~~ 0 * ", names(factors)), collapse = "\n"),
        paste(unlist(m.factors.covariances), collapse = "\n"),
        sep = "\n")

    fit <- switch(type,
                  onefactor = sem(onefactor, data = object$data),
                  secondorderfactor = sem(secondorderfactor, data = object$data),
                  bifactor = sem(bifactor, data = object$data))

    out <- c(list(
      Scores = as.data.frame(predict(fit)),
      ScoreHistogram = ldensity(as.data.frame(predict(fit)), x = "Comp", hist = TRUE),
      Fit = fit), object)

  class(out) <- c("factorcomposite", "list")

  return(out)
}
