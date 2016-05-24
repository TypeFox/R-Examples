
#---------------
# Configure functions
#---------------

#' Get a list of available transforms.
#'
#' @return A character vector of transform names.
#' @examples
#'  availableTransforms()
#' @export
availableTransforms <- function() c("none", "inormal", "log", "out")

#---------------
# Transform functions
#---------------

#' Apply transforms to a data set.
#'
#' The function looks for traits in colulms of input data frame,
#' call \code{\link{transformTrait}} function per trait,
#' and rename the traits.
#' By default, the name of a transformed trait is updated 
#' to one with a prefixed given in \code{transform.prefix} argument
#' (the default value is \code{"tr_"}).
#' Such renaming is assumed to make the user aware that the trait is transformed.
#'
#' This function is internally called in \code{solarPolygenic}
#' if \code{transforms} argument is specified.
#' In this case of the polygenic analysis,
#' the transform operation is invisible to the user.
#' However, it is recommended to manually transform traits
#' in other linkage and association analyses.
#'
#' @param transforms 
#'    A named character vector of transforms,
#'    where traits are given in the names of the vector.
#' @param data
#'    A matrix or a data.frame, where the column names represent the names of traits.
#' @param transform.prefix
#'    A character vector, that is a prefix to be added to the name of transformed trait.
#'    The default value is \code{"tr_"}.
#' @param ...
#'    Additional parameters to be passed to \code{transformTrait} function called inside of \code{transformData}.
#'    For example, it might be a parameter \code{log.base} for \code{\link{transformTrait}} function
#'    in the case \code{transform} is equal to \code{"log"}.
#' @return
#'    A matrix or a data.frame of the transformed data.
#'
#' @seealso \code{\link{availableTransforms}}, \code{\link{transformTrait}}
#' @export
transformData <- function(transforms, data, transform.prefix = "tr_", ...)
{
  ### arg
  stopifnot(!missing(transforms))
  stopifnot(!missing(data))
  
  stopifnot(class(transforms) == "character")
  stopifnot(all(transforms %in% availableTransforms()))
  stopifnot(!is.null(names(transforms)))
  
  traits <- names(transforms)
  stopifnot(length(traits) == length(transforms))
  stopifnot(all(traits %in% colnames(data)))
  
  for(i in 1:length(transforms)) {
    var <- traits[i]
    tr <- transforms[i]
    
    x <- data[, var]    
    xt <- transformTrait(x, tr, ...)
    
    vart <- paste0(transform.prefix, var)
    data[, vart] <- xt
  }
  
  return(data)
}

#' Transform a trait.
#'
#' @param x 
#'    a numeric vector (of a trait).
#' @param transform
#'    a character vector, the name of transformation.
#'    Possible values are returned by \code{\link{availableTransforms}} function.
#' @param mult
#'    A numeric, the multiplicator for the transformed value of a trait.
#'    The default value is \code{1}.
#' @param ...
#'    additional parameters passed to internal \code{transform_trait_*} functions.
#'    Possible parameters might be \code{log.base}, \code{log.intercept} (\code{"log"} transformation).
#'
#' @return A numeric vector, which contains the transformed values (of a trait).
#'
#' @examples
#' library(plyr)
#' library(ggplot2)
#'
#' data(dat30)
#' dat <- mutate(dat30,
#'    inormal_trait1 = transformTrait(trait1, "inormal"))
#' 
#' ggplot(dat, aes(trait1)) + geom_histogram()
#' ggplot(dat, aes(inormal_trait1)) + geom_histogram()
#'
#' @seealso \code{\link{availableTransforms}}, \code{\link{transformData}}
#'
#' @export
transformTrait <- function(x, transform, mult = 1, ...)
{
  ### arg
  stopifnot(!missing(x))
  stopifnot(!missing(transform))  
  
  stopifnot(class(transform) == "character")
  stopifnot(length(transform) == 1)
  stopifnot(transform %in% availableTransforms())
  
  stopifnot(class(x) %in% c("integer", "numeric"))
  
  res <- switch(transform,
    "none" = list(x = x),
    "log" = list(x = transform_trait_log(x, ...)),
    "inormal" = list(x = transform_trait_inormal(x, ...)),
    "out" = list(x = transform_trait_out(x, ...)),
    stop("error in switch (unknown transform)"))

  xt <- res$x
  
  # mult
  xt <- mult * xt
  
  return(xt)
}

transform_trait_log <- function(x, log.base, log.intercept, log.xintercept)
{
  stopifnot(!missing(x))
  
  ### case 1
  if(missing(log.intercept) & missing(log.xintercept)) {
    x.min <- min(x, na.rm = TRUE)
    if(x.min <= 0) {
      x.min <- min(x, na.rm = TRUE)
      x <- x - x.min + 0.1
    }
  } 
  
  ### case 2
  if(!missing(log.intercept)) {
    x <- x - log.intercept

    x.min <- min(x, na.rm = TRUE)
    stopifnot(x.min > 0)
  }
  
  ### case 3
  if(!missing(log.xintercept)) {
    x.min <- min(x, na.rm = TRUE)
    stopifnot(log.intercept > x.min)
    
    # y = kx + b
    k <- 0.9 / (log.intercept - x.min)
    b <- 1 - 0.9 * log.intercept / (log.intercept - x.min)
    
    x <- k * x + b

    x.min <- min(x, na.rm = TRUE)
    stopifnot(x.min > 0)
  }
  
  if(missing(log.base)) {
    xt <- log(x)
  } else {
    xt <- log(x, log.base)
  }

  return(xt)
}  
  
# test data: x <- c(NA, 10:1, NA)
transform_trait_inormal <- function(x, mean = 0, sd = 1)
{
  stopifnot(!missing(x))
  
  df0 <- data.frame(sample = 1:length(x), x = x) # data.frame with NA
  n0 <- nrow(df0)
  
  df <- subset(df0, !is.na(x)) # data.frame with NA removed

  df <- df[order(df$x), ]

  n <- nrow(df)
  df$y <- qnorm((1:n) / (n + 1), mean = mean, sd = sd)

  of <- join(df0, df[, c("sample", "y")], by = "sample") # join input and output data.frames
  of <- of[with(of, order(sample)), ]
  
  ### duplicated values
  y <- NULL # due to R CMD check: no visible binding
  of <- ddply(of, .(x), mutate, 
    y = median(y, na.rm = TRUE))

  of <- arrange(of, sample)
  
  return(of$y)
}

transform_trait_out <- function(x, threshold = 4)
{
    repeat {
        sd.pre <- sd(x,na.rm=TRUE)
        mean.pre <- mean(x,na.rm=TRUE)
        
        x[x > (mean.pre + threshold * sd.pre)] <- NA
        x[x < (mean.pre - threshold * sd.pre)] <- NA
        
        sd.post <- sd(x,na.rm=T)
        if (sd.pre == sd.post) break
    }
    return(x)
}
