GenerateFactorData <- function(factor.levels.list, sample.size) {
  ## This function is intended to be used for generating fake data to
  ## test other functions.
  ## Args:
  ##   factor.levels.list: A named list.  Each list element is a
  ##     character vector giving the levels of the factors to be
  ##     sampled.  The list names will become the variable names.
  ##   sample.size:  The number of rows in the generated data frame.
  ##
  ## Returns:
  ##   A data frame  with sample.size rows and length(variable.names)
  ##   columns.  Each column is a factor.

  stopifnot(is.numeric(sample.size) &&
            length(sample.size) == 1 &&
            sample.size > 0)
  stopifnot(is.list(factor.levels.list))
  repeated.levels <- any(sapply(factor.levels.list, anyDuplicated))
  if (repeated.levels) {
    stop("All factor levels in factor.levels.list must be unique")
  }

  if (is.null(names(factor.levels.list))) {
    names(factor.levels.list) <- make.names(seq_along(factor.levels.list))
  }
  variable.names <- names(factor.levels.list)

  ans <- list()
  for (v in 1:length(factor.levels.list)) {
    level.names <- as.character(factor.levels.list[[v]])
    values <- factor(sample(level.names,
                            replace = TRUE,
                            size = sample.size),
                     levels = level.names)
    ans[[variable.names[v]]] <- values
  }
  return(as.data.frame(ans))
}
