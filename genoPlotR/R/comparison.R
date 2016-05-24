################################
# Comparison class and methods
################################
# A comparison is mostly a bunch of similarities
comparison <- function(x){
  # if data.frame
  if (is.data.frame(x)){
    return(as.comparison(x))
  } else if (is.list(x)){
    return(as.comparison(as.data.frame(x, stringsAsFactors=FALSE)))
  } else {
    stop(paste("Cannot coerce class", class(x), "to comparison"))
  }
}
as.comparison <- function(df){
  # check for class comparison, list, df
  if (is.comparison(comparison)) return(df)
  if (is.list(df) && !is.data.frame(df)) df <- comparison(df)
  if (is.data.frame(df)) {
    # check for columns
    names <- c("start1", "end1", "start2", "end2")
    if (!identical(names(df)[1:4], names))
      stop("Col names should start with start1, end1, start2 and end2")
    if (!all(sapply(df, function(x) all(!is.null(x)))))
      stop("NULL values not allowed in data.frame")
    if (!is.numeric(df$start1) | !is.numeric(df$end1) |
        !is.numeric(df$start2) | !is.numeric(df$end2)){
      stop("Starts and ends must be numeric")
    }
    # add direction column
    df$direction <- ifelse(sign(df$start1-df$end1)
                           * sign(df$start2-df$end2) > 0, 1, -1)
    # check color (character and default)
    if (is.factor(df$col)) df$col <- as.character(df$col)
  }
  else {
    stop("Unable to handle this format")
  }
  class(df) <- c("comparison", "data.frame")
  df
}
is.comparison <- function(comparison){
  inherits(comparison, "comparison")
}
# gives the range of a comparison
range.comparison <- function(x, overall=TRUE, ...){
  comparison <- x
  if (overall){
    range <- range(comparison$start1, comparison$end1,
                   comparison$start2, comparison$end2, na.rm=FALSE)
  } else {
    xlim1 <- range(comparison$start1, comparison$end1, na.rm=FALSE)
    xlim2 <- range(comparison$start2, comparison$end2, na.rm=FALSE)
    range <- data.frame(xlim1=xlim1, xlim2=xlim2)
  }
  range
}
# trim comparison given x limits
trim.comparison <- function(x, xlim1=c(-Inf, Inf), xlim2=c(-Inf, Inf), ...){
  comparison <- x
  if (!is.null(xlim1) && !is.null(xlim2)){
    if (!is.numeric(xlim1) || !is.numeric(xlim2)) stop("xlims must be numeric")
    if (length(xlim1) != 2 || length(xlim2) != 2) stop("xlims must be length 2")
    ## testing to include overlapping comps
    ## direction 1
    comparison$start1[comparison$start1 < xlim1[1] &
                      comparison$end1 > xlim1[1]] <- xlim1[1]
    comparison$end1[comparison$start1 < xlim1[2] &
                    comparison$end1 > xlim1[2]] <- xlim1[2]
    comparison$start2[comparison$start2 < xlim2[1] &
                      comparison$end2 > xlim2[1]] <- xlim2[1]
    comparison$end2[comparison$start2 < xlim2[2] &
                    comparison$end2 > xlim2[2]] <- xlim2[2]
    ## direction -1
    comparison$start1[comparison$start1 > xlim1[2] &
                      comparison$end1 < xlim1[2]] <- xlim1[2]
    comparison$end1[comparison$start1 > xlim1[1] &
                    comparison$end1 < xlim1[1]] <- xlim1[1]
    comparison$start2[comparison$start2 > xlim2[2] &
                      comparison$end2 < xlim2[2]] <- xlim2[2]
    comparison$end2[comparison$start2 > xlim2[1] &
                    comparison$end2 < xlim2[1]] <- xlim2[1]
    comparison <-
      comparison[comparison$start1 >= xlim1[1] & comparison$end1 <= xlim1[2] &
                 comparison$start2 >= xlim2[1] & comparison$end2 <= xlim2[2],]
  }
  comparison
}
# reverses a comparison. side <1 for no sides, 1 for first,
# 2 for second, >2 for bth
reverse.comparison <- function(x, side=0, ...){
  comparison <- x
  if (side > 1){
    comparison$start2 <- -comparison$start2
    comparison$end2 <- -comparison$end2
  }
  if (side == 1 || side >2){
    comparison$start1 <- -comparison$start1
    comparison$end1 <- -comparison$end1    
  }
  # if only one side changed, change direction
  if (side == 1 || side == 2){
    comparison$direction <- -comparison$direction
  }
  comparison
}
