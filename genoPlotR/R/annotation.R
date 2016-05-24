################################
# Annotation class and methods
################################
# annotation is a set of text and one or two positions for each
# if one, other must be NA
annotation <- function(x1, x2=NA, text, rot=0, col="black"){
  if (missing(x1) | missing(text)) stop("Args x1 and text must be provided")
  if (!is.numeric(x1)) stop("x1 must be numeric")
  if (!is.na(x2) && !is.numeric(x2)) stop("x2 must be numeric")
  if (!is.character(text)) stop("text must be character")
  as.annotation(data.frame(x1=x1, x2=x2, text=text, stringsAsFactors=FALSE),
                rot=rot, col=col)
}
as.annotation <- function(df, x2=NA, rot=0, col="black"){
  if (is.annotation(df)) return(df)
  if (!all(c("x1", "text") %in% names(df)))
    stop("Data frame should have at least a x1 and text column")
  # attributes x2, col and arg to all rows if not defined
  if (is.null(df$x2)) df$x2 <- x2
  if (is.null(df$color)) df$color <- col
  if (is.null(df$rot)) df$rot <- rot
  class(df) <- c("annotation", "data.frame")
  df
}
is.annotation <- function(annotation){
  inherits(annotation, "annotation")  
}
range.annotation <- function(x, ...){
  annotation <- x
  range(annotation$x1, annotation$x2, na.rm=TRUE)
}
trim.annotation <- function(x, xlim=NULL, ...){
  annotation <- x
  xlim <- as.numeric(xlim)
  if (!is.null(xlim)){
    if (!is.numeric(xlim)) stop("xlim must be numeric")
    if (length(xlim) != 2) stop("xlim must be length 2")
    # to be accepted, x1 > xlim1 and, if x2=NA, xlim1 also < xlim1 or,
    # x2 < xlim2
    annotation <- annotation[annotation$x1 >= xlim[1] &
                             ((is.na(annotation$x2) &
                               annotation$x1 <= xlim[2]) |
                              (!is.na(annotation$x2) &
                               annotation$x2 <= xlim[2])),]
  }
  annotation
}
