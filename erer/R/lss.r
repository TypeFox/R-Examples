lss <- function(n = 5, pos = 1, decreasing = TRUE, 
  order.by = c("Size", "Type")) 
{
  order.by <- match.arg(order.by)
  napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos = pos)))
  names <- ls(pos = pos)
  if (length(names) == 0) {
    return("No object available.")
  } else {
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
      capture.output(print(object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    out$Name <- rownames(out)
    out <- out[, c(6, 1:5)]
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
    out <- head(out, n)
    rownames(out) <- 1:nrow(out)
    return(out)
  }
}