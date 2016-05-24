##' @title Print \code{spotDist} objects
##' @description \code{\link{print}} function for \code{\link[=spot.dist]{spotDist}} objects
##' @param x \code{spotDist} object to be printed.
##' @param ... further arguments passed to \code{\link{print}}.
##' @method print spotDist
##' @export
##' @author Mikko Vihtakari
##' @seealso \code{\link{spot.dist}}

print.spotDist <- function(x, ...) {
  if(is.null(x$sample.name)) {
    title <- paste0("\"", class(x)[1], "\"", " object, sample name is not defined")
  } else title <- paste0("\"", class(x)[1], "\"", " object for ", x$sample.name)
  
  title2 <- paste0("Main axis type: ", x$main.type)
  
  if(is.null(x$unit)) {
    title3 <- paste0("Scaling factor: ", x$scaling.factor, ", units not defined")
    } else title3 <- paste0("Scaling factor: ", x$scaling.factor, " pixels/", x$unit)
  
  title4 <- "Including following elements:"
  title5 <- "Sample spot output:"
  title6 <- "Growth line output:"
  type <- unlist(lapply(x, function(k) class(k)[1]))
  NCOL <- unlist(lapply(x, function(k) {GA <- ncol(k); if(is.null(GA)) {GA <- 0}; GA}))
  NROW <- unlist(lapply(x, function(k) {GA <- nrow(k); if(is.null(GA)) {GA <- 0}; GA}))
  LEN <- unlist(lapply(x, function(k) length(k)))
  DA <- data.frame(list.elements = names(x), type = type, ncol = NCOL, nrow = NROW, length = LEN)
  row.names(DA) <- 1:nrow(DA)
  
  cat(title, sep = "\n")
  cat(title2, sep = "\n")
  cat(title3, sep = "\n")
  cat(NULL, sep = "\n")
  cat(title4, sep = "\n")
  print(format(DA, justify = "left"), row.names = FALSE)
  cat(NULL, sep = "\n")
  cat(title5, sep = "\n")
  lapply(x$output, function(k) print(format(k, justify = "left"), row.names = FALSE))
  cat(NULL, sep = "\n")
  cat(title6, sep = "\n")
  print(format(x$gb.projections.dist, justify = "left"), row.names = FALSE)
  }

