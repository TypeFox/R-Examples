#' @encoding UTF-8
#' @title Stratified Sampling
#'
#' @description A handy function for sampling row values of a data.frame conditional to some strata.
#'
#' @param .data The data.frame from which the sample is desired.
#' @param group The grouping factor, may be a list.
#' @param size The sample size.
#' @param select If sampling from a specific group or list of groups.
#' @param replace Should sampling be with replacement?
#' @param both.sets If TRUE, both `sample` and `.data` are returned.
#'
#' @keywords Manipulation
#'
#' @examples
#' # Generate a couple of sample data.frames to play with
#'
#' set.seed(51)
#' dat1 <- data.frame(ID = 1:100, A = sample(c("AA", "BB", "CC", "DD", "EE"),
#' 100, replace = TRUE), B = rnorm(100), C = abs(round(rnorm(100), digits = 1)),
#' D = sample(c("CA", "NY", "TX"), 100, replace = TRUE), E = sample(c("M","F"),
#' 100, replace = TRUE))
#'
#' # Let's take a 10% sample from all -A- groups in dat1
#'  stratified(dat1, "A", 0.1)
#'
#'  # Let's take a 10% sample from only 'AA' and 'BB' groups from -A- in dat1
#'  stratified(dat1, "A", 0.1, select = list(A = c("AA", "BB")))
#'
#'  # Let's take 5 samples from all -D- groups in dat1, specified by column
#' stratified(dat1, group = 5, size = 5)
#'
#' # Let's take a sample from all -A- groups in dat1, where we specify the
#' # number wanted from each group
#' stratified(dat1, "A", size = c(3, 5, 4, 5, 2))
#'
#' # Use a two-column strata (-E- and -D-) but only interested in cases where
#' # -E- == 'M'
#' stratified(dat1, c("E", "D"), 0.15, select = list(E = "M"))
#'
#' @export
`stratified` <- function(.data, group, size, select = NULL,
                       replace = FALSE, both.sets = FALSE) {
  if (is.null(select)) {
    .data <- .data
  } else {
    if (is.null(names(select))) stop("'select' must be a named list")
    if (!all(names(select) %in% names(.data)))
      stop("Please verify your 'select' argument")
    temp <- sapply(names(select),
                   function(x) .data[[x]] %in% select[[x]])
    .data <- .data[rowSums(temp) == length(select), ]
  }
  .data.interaction <- interaction(.data[group], drop = TRUE)
  .data.table <- table(.data.interaction)
  .data.split <- split(.data, .data.interaction)
  if (length(size) > 1) {
    if (length(size) != length(.data.split))
      stop("Number of groups is ", length(.data.split),
           " but number of sizes supplied is ", length(size))
    if (is.null(names(size))) {
      n <- stats::setNames(size, names(.data.split))
      message(sQuote("size"), " vector entered as:\n\nsize = structure(c(",
              paste(n, collapse = ", "), "),\n.Names = c(",
              paste(shQuote(names(n)), collapse = ", "), ")) \n\n")
    } else {
      ifelse(all(names(size) %in% names(.data.split)),
             n <- size[names(.data.split)],
             stop("Named vector supplied with names ",
                  paste(names(size), collapse = ", "),
                  "\n but the names for the group levels are ",
                  paste(names(.data.split), collapse = ", ")))
    }
  } else if (size < 1) {
    n <- round(.data.table * size, digits = 0)
  } else if (size >= 1) {
    if (all(.data.table >= size) || isTRUE(replace)) {
      n <- stats::setNames(rep(size, length.out = length(.data.split)),
                    names(.data.split))
    } else {
      message(
        "Some groups\n---",
        paste(names(.data.table[.data.table < size]), collapse = ", "),
        "---\ncontain fewer observations",
        " than desired number of samples.\n",
        "All observations have been returned from those groups.")
      n <- c(sapply(.data.table[.data.table >= size], function(x) x = size),
             .data.table[.data.table < size])
    }
  }
  temp <- lapply(
    names(.data.split),
    function(x) .data.split[[x]][sample(.data.table[x],
                                     n[x], replace = replace), ])
  set1 <- do.call("rbind", temp)

  if (isTRUE(both.sets)) {
    set2 <- .data[!rownames(.data) %in% rownames(set1), ]
    list(SET1 = set1, SET2 = set2)
  } else {
    set1
  }
}### end -- stratified function
NULL
