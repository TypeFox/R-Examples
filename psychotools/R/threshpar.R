## threshold generic
threshpar <- function (object, ...) {
  UseMethod("threshpar")
}


## methods for class 'thrsldpar'
print.threshpar <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  if (attr(x, "relative")) cat("Relative item ") else cat("Item ")
  cat("response threshold parameters (", attr(x, "model"), ", type ", attr(x, "type"), "):\n", sep = "")
  print(coef(x), digits = digits, ...)
  invisible(x)
}

coef.threshpar <- function (object, type = c("vector", "matrix", "list"), ...) {

  ## check input
  type <- match.arg(type)

  ## create requested result value
  if (type == "vector") {
    rv <- unlist(object)
    if (!is.null(rv)) names(rv) <- gsub("(.*)\\.(.*)", "\\1-\\2", names(rv))
  } else if (type == "matrix") {
    oj <- sapply(object, length)
    if (is.logical(attr(object, "alias"))) {
      rv <- matrix(NA, nrow = length(object), ncol = max(oj), dimnames = list(names(object), paste0("C", 1:max(oj))))
      for (j in 1:length(object)) rv[j, 1:oj[j]] <- object[[j]]
    } else {
      alias <- attr(object, "alias")
      if (attr(object, "relative")) {
        oj <- oj + 1
        rv <- matrix(NA, nrow = length(object), ncol = max(oj), dimnames = list(names(object), paste0("C", 1:max(oj))))
        for (j in 1:length(object)) rv[j, setdiff(1:oj[j], alias[[j]])] <- object[[j]]
      } else {
        alias <- strsplit(alias, "-")[[1]]
        item <- as.integer(gsub(pattern = "I", replacement = "", x = alias[1]))
        cat <- as.integer(gsub(pattern = "C", replacement = "", x = alias[2]))

        if (attr(object, "model") != "RM") {
          oj[item] <- oj[item] + 1
          rv <- matrix(NA, nrow = length(object), ncol = max(oj), dimnames = list(names(object), paste0("C", 1:max(oj))))
          for (j in 1:length(object)) {
            if (j == item) rv[j, setdiff(1:oj[j], cat)] <- object[[j]]
            else rv[j, 1:oj[j]] <- object[[j]]
          }
        } else { ## in RM: ref for absolute threshold parameters means one item is dismissed
          m <- length(object) + 1
          nms <- vector(mode = "character", length = m)
          nms[setdiff(1:m, item)] <- names(object)
          nms[item] <- names(attr(object, "alias"))
          rv <- matrix(NA, nrow = m, ncol = 1, dimnames = list(nms, "C1"))
          rv[setdiff(1:m, item), 1] <- unlist(object)
        }
      }
    }
  } else {
    rv <- object
  }

  ## remove all attributes beside names, then return named threshold parameters
  return(rv)
}

vcov.threshpar <- function (object, ...) attr(object, "vcov")
