# From Hmisc package

##' @export
##' @method ascii describe.single
ascii.describe.single <- function (x, condense = TRUE, ...) {
   wide <- .Options$width
    # des : le titre
    des <- x$descript
    if (length(x$units))
      des <- paste(des, " [", x$units, "]", sep = "")
        if (length(x$format))
          des <- paste(des, "  Format:", x$format, sep = "")
            dim.counts <- dim(x$count)
            if (is.null(dim.counts)) {counts <-as.character(x$count)} else {
            counts <- matrix(as.character(x$count), dim.counts[1], dim.counts[2])}
            names(counts) <- names(x$count)
            counts <- ascii(counts, include.colnames = TRUE, caption = des, caption.level = "s")
            val <- x$values
            if (length(val)) {
              if (!is.matrix(val)) {
                if (length(val) != 10 || !all(names(val) == c("L1",
                        "L2", "L3", "L4", "L5", "H5", "H4", "H3", "H2",
                        "H1"))) {
                  cat("\n")
                    val <- paste(names(val), ifelse(val > 1, paste(" (",
                            val, ")", sep = ""), ""), sep = "")
                    val <- strwrap(val, exdent = 4)
                    val <- as.list(sub('(^    )(.*)', '\t\\2', val))
                    val <- ascii(val, list.type = "none")
                }
                else {
                  if (condense) {
                    low <- paste("lowest:", paste(val[1:5], collapse = " "))
                      hi <- paste("highest:", paste(val[6:10], collapse = " "))
                      if (nchar(low) + nchar(hi) + 2 > wide)
                        val <- as.list(c(low, hi))
                      else val <- as.list(paste(low, hi, sep = ", "))
                    val <- ascii(val, list.type = "none")
                  }
                  else {
                    dim.val <- dim(val)
                    if (is.null(dim.val)) {val <- as.character(val)} else {
                    val <- matrix(as.character(val), dim.val[1], dim.val[2])}
                    names(val) <- names(x$values)
                    val <- ascii(val, include.colnames = TRUE)
                  }
                }
              }
              else {
                lev <- dimnames(val)[[2]]
                  if (condense && (mean(nchar(lev)) > 10 | length(lev) <
                        5)) {
                    z <- ""
                      len <- 0
                      for (i in 1:length(lev)) {
                        w <- paste(lev[i], " (", val[1, i], ", ", val[2,
                            i], "%)", sep = "")
                          if (i == 1) z <- w
                          else z <- paste(z, w, sep = ", ")
                  }
                  val <- ascii(as.list(z), list.type = "none")
                  }
                  else {
                    dim.val <- dim(val)
                    if (is.null(dim.val)) {val <- as.character(val)} else {
                    val <- matrix(as.character(val), dim.val[1], dim.val[2])}
                    rownames(val) <- rownames(x$values)
                    colnames(val) <- colnames(x$values)
                    val <- ascii(val, include.rownames = TRUE, include.colnames = TRUE)
                  }
              }
            }
#   if (length(x$mChoice)) {
#       print(x$mChoice, prlabel = FALSE)
#   }
  res <- asciiMixed$new(args = list(counts, val))
  return(res)
}

##' @param condense default is TRUE to condense the output with regard to the 5
##'   lowest and highest values and the frequency table (\code{describe()} in
##'   package \code{Hmisc}).
##' @export
##' @method ascii describe
##' @rdname ascii
ascii.describe <- function (x, condense = TRUE, ...) {
  at <- attributes(x)
  descrip <- ifelse(is.null(at$descript), "", at$descrip)
  if (is.null(at$dimensions[2])) {
    variable <- NULL
  } else {
    variable <- paste(at$dimensions[2], "Variable")
  }
  if (is.null(at$dimensions[1])) {
    observation <- NULL
  } else {
    observation <- paste(at$dimensions[1], "Observations")
  }
  if (!is.null(variable) | !is.null(observation) | descrip != "") {
    des <- ascii(list(variable, observation), caption = descrip, caption.level = NULL)
  } else {des <- NULL}
    if (length(at$dimensions)) {
      xx <- lapply(x, ascii.describe.single, condense = condense)
        res <- NULL
        for (z in 1:length(x)) {
          if (length(x[[z]]) == 0) next
          res <- asciiMixed$new(args = list(res, xx[[z]]))
        }
    }
    else res <- ascii.describe.single(x, condense = condense)

  if (length(at$naprint)) { na <- ascii(as.list(at$naprint)) ; res <- asciiMixed$new(args = list(des, na, res)) }
  else res <- asciiMixed$new(args = list(des, res))

  return(res)
}

# ascii.summary.formula.response <- function (x, vnames = c("labels", "names"), prUnits = TRUE, abbreviate.dimnames = FALSE,
#     prefix.width, min.colwidth, formatArgs = NULL, ...)
# {
#   stats <- x
#     stats <- oldUnclass(stats)
#     vnames <- match.arg(vnames)
#     ul <- vnames == "labels"
#     at <- attributes(stats)
#     ns <- length(at$strat.levels)
#     vlabels <- at$labels
#     if (prUnits) {
#       atu <- translate(at$units, "*", " ")
#         vlabels <- ifelse(atu == "", vlabels, paste(vlabels,
#               " [", atu, "]", sep = ""))
#     }
#   entete <- list(paste(at$ylabel, if (ns > 1)
#       paste(" by", if (ul)
#         at$strat.label
#         else at$strat.name), "    N=", at$n, if (at$nmiss)
#       paste(", ", at$nmiss, " Missing", sep = ""),
#       sep = ""))
#   entete <- ascii(entete, list.type = "none")
#     d <- dim(stats)
#     if (exists("print.char.matrix")) {
#       nr <- length(at$nlevels)
#         vlab <- if (ul)
#         at$vlabel
#         else at$vname
# #           z <- matrix("", nrow = nr, ncol = 1 + d[2], dimnames = list(vlab,
# #                 NULL))
#             dz <- dimnames(stats)[[1]]
#             cstats <- matrix("", nrow = d[1], ncol = d[2])
#             for (j in 1:d[2]) {
#               ww <- c(list(stats[, j]), formatArgs)
#                 cstats[, j] <- do.call("format", ww)
#                 cstats[is.na(stats[, j]), j] <- ""
#             }
#             z <- cbind(vlab, at$dimnames[[1]], cstats)
#             dimnames(z) <- list(NULL, c("", "", dimnames(stats)[[2]]))
#
# #         is <- 1
# #           for (i in 1:nr) {
# #             ie <- is + at$nlevels[i] - 1
# #               z[i, 1] <- paste(dz[is:ie], collapse = "\n")
# #               for (j in 1:d[2]) z[i, j + 1] <- paste(cstats[is:ie,
# #                   j], collapse = "\n")
# #                 is <- ie + 1
# #           }
# #         if (missing(prefix.width))
# #           prefix.width <- max(nchar(dimnames(z)[[1]]))
# #             if (missing(min.colwidth))
# #               min.colwidth <- max(min(nchar(cstats)[nchar(cstats) >
# #                     0]), min(nchar(dimnames(stats)[[2]])))
# #                 z <- rbind(c("", dimnames(stats)[[2]]), z)
# #                  print.char.matrix(z, col.names = FALSE, ...)
#                  tab <- ascii(z, include.colnames = TRUE)
#     }
# #   dz <- if (length(at$strat.levels) == 1)
# #     dimnames(stats)[[2]]
# #     else paste(rep(at$strat.levels, length = d[2]), dimnames(stats)[[2]],
# #         sep = ":")
# #       z <- matrix("", ncol = d[2] + 2, nrow = d[1], dimnames = list(rep("",
# #               d[1]), c("", "", dz)))
# #         z[, 1] <- if (ul)
# #         vlabels
# #         else at$vname
# #           z[, 2] <- dimnames(stats)[[1]]
# #             for (i in 1:d[2]) {
# #               ww <- c(list(stats[, i]), formatArgs)
# #                 z[, i + 2] <- do.call("format", ww)
# #             }
# #         print(z, quote = FALSE)
# #           invisible()
#   res <- asciiMixed$new(entete, tab)
#   class(res) <- c("ascii", "proto", "environment")
#   return(res)
#
# }
#
