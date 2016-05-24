##' Compute a contingency table
##'
##' @importFrom Hmisc label
##' @import ascii
##' @param x factor
##' @param y factor
##' @param margin margin
##' @param useNA useNA
##' @param propNA propNA
##' @param addmargins addmargins
##' @author David Hajage
##' @keywords internal
tabular <- function(x, y, margin = 0:2, useNA = c("no", "ifany", "always"), propNA = TRUE, addmargins = FALSE, test = FALSE, test.tabular = test.tabular.auto, show.test = display.test, plim = 4, show.method = TRUE) {
  n <- n.table(x, y, useNA = useNA, margin = margin, addmargins = addmargins, test = test, test.tabular = test.tabular, show.test = show.test, plim = plim, show.method = show.method, na = propNA)
  p <- p.table(x, y, useNA = useNA, propNA = propNA, margin = margin, addmargins = addmargins)
  rn <- rownames(n)
  rn[is.na(rn)] <- "NA"
  results <- c(n = list(n), p)
  results <- do.call(interleave.matrix, results)
  # remove unnecessary rows (all NA)
  results <- results[apply(results, 1, function(x) any(!is.na(x))), ]

  attr(results, "lgroup") <- list(gsub("(^n|^cell|^row|^col)(\\.)", "\\1", gsub("(^n\\.|^cell\\.|^row\\.|^col\\.)(.+$)", "\\1", rownames(results))), rn)
  attr(results, "n.lgroup") <- list(1, table(gsub("(^n\\.|^cell\\.|^row\\.|^col\\.)(.+$)", "\\2", rownames(results)))[rn])
  attr(results, "tgroup") <- NULL
  attr(results, "n.tgroup") <- NULL
  attr(results, "rgroup") <- attr(n, "test")
  attr(results, "n.rgroup") <- nrow(results)
  
  class(results) <- c("tabular", "matrix")
  results
}

##' Compute a contingency table (data.frame input)
##'
##' @importFrom Hmisc label
##' @param dfx data.frame
##' @param dfy data.frame
##' @param margin margin
##' @param useNA useNA
##' @param propNA propNA
##' @author David Hajage
##' @keywords internal
tabular.data.frame <- function(dfx, dfy, margin = 0:2, useNA = c("no", "ifany", "always"), propNA = TRUE, addmargins = FALSE, test = FALSE, test.tabular = test.tabular.auto, show.test = display.test, plim = 4, show.method = TRUE, label = FALSE) {
  results <- lapply(dfy, function(y) lapply(dfx, tabular, y, margin = margin, useNA = useNA, propNA = propNA, addmargins = addmargins, test = test, test.tabular = test.tabular, show.test = show.test, plim = plim, show.method = show.method))

  ## noms <- names(results[[1]])
  if (!label)
    noms <- names(dfx)
  else
    noms <- sapply(dfx, Hmisc:::label.default)
  
  
  rgroup <- lapply(results, function(x) sapply(x, attr, "rgroup"))
  n.rgroup <- lapply(results, function(x) sapply(x, attr, "n.rgroup"))
  
  lgroup <- lapply(results[[1]], function(x) attr(x, "lgroup"))

  n.lgroup <- lapply(results[[1]], function(x) attr(x, "n.lgroup"))
  for (i in 1:length(n.lgroup)) {
    n.lgroup[[i]] <- c(n.lgroup[[i]], sum(n.lgroup[[i]][[2]]))
  }
  n.lgroup1 <- unlist(lapply(n.lgroup, function(x) x[[1]]))
  n.lgroup2 <- unlist(lapply(n.lgroup, function(x) x[[2]]))
  n.lgroup3 <- unlist(lapply(n.lgroup, function(x) x[[3]]))
  n.lgroup <- list(n.lgroup1, n.lgroup2, n.lgroup3)

  for (i in 1:length(lgroup)) {
    lgroup[[i]] <- c(lgroup[[i]], noms[i])
  }
  lgroup1 <- unlist(lapply(lgroup, function(x) x[[1]]))
  lgroup2 <- unlist(lapply(lgroup, function(x) x[[2]]))
  lgroup3 <- unlist(lapply(lgroup, function(x) x[[3]]))
  lgroup <- list(lgroup1, lgroup2, lgroup3)

  ## tgroup <- names(results)
  if (!label)
    tgroup <- names(dfy)
  else
    tgroup <- sapply(dfy, Hmisc:::label.default)
  
  n.tgroup <- unlist(lapply(results, function(x) ncol(x[[1]])))

  results <- lapply(results, rbind.list)

  attr(results[[1]], "lgroup") <- lgroup
  attr(results[[1]], "n.lgroup") <- n.lgroup

  if(test) {
    for (i in 1:length(results)) {
      attr(results[[i]], "rgroup") <- rgroup[[i]]
      attr(results[[i]], "n.rgroup") <- n.rgroup[[i]]
    }
  }
  
  for (i in 1:length(results)) {
    attr(results[[i]], "tgroup") <- tgroup[i]
    attr(results[[i]], "n.tgroup") <- n.tgroup[i]
  }
  attr(results, "dfx") <- dfx
  attr(results, "dfy") <- dfy
  
  class(results) <- "tabular"
  results
}

##' Ascii for tabular object.
##'
##' Ascii method for tabular object (internal).
##'
##' @export
##' @method ascii tabular
##' @import ascii 
##' @param x a tabular object
##' @param format see \code{?ascii} in \code{ascii} package
##' @param digits see \code{?ascii} in \code{ascii} package
##' @param include.rownames see \code{?ascii} in \code{ascii} package
##' @param include.colnames see \code{?ascii} in \code{ascii} package
##' @param header see \code{?ascii} in \code{ascii} package
##' @param rstyle see \code{?ascii} in \code{ascii} package
##' @param caption see \code{?ascii} in \code{ascii} package
##' @param caption.level see \code{?ascii} in \code{ascii} package
##' @param ... other arguments passed to \code{ascii}
##' @author David Hajage
##' @keywords univar
ascii.tabular <- function(x, format = "nice", digits = 5, include.rownames = FALSE, include.colnames = TRUE, header = TRUE, rstyle = "d", caption = NULL, caption.level = NULL, ...) {
  do.call(cbind.ascii, c(lapply(x, function(x) {
    ascii(x, format = format, digits = digits, include.rownames = include.rownames, include.colnames = include.colnames, header = header, lgroup = attr(x, "lgroup"), n.lgroup = attr(x, "n.lgroup"), tgroup = attr(x, "tgroup"), n.tgroup = attr(x, "n.tgroup"), rgroup = attr(x, "rgroup"), n.rgroup = attr(x, "n.rgroup"), rstyle = rstyle, ...)}), caption = caption, caption.level = caption.level))
}

##' Print tabular object.
##'
##' Print tabular object (internal).
##'
##' @export
##' @method print tabular
##' @import ascii
##' @param x a tabular object
##' @param type type of output (see \code{?ascii} in \code{ascii}
##' package)
##' @param lstyle see \code{?ascii} in \code{ascii} package
##' @param tstyle see \code{?ascii} in \code{ascii} package
##' @param ... other arguments passed to \code{ascii}
##' @author David Hajage
##' @keywords univar
print.tabular <- function(x, type = "rest", lstyle = "", tstyle = "", ...) {
  print(ascii.tabular(x, lstyle = lstyle, tstyle = tstyle, ...), type = type)
  ## invisible(x)
}

##' as.data.frame for tabular object.
##'
##' as.data.frame for tabular object (internal).
##'
##' @export
##' @param x a tabular object
##' @param ... not used
##' @author David Hajage
##' @keywords internal
as.data.frame.tabular <- function(x, ...) {
  xx <- do.call("cbind", x)  
  stat <- attr(x[[1]], "lgroup")[[1]]
  levels <- unlist(mapply(rep, attr(x[[1]], "lgroup")[[2]], attr(x[[1]], "n.lgroup")[[2]], SIMPLIFY = FALSE))
  var <- rep(attr(x[[1]], "lgroup")[[3]], attr(x[[1]], "n.lgroup")[[3]])
  data.frame(var = var, levels = levels, stat = stat, xx, row.names = NULL, check.names = FALSE)
}

##' Test if \code{x} is an tabular object
##'
##' @param x a tabular object
##' @author David Hajage
##' @keywords internal
is.tabular <- function(x)
  inherits(x, "tabular")
