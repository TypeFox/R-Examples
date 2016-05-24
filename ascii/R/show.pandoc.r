##' beauty.pandoc
##'
##' @keywords internal
##' @param x x
##' @param beauti beauti
beauty.pandoc <- function(x, beauti = c("e", "m", "s")) {
  x[is.na(x)] <- "NA"
  if (beauti == "s") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("\\*\\*.*\\*\\*", x)+1)/2) # bold seulement si != de "" et si pas de bold
    if (length(x[!y]) != 0) x[!y] <- sub("(^ *)([:alpha]*)", "\\1\\*\\*\\2", sub("([:alpha:]*)( *$)", "\\1\\*\\*\\2", x[!y]))
    if (length(x[y]) != 0) x[y] <- sub("(^  *$)", "\\1 ", x[y]) # rajouter suffisamment d'espaces lorsque la case est vide pour l'alignement globale
  }
  if (beauti == "e") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("\\*.*\\*", x)+1)/2) # it seulement si != de "" et si pas de it
    if (length(x[!y]) != 0) x[!y] <-sub("(^ *)([:alpha]*)", "\\1\\*\\2", sub("([:alpha:]*)( *$)", "\\1\\*\\2", x[!y]))
    if (length(x[y]) != 0) x[y] <- sub("(^ *$)", "\\1 ", x[y]) # rajouter suffisamment d'espaces lorsque la case est vide pour l'alignement globale
  }
  if (beauti == "m") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("`.*`", x)+1)/2) # it seulement si != de "" et si pas de mono
    if (length(x[!y]) != 0) x[!y] <-sub("(^ *)([:alpha]*)", "\\1`\\2", sub("([:alpha:]*)( *$)", "\\1`\\2", x[!y]))
    if (length(x[y]) != 0) x[y] <- sub("(^ *$)", "\\1 ", x[y]) # rajouter suffisamment d'espaces lorsque la case est vide pour l'alignement globale
  }
  return(x)
}

##' header.pandoc
##'
##' @keywords internal
##' @param caption caption
##' @param caption.level caption.level
header.pandoc <- function(caption = NULL, caption.level = NULL) {
  res <- ""
  if (is.null(caption.level))
    caption.level <- "s"
  if (!is.null(caption)) {
    if (is.numeric(caption.level) & caption.level > 0) {
      res <- paste(paste(rep("#", caption.level), collapse = ""), caption, "\n")
    } else if (is.character(caption.level) & caption.level %in% c("s", "e", "m")) {
      if (caption.level == "s")
        res <- paste(beauty.pandoc(caption, "s"), "\n", sep = "")
      else if (caption.level == "e")
        res <- paste(beauty.pandoc(caption, "e"), "\n", sep = "")
      else if (caption.level == "m")
        res <- paste(beauty.pandoc(caption, "m"), "\n", sep = "")
    } else if (caption.level == "none")
      res <- paste(caption, "\n", sep = "")
  }
  return(res)
}


## ##' escape.pandoc
## ##'
## ##' @keywords internal
## ##' @param x x
## escape.pandoc <- function(x) {
##   xx <- gsub("\\|", " \\\\vert ", x)
##   xx
## }


##' show.pandoc.table
##'
##' @keywords internal
##' @param x x
##' @param include.rownames include.rownames 
##' @param include.colnames include.colnames 
##' @param rownames rownames 
##' @param colnames colnames 
##' @param format format 
##' @param digits digits 
##' @param decimal.mark decimal.mark 
##' @param na.print na.print 
##' @param caption caption 
##' @param caption.level 
##' @param width width 
##' @param frame frame 
##' @param grid grid 
##' @param valign valign 
##' @param header header 
##' @param footer footer 
##' @param align align 
##' @param col.width col.width 
##' @param style style 
##' @param lgroup lgroup 
##' @param n.lgroup n.lgroup 
##' @param lalign lalign 
##' @param lvalign lvalign 
##' @param lstyle lstyle 
##' @param rgroup rgroup 
##' @param n.rgroup n.rgroup 
##' @param ralign ralign 
##' @param rvalign rvalign 
##' @param rstyle rstyle 
##' @param tgroup tgroup 
##' @param n.tgroup n.tgroup 
##' @param talign talign 
##' @param tvalign tvalign 
##' @param tstyle tstyle 
##' @param bgroup bgroup
##' @param n.bgroup n.bgroup 
##' @param balign balign 
##' @param bvalign bvalign 
##' @param bstyle bstyle 
##' @param ... ...
show.pandoc.table <- function(x, include.rownames = FALSE, include.colnames = FALSE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = FALSE, footer = FALSE, align = NULL, col.width = 1, style = NULL, lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", ...) {

  x <- tocharac(x, include.rownames, include.colnames, rownames, colnames, format, digits, decimal.mark, na.print)
  nrowx <- nrow(x)
  ncolx <- ncol(x)
  
  if (!is.null(style)) {
    style <- expand(style, nrowx, ncolx)
    style[!(style %in% c("s", "e", "m"))] <- ""
    style[style == "s"] <- "**"
    style[style == "e"] <- "*"
    style[style == "m"] <- "`"
  } else {
    style <- ""
    style <- expand(style, nrowx, ncolx)
  }
  if ((is.logical(header) & header) | header >= 1) {
    style[1:min(nrow(style), header), ] <- "**"
    header <- 1
  }
  if (footer >= 1) {
    style[nrow(style):max((nrow(style)-footer+1), 1), ] <- "**"
  }
  if (include.rownames & include.colnames) {
    style[1, 1] <- ""
  }

  before_cell_content <- after_cell_content <- style

  x <- paste.matrix(before_cell_content, x, after_cell_content, sep = "", transpose.vector = TRUE)
  
  if (tstyle == "h")
    tstyle <- "s"
  if (bstyle == "h")
    bstyle <- "s"
  if (rstyle == "h")
    rstyle <- "s"
  if (lstyle == "h")
    lstyle <- "s"

  # groups
  if (!is.null(lgroup)) {
    if (!is.list(lgroup))
      lgroup <- list(lgroup)
    n.lgroup <- groups(lgroup, n.lgroup, nrowx-include.colnames)[[2]]
    linelgroup <- linegroup(lgroup, n.lgroup)
  }
  if (!is.null(rgroup)) {
    if (!is.list(rgroup))
      rgroup <- list(rgroup)
    n.rgroup <- groups(rgroup, n.rgroup, nrowx-include.colnames)[[2]]
    linergroup <- linegroup(rgroup, n.rgroup)
  }
  if (!is.null(tgroup)) {
    if (!is.list(tgroup))
      tgroup <- list(tgroup)
    n.tgroup <- groups(tgroup, n.tgroup, ncolx-include.rownames)[[2]]
    linetgroup <- linegroup(tgroup, n.tgroup)
  }
  if (!is.null(bgroup)) {
    if (!is.list(bgroup))
      bgroup <- list(bgroup)
    n.bgroup <- groups(bgroup, n.bgroup, ncolx-include.rownames)[[2]]
    linebgroup <- linegroup(bgroup, n.bgroup)
  }

  if (!is.null(lgroup)) {
    for (i in 1:length(lgroup)) {
      x <- cbind(c(rep("", include.colnames), beauty.pandoc(linelgroup[[i]], lstyle)), x)
    }
  }
  if (!is.null(rgroup)) {
    for (i in 1:length(rgroup)) {
      x <- cbind(x, c(rep("", include.colnames), beauty.pandoc(linergroup[[i]], rstyle)))
    }
  }
  if (!is.null(tgroup)) {
    for (i in 1:length(tgroup)) {
      x <- rbind(c(rep("", include.rownames + length(lgroup)), beauty.pandoc(linetgroup[[i]], tstyle), rep("", length(rgroup))), x)
    }
  }
  if (!is.null(bgroup)) {
    for (i in 1:length(bgroup)) {
      x <- rbind(x, c(rep("", include.rownames + length(lgroup)), beauty.pandoc(linebgroup[[i]], bstyle), rep("", length(rgroup))))
    }
  }
  
  line_separator <- FALSE
  line_separator_pos <- NULL

  line_sep <- sapply(apply(x, 2, function(x) max(nchar(x))), function(x) paste(rep("-", x), collapse = ""))

  if (is.null(align)) {
    align <- "l"
  }
  if (is.null(lalign) & length(lgroup) > 0) {
    lalign <- "c"
  } else if (length(lgroup) == 0) {
    lalign <- NULL
  }
  if (is.null(ralign) & length(rgroup) > 0) {
    ralign <- "c"
  } else if (length(rgroup) == 0) {
    ralign <- NULL
  }
  
  align <- c(rep(lalign, length(lgroup)), expand(align, 1, ncolx), rep(ralign, length(rgroup)))
  justify <- align
  justify[justify == "l"] <- "left"
  justify[justify == "r"] <- "right"
  justify[justify == "c"] <- "centre"
  for (i in seq_along(align)) {
    if (align[i] == "l") {
      line_sep[i] <- paste(line_sep[i], "--", sep = "")
    }
    if (align[i] == "r") {
      line_sep[i] <- paste("--", line_sep[i], sep = "")
    }
    if (align[i] == "c") {
      line_sep[i] <- paste("--", line_sep[i], "--", sep = "")
    }
  }
    
  if (header > 0) {
    x <- rbind(x[1, ], line_sep, x[-1, ])
  } else {
    x <- rbind(line_sep, x)
  }
  
  x <- rbind(x, line_sep)

  for (i in 1:ncol(x)) {
    x[, i] <- as.character(format(x[, i], trim = TRUE, justify = justify[i]))
  }

  results <- paste.matrix(x, collapse = " ")

  cat("\n")
  cat(results, sep = "\n")
  cat("\n")
  if (!is.null(caption))
    cat("Table: ", caption, "\n\n", sep = "")
}


##' show.pandoc.list
##'
##' @keywords internal
##' @param x x
##' @param caption caption
##' @param caption.level caption.level
##' @param list.type list.type
##' @param ... ...
show.pandoc.list <- function(x, caption = NULL, caption.level = NULL, list.type = "bullet", ...) {
  indent.mark <- "    "
  if (list.type == "bullet") mark <- rep("-", length(x))
  if (list.type == "number") mark <- rep("#.", length(x)) #paste(seq(1, length(x), 1), "#.", sep = "")
  if (list.type == "none")  { mark <- rep("", length(x)); indent.mark = ""}
  if (list.type == "label") {
    if (is.null(names(x))) {
      namesx <- paste("[[", 1:length(x), "]]", sep = "")
    } else {
      namesx <- names(x)
    }
    mark <- paste(namesx, "\n  ~")
    x <- lapply(x, function(x) gsub("\t", "", x))
  }
  
  charac.x <- vector("character", length(x))
  for (i in 1:length(x)) {
    tmp <- x[[i]]
    tmp <- gsub('\t|(*COMMIT)(*FAIL)', indent.mark, tmp, perl = TRUE)
    charac.x[i] <- sub("(^ *)", paste("\\1", mark[i], " ", sep = ""), tmp)
  }
  cat(header.pandoc(caption = caption, caption.level = caption.level))
  cat("\n")
  cat(charac.x, sep = "\n")
}
