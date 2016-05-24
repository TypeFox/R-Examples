 
##' header.t2t
##'
##' @keywords internal
##' @param caption caption
##' @param caption.level caption.level
header.t2t <- function(caption = NULL, caption.level = "") {
  res <- ""
  if (is.null(caption.level))
    caption.level <- ""
  
  if (!is.null(caption)) {
    if (is.numeric(caption.level) & caption.level > 0) { lev <- paste(rep("=", caption.level), collapse = "") ; res <- paste(lev, " ", caption, " ", lev, "\n\n", sep = "") }
    else if (caption.level == "s") res <- paste(beauty.t2t(caption, "s"), "\n\n", sep = "")
    else if (caption.level == "e") res <- paste(beauty.t2t(caption, "e"), "\n\n", sep = "")
    else if (caption.level == "m") res <- paste(beauty.t2t(caption, "m"), "\n\n", sep = "")
    else res <- paste(caption, "\n\n", sep = "")
  }
  return(res)
}

##' beauty.t2t
##'
##' @keywords internal
##' @param x x
##' @param beauti beauti
beauty.t2t <- function(x, beauti = c("e", "m", "s")) {
  x[is.na(x)] <- "NA"
  if (beauti == "s" | beauti == "h") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("\\*\\*.*\\*\\*", x)+1)/2) # bold seulement si != de "" et si pas de bold
    if (length(x[!y]) != 0) x[!y] <- sub("(^ *)([:alpha]*)", "\\1\\*\\*\\2", sub("([:alpha:]*)( *$)", "\\1\\*\\*\\2", x[!y]))
  }
  if (beauti == "e") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("//.*//", x)+1)/2) # it seulement si != de "" et si pas de it
    if (length(x[!y]) != 0) x[!y] <-sub("(^ *)([:alpha]*)", "\\1//\\2", sub("([:alpha:]*)( *$)", "\\1//\\2", x[!y]))
  }
  if (beauti == "m") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("``.*``", x)+1)/2) # it seulement si != de "" et si pas de mono
    if (length(x[!y]) != 0) x[!y] <-sub("(^ *)([:alpha]*)", "\\1``\\2", sub("([:alpha:]*)( *$)", "\\1``\\2", x[!y]))
  }
  return(x)
}

##' escape.t2t
##'
##' @keywords internal
##' @param x x
escape.t2t <- function(x) {
  xx <- gsub("\\|", "\\\\|", x)
  xx
}

##' show.t2t.table
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
show.t2t.table <- function(x, include.rownames = FALSE, include.colnames = FALSE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = FALSE, footer = FALSE, align = NULL, col.width = 1, style = NULL, lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", ...) {

  x <- escape.t2t(tocharac(x, include.rownames, include.colnames, rownames, colnames, format, digits, decimal.mark, na.print))
  nrowx <- nrow(x)
  ncolx <- ncol(x)
  
  if (!is.null(style)) {
    style <- expand(style, nrowx, ncolx)
    style[!(style %in% c("s", "e", "m"))] <- ""
    style[style == "s"] <- "**"
    style[style == "e"] <- "//"
    style[style == "m"] <- "``"
  } else {
    style <- ""
    style <- expand(style, nrowx, ncolx)
  }
  if (include.rownames & include.colnames) {
    style[1, 1] <- ""
  }
  
  before_cell_content <- after_cell_content <- style
  before_cell_content <- paste.matrix(" ", before_cell_content, sep = "", transpose.vector = TRUE)
  after_cell_content <- paste.matrix(after_cell_content, " ", sep = "", transpose.vector = TRUE)
  
  if (is.null(align)) {
    align <- "l"
  } 
  
  if (!is.null(frame)) {
    if (frame == "none") {
      frame <- 0
    } else {
      frame <- 1
    }
  } else {
    frame <- 0
  }
  
  line_separator <- FALSE
  
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
      x <- cbind(c(rep("", include.colnames), beauty.t2t(linelgroup[[i]], lstyle)), x)
    }
  }
  if (!is.null(rgroup)) {
    for (i in 1:length(rgroup)) {
      x <- cbind(x, c(rep("", include.colnames), beauty.t2t(linergroup[[i]], rstyle)))
    }
  }
  if (!is.null(tgroup)) {
    for (i in 1:length(tgroup)) {
      x <- rbind(c(rep("", include.rownames + length(lgroup)), beauty.t2t(linetgroup[[i]], tstyle), rep("", length(rgroup))), x)
    }    
  }
  if (!is.null(bgroup)) {
    for (i in 1:length(bgroup)) {
      x <- rbind(x, c(rep("", include.rownames + length(lgroup)), beauty.t2t(linebgroup[[i]], bstyle), rep("", length(rgroup))))
    }
  }

  vsep <- expand("|", nrowx+length(tgroup)+length(bgroup), ncolx+length(lgroup)+length(rgroup)+frame)
  
  if (!is.null(tgroup) | !is.null(bgroup)) {
    if (frame == 0) {
      vsep <- expand(vsep, nrow(vsep), ncol(vsep)+1, what = "")
    }
  }

  if (!is.null(tgroup)) {
    for (i in length(tgroup):1) {
      vsep[i, length(lgroup)+include.rownames+cumsum(n.tgroup[[i]])+1] <- sapply(n.tgroup[[i]], function(x) paste(rep("|", x), collapse = ""))
      vsep[i, -(length(lgroup)+include.rownames+cumsum(n.tgroup[[i]])+1)] <- ""
      vsep[i, 1] <- "|"
    }
  }
  if (!is.null(bgroup)) {
    for (i in 1:length(bgroup)) {
      vsep[nrowx+length(tgroup)+i, length(lgroup)+include.rownames+cumsum(n.bgroup[[i]])+1] <- sapply(n.bgroup[[i]], function(x) paste(rep("|", x), collapse = ""))
      vsep[nrowx+length(tgroup)+i, -(length(lgroup)+include.rownames+cumsum(n.bgroup[[i]])+1)] <- ""
      vsep[nrowx+length(tgroup)+i, 1] <- ""
    }
  }

  topleftrow <- 0 + include.colnames + length(tgroup)
  topleftcol <- 0 + include.rownames + length(lgroup)
  if (topleftrow > 0 & topleftcol > 1) {
    vsep[1:topleftrow, topleftcol+1] <- paste(rep("|", topleftcol), collapse = "")
    vsep[1:topleftrow, 2:topleftcol] <- ""
  }

  bottomleftrow <- 0 + length(bgroup)
  bottomleftcol <- 0 + include.rownames + length(lgroup)
  if (bottomleftrow > 0 & bottomleftcol > 0) {
    vsep[nrow(vsep):(nrow(vsep)-bottomleftrow+1), bottomleftcol+1] <- paste(rep("|", bottomleftcol), collapse = "")
    vsep[nrow(vsep):(nrow(vsep)-bottomleftrow+1), 2:bottomleftcol] <- ""
  }

  toprightrow <- 0 + include.colnames + length(tgroup)
  toprightcol <- 0 + length(rgroup)
  if (toprightrow > 0 & toprightcol > 1) {
    vsep[1:toprightrow, ncol(vsep)] <- paste(rep("|", toprightcol), collapse = "")
    vsep[1:toprightrow, (ncol(vsep)-1):(ncol(vsep)-toprightcol+1)] <- ""
  }

  bottomrightrow <- 0 + length(bgroup)
  bottomrightcol <- 0 + length(rgroup)
  if (bottomrightrow > 0 & bottomrightcol > 1) {
    if (frame == 0) {
      vsep <- expand(vsep, nrow(vsep), ncol(vsep)+1, what = "")
    }
    vsep[nrow(vsep):(nrow(vsep)-bottomleftrow+1), ncol(vsep)] <- paste(rep("|", bottomrightcol), collapse = "")
    vsep[nrow(vsep):(nrow(vsep)-bottomleftrow+1), (ncol(vsep)-1):(ncol(vsep)-bottomrightcol+1)] <- ""
  }
  
  results <- print.character.matrix(x, line_separator = line_separator, vsep = vsep, before_cell_content = before_cell_content, after_cell_content = after_cell_content, justify = align, space = 1, print = FALSE)

  headfoot <- NULL
  if (header + footer > nrowx) {
    headfoot <- 1:nrowx
  } 
  if (is.logical(header) & header)
    header <- 1
  if (header > 0) {
    header <- 1:min(c(header, nrowx)) + length(tgroup)
  }
  if (is.logical(footer) & footer)
    footer <- nrowx
  if (footer > 0) {
    footer <- nrowx:(nrowx+1-min(c(footer, nrowx))) + length(tgroup)
  }
  headfoot <- unique(c(header, footer))
  results[headfoot] <- paste("|", results[headfoot], sep = "")
  
  cat(header.t2t(caption = caption, caption.level = caption.level))
  cat(results, sep = "\n")
}

##' show.t2t.list
##'
##' @keywords internal
##' @param x x
##' @param caption caption
##' @param caption.level caption.level
##' @param list.type list.type
##' @param ... ...
show.t2t.list <- function(x, caption = NULL, caption.level = NULL, list.type = "bullet", ...) {
  indent.mark <- " "
  if (list.type == "bullet") mark <- rep("-", length(x))
  if (list.type == "number") mark <- rep("+", length(x))
  if (list.type == "none")  { mark <- rep("", length(x)); indent.mark = "" }
  if (list.type == "label") {
    if (is.null(names(x))) {
      namesx <- paste("[[", 1:length(x), "]]", sep = "")
    } else {
      namesx <- names(x)
    }
    mark <- paste(": ", namesx, "\n", sep = "")
    indent.mark = ""
  }
  
  charac.x <- vector("character", length(x))
  for (i in 1:length(x)) {
    tmp <- x[[i]]
    if (list.type == "label") tmp <- sub("^\t*", "", tmp)
    tmp <- gsub('\t|(*COMMIT)(*FAIL)', indent.mark, tmp, perl = TRUE)
    charac.x[i] <- sub("(^ *)", paste("\\1", mark[i], indent.mark, sep = ""), tmp)
  }
  cat(header.t2t(caption = caption, caption.level = caption.level))
  cat(charac.x, sep = "\n")
}
