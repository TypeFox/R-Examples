##' vsep.asciidoc
##'
##' @keywords internal
##' @param align align
##' @param valign valign
##' @param style style
vsep.asciidoc <- function (align = NULL, valign = NULL, style = NULL) {
  if (is.null(align) & is.null(valign) & is.null(style)) {
    res <- ""
  } else {
    if (is.null(align))
      align <- ""
    if (is.null(valign))
      valign <- ""
    if (is.null(style))
      style <- ""
    align[align == "l"] <- "<"
    align[align == "c"] <- "^"
    align[align == "r"] <- ">"
    
    valign[valign == "top"] <- ".<"
    valign[valign == "middle"] <- ".^"
    valign[valign == "bottom"] <- ".>"
    
    res <- paste.matrix(align, valign, style, sep = "")
  }
  return(res)
}


##' beauty.asciidoc
##'
##' @keywords internal
##' @param x x
##' @param beauti beauti
beauty.asciidoc <- function(x, beauti = c("e", "m", "s")) {
  x[is.na(x)] <- "NA"
  if (beauti == "s") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("\\*.*\\*", x)+1)/2) # bold seulement si != de "" et si pas de bold
    if (length(x[!y]) != 0) x[!y] <- sub("(^ *)([:alpha]*)", "\\1\\*\\2", sub("([:alpha:]*)( *$)", "\\1\\*\\2", x[!y]))
  }
  if (beauti == "e") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("'.*'", x)+1)/2) # it seulement si != de "" et si pas de it
    if (length(x[!y]) != 0) x[!y] <-sub("(^ *)([:alpha]*)", "\\1'\\2", sub("([:alpha:]*)( *$)", "\\1'\\2", x[!y]))
  }
  if (beauti == "m") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("\\+.*\\+", x)+1)/2) # it seulement si != de "" et si pas de mono
    if (length(x[!y]) != 0) x[!y] <-sub("(^ *)([:alpha]*)", "\\1\\+\\2", sub("([:alpha:]*)( *$)", "\\1\\+\\2", x[!y]))
  }
  return(x)
}


##' header.asciidoc
##'
##' @keywords internal
##' @param caption caption
##' @param caption.level caption.level
##' @param frame frame
##' @param grid grid
##' @param col.width col.width
##' @param width width
header.asciidoc <- function (caption = NULL, caption.level = NULL, frame = NULL, grid = NULL, col.width = 1, width = 0) {
  if (!is.null(frame)) {
    frame <- paste("frame=\"", switch(frame, topbot = "topbot", sides = "sides", all = "all", none = "none"), "\"", sep = "")
  }
  if (!is.null(grid)) {
    grid <- paste("grid=\"", switch(grid, all = "all", rows = "rows", cols = "cols", none = "none"), "\"", sep = "")
  }
  if (width != 0) {
    width <- paste("width=\"", width, "%\"", sep = "")
  } else {
    width <- ""
  }
  if (sum(col.width) > length(col.width)) {
    col.width <- paste("cols=\"", paste(col.width, collapse = ","), "\"", sep = "")
  } else {
    col.width <- ""
  }
  listarg <- c(frame, grid, col.width, width)
  listarg <- listarg[listarg != ""]
  if (length(listarg) != 0) {
    res <- paste("[", paste(listarg, collapse = ","), "]\n", sep = "")
  } else {
    res <- ""
  }
  if (!is.null(caption)) {
    if (is.null(caption.level)) {
      res <- paste(".", caption, "\n", res, sep = "")
    } else if (caption.level == ".") {
      res <- paste(".", caption, "\n", res, sep = "")
    } else if (is.numeric(caption.level) & caption.level > 0) {
      lev <- paste(rep("=", caption.level), collapse = "")
      res <- paste(lev, " ", caption, " ", "\n\n", res, sep = "")
    } else if (caption.level == "s") {
      res <- paste(beauty.asciidoc(caption, "s"), "\n\n", res, sep = "")
    }
    else if (caption.level == "e") {
      res <- paste(beauty.asciidoc(caption, "e"), "\n\n", res, sep = "")
    }
    else if (caption.level == "m") {
      res <- paste(beauty.asciidoc(caption, "m"), "\n\n", res, sep = "")
    }
    else if (caption.level == "none" | caption.level == "") {
      res <- paste(caption, "\n\n", res, sep = "")
    }
  }
  return(res)
}


##' escape.asciidoc
##'
##' @keywords internal
##' @param x x
escape.asciidoc <- function(x) {
  xx <- gsub("\\|", "\\\\|", x)
  xx
}


##' show.asciidoc.table
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
show.asciidoc.table <- function(x, include.rownames = FALSE, include.colnames = FALSE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = FALSE, footer = FALSE, align = NULL, col.width = 1, style = NULL, lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", ...) {

  x <- escape.asciidoc(tocharac(x, include.rownames, include.colnames, rownames, colnames, format, digits, decimal.mark, na.print))
  nrowx <- nrow(x)
  ncolx <- ncol(x)
  col.width <- rep(col.width, length.out = ncolx)
  
  line_separator <- FALSE
  vsep = matrix("|", nrow(x), ncol(x))

  if (!is.null(style)) {
    style <- expand(style, nrow(x), ncol(x))
  } else {
    style <- expand("", nrow(x), ncol(x))
  }
  if (is.logical(header) & header)
    header <- 1
  if (header > 0) {
    style[1:min(c(header, nrowx)),] <- "h"
  }
  if (is.logical(footer) & footer)
    footer <- 1
  if (footer > 0) {
    style[nrow(x):(nrow(x)+1-min(c(footer, nrowx))),] <- "h"
  }
  
  if (!is.null(align)) {
    align <- expand(align, nrow(x), ncolx)
  }
  if (!is.null(valign)) {
    valign <- expand(valign, nrow(x), ncolx)
  }
  before_vsep <- vsep.asciidoc(align, valign, style)

  if (include.rownames & include.colnames)
    before_vsep[1, 1] <- ""
  
  head <- header.asciidoc(caption, caption.level, frame, grid, col.width, width)
  results <- print.character.matrix(x, line_separator = line_separator, vsep = vsep, before_vsep = before_vsep, print = FALSE)
  if (include.rownames & include.colnames)
    results[1] <- sub("^ +", "", substr(results[1], 5, nchar(results[1])))
  
  # groups
  if (!is.null(lgroup)) {
    if (!is.list(lgroup))
      lgroup <- list(lgroup)
    n.lgroup <- groups(lgroup, n.lgroup, nrowx-include.colnames)[[2]]
  }
  if (!is.null(rgroup)) {
    if (!is.list(rgroup))
      rgroup <- list(rgroup)
    n.rgroup <- groups(rgroup, n.rgroup, nrowx-include.colnames)[[2]]
  }
  if (!is.null(tgroup)) {
    if (!is.list(tgroup))
      tgroup <- list(tgroup)
    n.tgroup <- groups(tgroup, n.tgroup, ncolx-include.rownames)[[2]]
  }
  if (!is.null(bgroup)) {
    if (!is.list(bgroup))
      bgroup <- list(bgroup)
    n.bgroup <- groups(bgroup, n.bgroup, ncolx-include.rownames)[[2]]
  }
  
  if (!is.null(lgroup)) {
    for (i in 1:length(lgroup)) {  
      pos.lgroup <- ngroups(lgroup[[i]], n.lgroup[[i]], n = nrowx)
      results[pos.lgroup[, 2]+include.colnames] <- paste(paste(paste(paste(".", pos.lgroup[, 3], "+", sep = ""), vsep.asciidoc(align = lalign, valign = lvalign, style = lstyle), sep = ""), pos.lgroup[, 1], sep = "| "), results[pos.lgroup[, 2]+include.colnames], sep = " ")
    }
  }
  
  if (!is.null(rgroup)) {
    for (i in 1:length(rgroup)) {
      pos.rgroup <- ngroups(rgroup[[i]], n.rgroup[[i]], n = nrowx)
      results[pos.rgroup[, 2]+include.colnames] <- paste(results[pos.rgroup[, 2]+include.colnames], paste(paste(paste(".", pos.rgroup[, 3], "+", sep = ""), vsep.asciidoc(align = ralign, valign = rvalign, style = rstyle), sep = ""), pos.rgroup[, 1], sep = "| "), sep = " ")
    }
  }

  if (!is.null(tgroup)) {
    for (i in 1:length(tgroup)) {
      pos.tgroup <- ngroups(tgroup[[i]], n.tgroup[[i]], n = ncolx)
      results <- c(paste(paste(paste(pos.tgroup[, 3], "+", sep = ""), vsep.asciidoc(align = talign, valign = tvalign, style = tstyle), sep = ""), pos.tgroup[, 1], sep = "| ", collapse = " "), results)
    }
  }
  
  if (!is.null(bgroup)) {
    for (i in 1:length(bgroup)) {
      pos.bgroup <- ngroups(bgroup[[i]], n.bgroup[[i]], n = ncolx)
      results <- c(results, paste(paste(paste(pos.bgroup[, 3], "+", sep = ""), vsep.asciidoc(align = balign, valign = bvalign, style = bstyle), sep = ""), pos.bgroup[, 1], sep = "| ", collapse = " "))
    }
  }
  
  topleftrow <- 0 + include.colnames + length(tgroup)
  topleftcol <- 0 + include.rownames + length(lgroup)
  topleft <- ""
  if (topleftrow > 0 & topleftcol > 0)
    topleft <- paste(topleftcol, ".", topleftrow, "+| ", sep = "")

  bottomleftrow <- 0 + length(bgroup)
  bottomleftcol <- 0 + include.rownames + length(lgroup)
  bottomleft <- ""
  if (bottomleftrow > 0 & bottomleftcol > 0)
    bottomleft <- paste(bottomleftcol, ".", bottomleftrow, "+| ", sep = "")

  toprightrow <- 0 + include.colnames + length(tgroup)
  toprightcol <- 0 + length(rgroup)
  topright <- ""
  if (toprightrow > 0 & toprightcol > 0)
    topright <- paste(" ", toprightcol, ".", toprightrow, "+| ", sep = "")

  bottomrightrow <- 0 + length(bgroup)
  bottomrightcol <- 0 + length(rgroup)
  bottomright <- ""
  if (bottomrightrow > 0 & bottomrightcol > 0)
    bottomright <- paste(" ", bottomrightcol, ".", bottomrightrow, "+| ", sep = "")

  results[1] <- paste(topleft, results[1], sep = "")
  results[1] <- paste(results[1], topright, sep = "")
  results[length(results)] <- paste(bottomleft, results[length(results)], sep = "")
  results[length(results)] <- paste(results[length(results)], bottomright, sep = "")
  
  topbot <- paste("|", paste(rep("=", max(nchar(results))), collapse = ""), sep = "")
  cat(head)
  cat(topbot, "\n")
  cat(results, sep = "\n")
  cat(topbot, "\n")
}


##' show.asciidoc.list
##'
##' @keywords internal
##' @param x x
##' @param caption caption
##' @param caption.level caption.level
##' @param list.type list.type
##' @param ... ...
show.asciidoc.list <- function(x, caption = NULL, caption.level = NULL, list.type = "bullet", ...) {
  
  if (list.type == "bullet") mark <- rep("\\*", length(x))
  if (list.type == "number") mark <- rep("\\.", length(x))
  if (list.type == "none")   mark <- rep("", length(x))
  if (list.type == "label") {
    if (is.null(names(x))) {
      namesx <- paste("[[", 1:length(x), "]]", sep = "")
    } else {
      namesx <- names(x)
    }
    mark <- paste(namesx, "::\n  ", sep = "")
  }
  
  charac.x <- vector("character", length(x))
  for (i in 1:length(x)) {
    if (is.null(x[[i]])) next
    tmp <- x[[i]]
    if (list.type == "label") tmp <- sub("^\t*", "", tmp)
    tmp <- sub("(^.*)", paste(mark[i], "\\1", sep = ""), gsub('\t|(*COMMIT)(*FAIL)', mark[i], tmp, perl = TRUE))
    if (list.type != "none")
      tmp <- sub(paste('(^', mark[i], '+)(.*)', sep = ""), '\\1 \\2', tmp)
    charac.x[i] <- tmp
  }
  cat(header.asciidoc(caption = caption, caption.level = caption.level))
  cat(charac.x, sep = "\n")
}
