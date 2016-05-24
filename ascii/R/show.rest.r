
##' header.rest
##'
##' @keywords internal
##' @param caption caption
##' @param caption.level caption.level
header.rest <- function(caption = NULL, caption.level = NULL) {
  niv <- c("=", "-", "~", "^", "+")
  ncharcap <- nchar(caption)
  if (is.null(caption.level))
    caption.level <- ""
  res <- ""
  if (!is.null(caption)) {
    if (is.numeric(caption.level) & caption.level > 0) {
      res <- paste(caption, paste(paste(rep(niv[caption.level], ncharcap), collapse = ""), "\n", sep = ""), sep = "\n")
    } else if (is.character(caption.level) & caption.level %in% c("s", "e", "m")) {
      if (caption.level == "s")
        res <- paste(beauty.rest(caption, "s"), "\n\n", sep = "")
      else if (caption.level == "e")
        res <- paste(beauty.rest(caption, "e"), "\n\n", sep = "")
      else if (caption.level == "m")
        res <- paste(beauty.rest(caption, "m"), "\n\n", sep = "")
    } else if (is.character(caption.level) & caption.level != "" & caption.level != "none") {
      res <- paste(caption, paste(paste(rep(caption.level, ncharcap), collapse = ""), "\n", sep = ""), sep = "\n")
    } else if (caption.level == "" | caption.level == "none") {
      res <- paste(caption, "\n\n", sep = "")
    }
  }
  return(res)
}
  

##' beauty.rest
##'
##' @keywords internal
##' @param x x
##' @param beauti beauti
beauty.rest <- function(x, beauti = c("e", "m", "s")) {
  x[is.na(x)] <- "NA"
  if (beauti == "s") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("\\*\\*.*\\*\\*", x)+1)/2) # bold seulement si != de "" et si pas de bold
    if (length(x[!y]) != 0) x[!y] <- sub("(^ *)([:alpha]*)", "\\1\\*\\*\\2", sub("([:alpha:]*)( *$)", "\\1\\*\\*\\2", x[!y]))
    if (length(x[y]) != 0) x[y] <- sub("(^ *$)", "\\1 ", x[y]) # rajouter suffisamment d'espaces lorsque la case est vide pour l'alignement globale
  }
  if (beauti == "e") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("\\*.*\\*", x)+1)/2) # it seulement si != de "" et si pas de it
    if (length(x[!y]) != 0) x[!y] <-sub("(^ *)([:alpha]*)", "\\1\\*\\2", sub("([:alpha:]*)( *$)", "\\1\\*\\2", x[!y]))
    if (length(x[y]) != 0) x[y] <- sub("(^ *$)", "\\1 ", x[y]) # rajouter suffisamment d'espaces lorsque la case est vide pour l'alignement globale
  }
  if (beauti == "m") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("``.*``", x)+1)/2) # it seulement si != de "" et si pas de mono
    if (length(x[!y]) != 0) x[!y] <-sub("(^ *)([:alpha]*)", "\\1``\\2", sub("([:alpha:]*)( *$)", "\\1``\\2", x[!y]))
    if (length(x[y]) != 0) x[y] <- sub("(^ *$)", "\\1 ", x[y]) # rajouter suffisamment d'espaces lorsque la case est vide pour l'alignement globale
  }
  return(x)
}


##' escape.rest
##'
##' @keywords internal
##' @param x x
escape.rest <- function(x) {
  xx <- gsub("\\|", "\\\\|", x)
  xx
}

##' show.rest.table
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
show.rest.table <- function(x, include.rownames = FALSE, include.colnames = FALSE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = FALSE, footer = FALSE, align = NULL, col.width = 1, style = NULL, lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", ...) {

  x <- escape.rest(tocharac(x, include.rownames, include.colnames, rownames, colnames, format, digits, decimal.mark, na.print))
  nrowx <- nrow(x)
  ncolx <- ncol(x)

  if (!is.null(style)) {
    style <- expand(style, nrowx, ncolx)
    style[!(style %in% c("s", "e", "m"))] <- ""
    style[style == "s"] <- "**"
    style[style == "e"] <- "*"
    style[style == "m"] <- "``"
  } else {
    style <- ""
    style <- expand(style, nrowx, ncolx)
  }
  if (include.rownames & include.colnames) {
    style[1, 1] <- ""
  }
  x <- paste.matrix(style, x, style, sep = "", transpose.vector = TRUE)

  if (tstyle == "h")
    tstyle <- "s"
  if (bstyle == "h")
    bstyle <- "s"
  if (rstyle == "h")
    rstyle <- "s"
  if (lstyle == "h")
    lstyle <- "s"
  
  if (!is.null(lgroup)) {
    if (!is.list(lgroup))
      lgroup <- list(lgroup)
    n.lgroup <- groups(lgroup, n.lgroup, nrowx-include.colnames)[[2]]
    lgroup <- lapply(lgroup, function(x) beauty.rest(x, lstyle))
  }
  if (!is.null(rgroup)) {
    if (!is.list(rgroup))
      rgroup <- list(rgroup)
    n.rgroup <- groups(rgroup, n.rgroup, nrowx-include.colnames)[[2]]
    rgroup <- lapply(rgroup, function(x) beauty.rest(x, rstyle))
  }
  if (!is.null(tgroup)) {
    if (!is.list(tgroup))
      tgroup <- list(tgroup)
    n.tgroup <- groups(tgroup, n.tgroup, ncolx-include.rownames)[[2]]
    tgroup <- lapply(tgroup, function(x) beauty.rest(x, tstyle))
  }
  if (!is.null(bgroup)) {
    if (!is.list(bgroup))
      bgroup <- list(bgroup)
    n.bgroup <- groups(bgroup, n.bgroup, ncolx-include.rownames)[[2]]
    bgroup <- lapply(bgroup, function(x) beauty.rest(x, bstyle))
  }
  if (!is.null(tgroup)) {
    tgroup.vsep <- NULL
    for (i in 1:length(tgroup)) {
      line.tgroup <- unlist(interleave(as.list(tgroup[[i]]), lapply(n.tgroup[[i]], function(x) rep("", x-1))))
      if (include.rownames) {
        line.tgroup <- c("", line.tgroup)
      }
      x <- rbind(line.tgroup, x)
      tvsep <- c(unlist(interleave(as.list(rep("|", length(tgroup[[i]]))), lapply(n.tgroup[[i]], function(x) rep(" ", x-1)))), "|")
      tgroup.vsep <- rbind(tvsep, tgroup.vsep)
    }
  }
  if (!is.null(bgroup)) {
    bgroup.vsep <- NULL
    for (i in 1:length(bgroup)) {
      line.bgroup <- unlist(interleave(as.list(bgroup[[i]]), lapply(n.bgroup[[i]], function(x) rep("", x-1))))
      if (include.rownames) {
        line.bgroup <- c("", line.bgroup)
      }
      x <- rbind(x, line.bgroup)
      bvsep <- c(unlist(interleave(as.list(rep("|", length(bgroup[[i]]))), lapply(n.bgroup[[i]], function(x) rep(" ", x-1)))), "|")
      bgroup.vsep <- rbind(bgroup.vsep, bvsep)
    }
  }
  if (!is.null(lgroup)) {
    lgroup.hsep <- NULL
    for (i in 1:length(lgroup)) {
      line.lgroup <- unlist(interleave(as.list(lgroup[[i]]), lapply(n.lgroup[[i]], function(x) rep("", x-1))))
      line.lgroup <- c(rep("", include.colnames + length(tgroup)), line.lgroup)
      line.lgroup <- c(line.lgroup, rep("", length(bgroup)))
      x <- cbind(line.lgroup, x)
      lhsep <- c(unlist(interleave(as.list(rep("-", length(lgroup[[i]]))), lapply(n.lgroup[[i]], function(x) rep(" ", x-1)))), "-")
      lgroup.hsep <- cbind(lhsep, lgroup.hsep)
    }
  }
  if (!is.null(rgroup)) {
    rgroup.hsep <- NULL
    for (i in 1:length(rgroup)) {
      line.rgroup <- unlist(interleave(as.list(rgroup[[i]]), lapply(n.rgroup[[i]], function(x) rep("", x-1))))
      line.rgroup <- c(rep("", include.colnames + length(tgroup)), line.rgroup)
      line.rgroup <- c(line.rgroup, rep("", length(bgroup)))
      x <- cbind(x, line.rgroup)
      rhsep <- c(unlist(interleave(as.list(rep("-", length(rgroup[[i]]))), lapply(n.rgroup[[i]], function(x) rep(" ", x-1)))), "-")
      rgroup.hsep <- cbind(rgroup.hsep, rhsep)
    }
  }

  vsep <- expand("|", nrowx+length(tgroup)+length(bgroup), ncolx+1+length(lgroup) + length(rgroup))
  if (!is.null(tgroup)) {
    vsep[1:length(tgroup), (length(lgroup)+include.rownames+1):(ncol(vsep)-length(rgroup))] <- tgroup.vsep
  }
  if (!is.null(bgroup)) {
    vsep[(length(tgroup)+nrowx+1):(length(tgroup)+nrowx+length(bgroup)), (length(lgroup)+include.rownames+1):(ncol(vsep)-length(rgroup))] <- bgroup.vsep
  }
  if ((length(lgroup)+include.rownames >= 1) & (length(tgroup)+include.colnames >= 1)) {
    topleft <- matrix(" ", length(tgroup)+include.colnames, length(lgroup)+include.rownames)
    topleft[, 1] <- "|"
    vsep[1:nrow(topleft), 1:ncol(topleft)] <- topleft
  }
  if ((length(lgroup)+include.rownames >=1) & (length(bgroup) >= 1)) {
    bottomleft <- matrix(" ", length(bgroup), length(lgroup)+include.rownames)
    bottomleft[, 1] <- "|"
    vsep[(nrow(vsep)-length(bgroup)+1):(nrow(vsep)), 1:ncol(bottomleft)] <- bottomleft
  }
  if ((length(rgroup) >= 1) & (include.colnames+length(tgroup) >= 1)) {
    topright <- matrix(" ", length(tgroup)+include.colnames, length(rgroup))
    topright[, ncol(topright)] <- "|"
    vsep[1:nrow(topright), (ncol(vsep)-length(rgroup)+1):ncol(vsep)] <- topright
  }
  if ((length(rgroup) >=1 ) & (length(bgroup) >= 1)) {
    bottomright <- matrix(" ", length(bgroup), length(rgroup))
    bottomright[, ncol(bottomright)] <- "|"
    vsep[(nrow(vsep)-nrow(bottomright)+1):nrow(vsep), (ncol(vsep)-length(rgroup)+1):ncol(vsep)] <- bottomright
  }

  hsep <- expand("-", nrowx+1+length(tgroup)+length(bgroup), ncolx+length(lgroup)+length(rgroup))
  if (!is.null(lgroup)) {
    hsep[(length(tgroup)+include.colnames+1):(nrow(hsep)-length(bgroup)), 1:length(lgroup)] <- lgroup.hsep
  }
  if (!is.null(rgroup)) {
    hsep[(length(tgroup)+include.colnames+1):(nrow(hsep)-length(bgroup)), (ncol(hsep)-length(rgroup)+1):(ncol(hsep))] <- rgroup.hsep
  }
  if ((length(lgroup)+include.rownames >= 1) & (length(tgroup)+include.colnames >= 1)) {
    topleft <- matrix(" ", length(tgroup)+include.colnames, length(lgroup)+include.rownames)
    topleft[1, ] <- "-"
    hsep[1:nrow(topleft), 1:ncol(topleft)] <- topleft
  }
  if ((length(lgroup)+include.rownames >= 1)&(length(bgroup) >= 1)) {
    bottomleft <- matrix(" ", length(bgroup), length(lgroup)+include.rownames)
    bottomleft[nrow(bottomleft), ] <- "-"
    hsep[(nrow(hsep)-length(bgroup)+1):(nrow(hsep)), 1:ncol(bottomleft)] <- bottomleft
  }
  if ((length(rgroup) >= 1) & (include.colnames+length(tgroup) >= 1)) {
    topright <- matrix(" ", length(tgroup)+include.colnames, length(rgroup))
    topright[1, ] <- "-"
    hsep[1:nrow(topright), (ncol(hsep)-length(rgroup)+1):ncol(hsep)] <- topright
  }
  if ((length(rgroup) >= 1) & (length(bgroup) >= 1)) {
    bottomright <- matrix(" ", length(bgroup), length(rgroup))
    bottomright[nrow(bottomright), ] <- "-"
    hsep[(nrow(hsep)-nrow(bottomright)+1):nrow(hsep), (ncol(hsep)-length(rgroup)+1):ncol(hsep)] <- bottomright
  }

  csep <- matrix("+", nrowx+1+length(tgroup)+length(bgroup), ncolx+length(lgroup)+length(rgroup)+1)
  if ((length(lgroup)+include.rownames >= 1) & (length(tgroup)+include.colnames >= 1)) {
    topleft <- matrix(" ", length(tgroup)+include.colnames, length(lgroup)+include.rownames)
    topleft[1, ] <- "+"
    topleft[, 1] <- "+"
    csep[1:nrow(topleft), 1:ncol(topleft)] <- topleft
  }
  if ((length(lgroup)+include.rownames >=1) & (length(bgroup) >= 1)) {
    bottomleft <- matrix(" ", length(bgroup), length(lgroup)+include.rownames)
    bottomleft[nrow(bottomleft), ] <- "+"
    bottomleft[, 1] <- "+"
    csep[(nrow(csep)-length(bgroup)+1):(nrow(csep)), 1:ncol(bottomleft)] <- bottomleft
  }
  if ((length(rgroup) >= 1) & (include.colnames+length(tgroup) >= 1)) {
    topright <- matrix(" ", length(tgroup)+include.colnames, length(rgroup))
    topright[1, ] <- "+"
    topright[, ncol(topright)] <- "+"
    csep[1:nrow(topright), (ncol(csep)-length(rgroup)+1):ncol(csep)] <- topright
  }
  if ((length(rgroup) >= 1) & (length(bgroup) >= 1)) {
    bottomright <- matrix(" ", length(bgroup), length(rgroup))
    bottomright[nrow(bottomright), ] <- "+"
    bottomright[, ncol(bottomright)] <- "+"
    csep[(nrow(hsep)-nrow(bottomright)+1):nrow(csep), (ncol(csep)-length(rgroup)+1):ncol(csep)] <- bottomright
  }

  if (header) {
    header <- min(c(header+length(tgroup), nrowx+length(tgroup)))+1
    hsep[header,] <- gsub("-", "=", hsep[header,])
  }
  results <- print.character.matrix(x, vsep = vsep, csep = csep, hsep = hsep, print = FALSE)
  cat(header.rest(caption = caption, caption.level = caption.level), sep = "\n")
  cat(results, sep = "\n")
}


##' show.rest.list
##'
##' @keywords internal
##' @param x x
##' @param caption caption
##' @param caption.level caption.level
##' @param list.type list.type
##' @param ... ...
show.rest.list <- function(x, caption = NULL, caption.level = NULL, list.type = "bullet", ...) {
  if (list.type == "bullet") mark <- rep("*", length(x))
  if (list.type == "number") mark <- rep("#.", length(x))
  if (list.type == "none")  mark <- rep("", length(x))
  if (list.type == "label") {
    if (is.null(names(x))) {
      namesx <- paste("[[", 1:length(x), "]]", sep = "")
    } else {
      namesx <- names(x)
    }
    mark <- paste(namesx, "\n ", sep = "")
  }
  y <- gsub("(^\t*)(.*)", "\\1", x)
  z <- NULL
  for (i in 2:length(y))
    z <- c(z, ifelse(y[i] != y[i-1], i-1, NA))
  
  cat(header.rest(caption = caption, caption.level = caption.level), sep = "\n")
  
  for (i in 1:length(x)) {
    tmp <- x[[i]]
    if (list.type == "label") tmp <- sub("^\t*", "", tmp)
    tmp <- gsub('\t|(*COMMIT)(*FAIL)', "  ", tmp, perl = TRUE)
    tmp <- sub("(^ *)", paste("\\1", mark[i], " ", sep = ""), tmp)
    cat(tmp, "\n")
    if (i %in% z)
      cat("\n")
  }
}

