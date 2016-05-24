contrasting.rainbow <- function (x, ...) {
  x <- as.integer(x)[1]
  stopifnot(x >= 1)
  ## make pretty colour cycle, by...
  ## ...finding non-unit divisors of ler...
  xnud <- setdiff(which(x %% 1:x == 0), 1)
  ## ...and thus largest coprime of x that is <=x/2...
  xcop <- max(which(sapply(1:(x/2), function(div) all(div %% xnud != 0))))
  ## ...and thus constrasting traversal of x
  xtrav <- (1:x*xcop) %% x + 1
  ## fallback; no coprime exists so use random shuffle
  if (!all(1:x %in% xtrav)) xtrav <- sample(x)
  return(rainbow(x, ...)[xtrav])
}


## plotpos.by.chr <- function (chr, pos, gap = 5e7, chrset = c(1:22, "X", "Y", "M")) {
##   stopifnot(length(pos) == length(chr))
##   chr <- toupper(as.character(chr))
##   chrset <- toupper(as.character(chrset))
##   plotpos <- rep(NA, length(pos))
##   ox <- 0
##   for (cx in chrset) {
##     this <- which(chr == cx | chr == paste("CHR", cx, sep = ""))
##     plotpos[this] <- pos[this] - min(pos[this], na.rm = TRUE) + ox + gap
##     ox <- max(plotpos[this], na.rm = TRUE)
##   }
##   return(plotpos)
## }

## plotcol.by.chr <- function (chr, cols = NULL, chrset = c(1:22, "X", "Y", "M")) {
##   chrset <- toupper(as.character(chrset))
##   chrnum <- match(sub("^CHR", "", toupper(as.character(chr))), chrset)
##   if (is.null(cols)) {
##     cdiv <- which.min(abs(length(chrset) %% (1:length(chrset)) - length(chrset)/3))
##     cols <- rainbow(length(chrset))[(1:length(chrset) * cdiv) %% length(chrset) + 1]
##   }
##   return(rep(cols, length.out = length(chrset))[chrnum])
## }

## mid.by.chr <- function (chr, plotpos) {
##   minpos <- rep(NA, 22)
##   medpt <- rep(NA, 22)
##   maxpos <- rep(NA, 22)
##   for (cx in 1:22) {
##     tmp <- quantile(plotpos[which(chr == cx | chr == paste("chr", cx, sep = ""))], probs = c(0, 0.5, 1), na.rm = TRUE)
##     minpos[cx] <- tmp[1]
##     medpt[cx] <- tmp[2]
##     maxpos[cx] <- tmp[3]
##   }
##   return (list(minpos = minpos, midpt = (minpos+maxpos)/2, medpt = medpt, maxpos = maxpos))
## }
