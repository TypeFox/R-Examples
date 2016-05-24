wtd.hist <- function (x, breaks = "Sturges", freq = NULL, probability = !freq, 
                      include.lowest = TRUE, right = TRUE, density = NULL, angle = 45, 
                      col = NULL, border = NULL, main = paste("Histogram of", xname), 
                      xlim = range(breaks), ylim = NULL, xlab = xname, ylab, axes = TRUE, 
                      plot = TRUE, labels = FALSE, nclass = NULL, weight=NULL, ...){ #weight command was added
  ## PLEASE NOTE THAT THIS SOFTWARE IS A COPIED AND MODIFIED VERSION OF THE hist.default FUNCTION PROVIDED IN THE R GRAPHICS PACKAGE AS OF THE TIME THIS SOFTWARE WAS MODIFIED.  THAT SOFTWARE WAS COPYRIGHT OF R-CORE (2010) AND WAS MODIFIED UNDER THE TERMS OF GNU LICENSE 2.  ALL CHANGES TO THE ORIGINAL SOFTWARE ARE NOTATED IN CODE.
 # require(Hmisc) # Requirement for Hmisc was added
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if(is.null(weight)) # added
    weight <- rep(1, length(x)) # added
  xname <- paste(deparse(substitute(x), 500), collapse = "\n")
  n <- sum(weight[is.finite(x)]) # modified
  weight <- weight[is.finite(x)] # added
  x <- x[is.finite(x)]
  use.br <- !missing(breaks)
  if (use.br) {
    if (!missing(nclass)) 
      warning("'nclass' not used when 'breaks' is specified")
  }
  else if (!is.null(nclass) && length(nclass) == 1) 
    breaks <- nclass
  use.br <- use.br && (nB <- length(breaks)) > 1
  if (use.br) 
    breaks <- sort(breaks)
  else {
    if (!include.lowest) {
      include.lowest <- TRUE
      warning("'include.lowest' ignored as 'breaks' is not a vector")
    }
    if (is.character(breaks)) {
      breaks <- match.arg(tolower(breaks), c("sturges", 
                                             "fd", "freedman-diaconis", "scott"))
      breaks <- switch(breaks, sturges = nclass.Sturges(x), 
                       `freedman-diaconis` = , fd = nclass.FD(x), scott = nclass.scott(x), 
                       stop("unknown 'breaks' algorithm"))
    }
    else if (is.function(breaks)) {
      breaks <- breaks(x)
    }
    if (!is.numeric(breaks) || !is.finite(breaks) || breaks < 
        1) 
      stop("invalid number of 'breaks'")
    breaks <- pretty(range(x), n = breaks, min.n = 1)
    nB <- length(breaks)
    if (nB <= 1) 
      stop("hist.default: pretty() error, breaks=", format(breaks))
  }
  h <- diff(breaks)
  equidist <- !use.br || diff(range(h)) < 1e-07 * mean(h)
  if (!use.br && any(h <= 0)) 
    stop("'breaks' are not strictly increasing")
  freq1 <- freq
  if (is.null(freq)) {
    freq1 <- if (!missing(probability)) 
      !as.logical(probability)
    else equidist
  }
  else if (!missing(probability) && any(probability == freq)) 
    stop("'probability' is an alias for '!freq', however they differ.")
  diddle <- 1e-07 * stats::median(diff(breaks))
  fuzz <- if (right) 
    c(if (include.lowest) -diddle else diddle, rep.int(diddle, 
                                                       length(breaks) - 1))
  else c(rep.int(-diddle, length(breaks) - 1), if (include.lowest) diddle else -diddle)
  fuzzybreaks <- breaks + fuzz
  h <- diff(fuzzybreaks)
  storage.mode(x) <- "numeric"
  storage.mode(fuzzybreaks) <- "numeric"
  counts <- as.numeric(xtabs(weight~cut(x, fuzzybreaks))) # modified
  if (any(counts < 0)) 
    stop("negative 'counts'. Internal Error in C-code for \"bincount\"")
  if (sum(counts) < n-.01) 
    stop("some 'x' not counted; maybe 'breaks' do not span range of 'x'")
  dens <- counts/(n * diff(breaks))
  mids <- 0.5 * (breaks[-1L] + breaks[-nB])
  r <- structure(list(breaks = breaks, counts = counts, intensities = dens, 
                      density = dens, mids = mids, xname = xname, equidist = equidist), 
                 class = "histogram")
  if (plot) {
    plot(r, freq = freq1, col = col, border = border, angle = angle, 
         density = density, main = main, xlim = xlim, ylim = ylim, 
         xlab = xlab, ylab = ylab, axes = axes, labels = labels, 
         ...)
    invisible(r)
  }
  else {
    nf <- names(formals())
    nf <- nf[is.na(match(nf, c("x", "breaks", "nclass", "plot", 
                               "include.lowest", "weight", "right")))]
    missE <- lapply(nf, function(n) substitute(missing(.), 
                                               list(. = as.name(n))))
    not.miss <- !sapply(missE, eval, envir = environment())
    if (any(not.miss)) 
      warning(sprintf(ngettext(sum(not.miss), "argument %s is not made use of", 
                               "arguments %s are not made use of"), paste(sQuote(nf[not.miss]), 
                                                                          collapse = ", ")), domain = NA)
    r
  }
}
