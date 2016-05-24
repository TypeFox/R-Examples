## main function
nops_scan <- function(
  images = dir(pattern = "\\.PNG$|\\.png$|\\.PDF|\\.pdf$"),
  file = NULL, dir = ".",
  verbose = TRUE, rotate = FALSE, cores = NULL, n = NULL,
  density = 300,
  size = 0.029, threshold = c(0.04, 0.42), minrot = 0.002,
  string = FALSE)
{
  ## required packages
  stopifnot(requireNamespace("png"))
  if(!is.null(cores)) {
    if(!requireNamespace("parallel")) cores <- NULL
  }

  ## check whether images exist
  if(!all(im <- file.exists(images))) {
    warning(paste("The following images cannot be found:", paste(images[!im], collapse = ", ")))
    images <- images[im]
  }
  if(length(images) < 0L) {
    stop("No images found.")
  }

  ## as an alternative to the R code for reading printed digits
  ## one could call the command line tool "tesseract" but this
  ## turned out to be less reliable for this special case and is
  ## hence never used here:
  tesseract <- FALSE

  ## directory handling
  dir <- file_path_as_absolute(dir)
  owd <- getwd()
  dir.create(tdir <- tempfile())
  on.exit(unlink(tdir))

  ## convert PDF to PNG (if necessary)
  pdfs <- grepl("\\.pdf$", tolower(images))
  if(any(pdfs)) {
    pngs <- pdfs2pngs(images[pdfs], density = density, cores = cores, dir = tdir,
      verbose = verbose, rotate = rotate, prefix = if(string) "T" else "S")
    images <- if(length(images) > sum(pdfs)) c(images[!pdfs], pngs) else pngs
  }
  
  file.copy(images, file.path(tdir, images <- basename(images)))
  setwd(tdir)
  on.exit(setwd(owd), add = TRUE)

  ## read nops images
  if(verbose) cat("Reading PNG files:\n")

  read_nops_all <- function(file)
  {
    err <- paste(file, "ERROR")

    if(verbose) {
      if(is.null(cores)) {
        cat(paste(file, ": Trimming PNG", sep = ""))
      } else {
        cat(paste("Reading ", file, ".\n", sep = ""))
      }
    }
    ss <- try(trim_nops_scan(file, verbose = verbose & is.null(cores), minrot = minrot))
    if(inherits(ss, "try-error")) {
      if(verbose) cat(", ERROR")
      return(err)
    }
    
    if(verbose & is.null(cores)) cat(", extracting information")
    ss <- if(!string) {
      try(paste(
        file,
        read_nops_digits(ss, "id", tesseract = tesseract),
        read_nops_digits(ss, "scrambling", tesseract = tesseract),
        ssty <- read_nops_digits(ss, "type", tesseract = tesseract),
        read_nops_backup(ss, threshold = threshold, size = size),
        read_nops_registration(ss, threshold = threshold, size = size),
        read_nops_answers(ss, threshold = threshold, size = size, n = if(is.null(n)) as.numeric(ssty) else n)
      ))
    } else {
      try(paste(
        file,
        read_nops_digits(ss, "id", tesseract = tesseract, adjust = TRUE),
	read_nops_digits(ss, "type", tesseract = tesseract, adjust = TRUE),
        substr(read_nops_answers(ss, threshold = threshold, size = size, n = 3L, adjust = TRUE), 1, 17)
      ))
    }
    
    if(inherits(ss, "try-error")) {
      if(verbose) cat(", ERROR")
      return(err)
    }
  
    if(verbose & is.null(cores)) cat(", done.\n")
    return(ss)
  }
  read_nops <- function(x) as.vector(sapply(x, read_nops_all))
  
  rval <- if(!is.null(cores)) {
    xi <- split(images, ceiling(seq_along(images) / (length(images) / cores)))
    unlist(parallel::mclapply(seq_along(xi), function(j) { read_nops(xi[[j]]) }, mc.cores = cores))
  } else {
    read_nops(images)
  }

  ## return output
  if(!identical(file, FALSE)) {
    if(verbose) cat("\nCreating ZIP file:\n")
    if(is.null(file) || !is.character(file)) file <- paste(if(string) "nops_string_scan" else "nops_scan",
      format(Sys.time(), "%Y%m%d%H%M%S"), sep = "_")
    if(substr(tolower(file), nchar(file) - 3L, nchar(file)) != ".zip") file <- paste(file, "zip", sep = ".")
    writeLines(rval, file.path(tdir, if(string) "Daten2.txt" else "Daten.txt"))
    zip(zipfile = file, files = list.files(tdir))
    file.copy(file, file.path(dir, file))
    invisible(rval)
  } else {
    return(rval)
  }
}


## auxiliary functions (not to be exported)

## conversion of PDF to PNG images
pdfs2pngs <- function(x, density = 300, dir = NULL, cores = NULL, verbose = TRUE, rotate = FALSE, prefix = "S")
{
  ## copy to temporary directory
  dir.create(tdir <- tempfile())
  if(!is.null(dir)) {
    dir <- file_path_as_absolute(dir)
    on.exit(unlink(tdir))
  } else {
    dir <- dirname(x)
  }
  owd <- getwd()  
  file.copy(x, file.path(tdir, x <- basename(x)))
  setwd(tdir)

  shsystem <- function(cmd) {
    sh <- Sys.getenv("COMSPEC")
    if(sh != "") sh <- paste(shQuote(sh), "/c")
    system(paste(sh, cmd))
  }

  ## if necessary: merge PDFs, otherwise rename only
  if(length(x) > 1L) {
    if(verbose) cat("Merging PDF files")
    shsystem(sprintf("pdftk %s cat output _NOPS_.pdf", paste(x, collapse = " ")))
    file.remove(x)
    if(verbose) cat(", done.\n")
  } else {
    file.rename(x, "_NOPS_.pdf")
  }

  ## if requested: rotate PDFs
  if(rotate) {
    if(verbose) cat("Rotating PDF files")
    shsystem("pdftk _NOPS_.pdf rotate 1-enddown output _NOPS_2_.pdf")
    file.remove("_NOPS_.pdf")
    file.rename("_NOPS_2_.pdf", "_NOPS_.pdf")
    if(verbose) cat(", done.\n")
  }

  ## burst PDF into individual pages
  if(verbose) cat("Splitting PDF files")
  shsystem(paste0("pdftk _NOPS_.pdf burst output ", prefix, "%07d.pdf"))
  file.remove(c("_NOPS_.pdf", "doc_data.txt"))
  if(verbose) cat(", done.\n")
  x <- dir(pattern = "\\.pdf$")

  ## actual conversion function
  pdf2png <- function(pdfs) {
    ## shell command on Windows
    for(i in pdfs) {
      if(verbose) cat(paste(i, ": Converting PDF to PNG.\n", sep = ""))
      cmd <- paste("convert -density", density, i, gsub(".pdf", ".PNG", i, fixed = TRUE))
      shsystem(cmd)
    }
  }

  if(!is.null(cores)) {
    xi <- split(x, ceiling(seq_along(x) / (length(x) / cores)))
    parallel::mclapply(1:cores, function(j) { pdf2png(xi[[j]]) }, mc.cores = cores)
  } else {
    pdf2png(x)
  }
  pngs <- gsub(".pdf", ".PNG", x, fixed = TRUE)  
  file.copy(pngs, pngs <- file.path(dir, pngs))
  setwd(owd)

  if(verbose) cat("\n")
  return(pngs)
}

## select sub-image from a pixel matrix
subimage <- function(x, center, prop = 0.01) {
  prop <- rep(prop, length.out = 2L)
  if(center[1L] < 1) center[1L] <- round(center[1L] * nrow(x))
  if(center[2L] < 1) center[2L] <- round(center[2L] * ncol(x))
  topleft  <- center - round(nrow(x) * prop/2)
  botright <- center + round(nrow(x) * prop/2)
  x[max(1L, topleft[1L]):min(nrow(x), botright[1L]), max(1L, topleft[2L]):min(ncol(x), botright[2L])]
}

"subimage<-" <- function(x, center, prop = 0.01, value) {
  prop <- rep(prop, length.out = 2L)
  if(center[1L] < 1) center[1L] <- round(center[1L] * nrow(x))
  if(center[2L] < 1) center[2L] <- round(center[2L] * ncol(x))
  topleft  <- center - round(nrow(x) * prop/2)
  botright <- center + round(nrow(x) * prop/2)
  x[max(1L, topleft[1L]):min(nrow(x), botright[1L]), max(1L, topleft[2L]):min(ncol(x), botright[2L])] <- value
  x
}

## shave (almost) white margins of a pixel matrix
shave <- function(x, zap = 0.07) {
  ix <- rowMeans(x) > zap
  if(any(ix)) {
    ix[min(which(ix)):max(which(ix))] <- TRUE
  } else {
    rep(TRUE, length(ix))
  }
  x <- x[ix, ]

  ix <- colMeans(x) > zap
  if(any(ix)) {
    ix[min(which(ix)):max(which(ix))] <- TRUE
  } else {
    rep(TRUE, length(ix))
  }
  x[, ix]
}

## shave box (and white margins) of a pixel matrix
shave_box <- function(x, border = 0.1, clip = TRUE)
{
  rm <- range(which(rowMeans(x) > 0.38))
  cm <- range(which(colMeans(x) > 0.38))
  x <- x[rm[1L]:rm[2L], cm[1L]:cm[2L]]
  n <- min(dim(x) * border)
  x <- x[n:(nrow(x) - n), n:(ncol(x) - n)]
  if(clip) shave(x) else x
}

## determine whether a pixel matrix has a check mark
has_mark <- function(x, threshold = c(0.04, 0.42), fuzzy = FALSE)
{
  rm <- which(rowMeans(x) > 0.38)
  cm <- which(colMeans(x) > 0.38)
  if(length(rm) < 2L || length(cm) < 2L || diff(range(rm)) < 5L || diff(range(cm)) < 5L) return(0L)
  rm <- range(rm)
  cm <- range(cm)
  x <- subimage(x[rm[1L]:rm[2L], cm[1L]:cm[2L]], c(0.5, 0.5), 0.75)
  if(mean(x) < threshold[1L]) return(0L)
  if(mean(x) < threshold[2L]) {
    if(fuzzy) return(mean(x)) else return(1L)
  } else {
    edges <- c(
      mean(subimage(x, c(0.5, 0.05), 0.1)),
      mean(subimage(x, c(0.5, 0.95), 0.1)),
      mean(subimage(x, c(0.05, 0.5), 0.1)),
      mean(subimage(x, c(0.95, 0.5), 0.1))
    )
    if(sort(edges)[2] <= 0.1) {
      if(fuzzy) return(mean(x)) else return(1L)
    } else {
      return(0L)
    }
  }
}

## read scanned PNG image into b/w pixel matrix and trim margins
trim_nops_scan <- function(x, verbose = FALSE, minrot = 0.002)
{
  ## read gray levels
  if(is.character(x)) {
    file <- x
    x <- png::readPNG(x)
  } else {
    file <- NULL
  }
  if(length(dim(x)) > 2L) {
    x <- pmin(x[, , 1L], x[, , 2L], x[, , 3L])
  }
  x <- matrix(as.integer(x < 0.7), nrow = nrow(x), ncol = ncol(x))
  d <- dim(x)

  ## force margins to be white
  x[, c(1L:round(0.02 * ncol(x)), round(0.98 * ncol(x)):ncol(x))] <- 0L
  x[c(1L:round(0.02 * nrow(x)), round(0.98 * nrow(x)):nrow(x)), ] <- 0L

  rot <- NULL
  while(is.null(rot) || (abs(rot) > minrot & abs(rot) < 0.05)) {

  ## rotate (if necessary)  
  if(!is.null(rot)) {
    rot <- 0.71 * rot ## try to avoid over-rotation
    if(verbose) cat(", rotating PNG")
    proj <- matrix(c(cos(rot), -sin(rot), sin(rot), cos(rot)), ncol = 2L)
    xcoord <- t(which(x > 0.5, arr.ind = TRUE))
    xcoord <- xcoord/d - 0.5
    xcoord <- proj %*% xcoord
    xcoord <- t(round(d * (xcoord + 0.5)))
    x[] <- 0L
    x[xcoord] <- 1L
  }

  ## find bottom markings
  xbl <- x[seq(round(0.93 * d[1L]), d[1L]), seq(1, round(0.15 * d[2L]))]
  xbr <- x[seq(round(0.93 * d[1L]), d[1L]), seq(round(0.85 * d[2L]), d[2L])]

  rb <- 0.93
  while(mean(xbl) < 0.0015 | mean(xbr) < 0.0015) {
    rb <- rb - 0.01
    xbl <- x[seq(round(rb * d[1L]), d[1L]), seq(1, round(0.15 * d[2L]))]
    xbr <- x[seq(round(rb * d[1L]), d[1L]), seq(round(0.85 * d[2L]), d[2L])]  
  }

  get_mark <- function(x, type = c("row", "col"), zap = 0.35)
  {
    x[rowMeans(x) >= zap,] <- 0
    x[,colMeans(x) >= zap] <- 0

    type <- match.arg(type)
    x <- if(type == "row") rowMeans(x) else colMeans(x)
    which(x > mean(range(x)))
  }
  get_mean <- function(x, maxdist = 10) {
    mean(x[abs(x - median(x)) < maxdist])
  }
  
  rbl <- get_mark(xbl, "row")
  rbr <- get_mark(xbr, "row")
  rb <- round(get_mean(unique(c(rbl, rbr))))
  rb <- as.vector(d[1L] - (nrow(xbl) - rb))

  cl <- round(get_mean(get_mark(xbl, "col")))
  cr <- round(get_mean(get_mark(xbr, "col")))
  cl <- as.vector(cl)
  cr <- as.vector(d[2L] - (ncol(xbr) - cr))

  ## rotation angle
  rot <- asin((get_mean(rbl) - get_mean(rbr)) / (cr - cl))
  }
  if(abs(rot) > 0.05) stop("image is too skewed, cannot be rotated")

  ## find top markings
  ctl <- round(cl + (cr - cl) * (30 - 20) / (190 - 20))
  ctr <- round(cl + (cr - cl) * (160 - 20) / (190 - 20))
  xtl <- x[seq(1L, round(0.13 * d[1L])), seq(1, round(1.15 * ctl))]
  xtr <- x[seq(1L, round(0.13 * d[1L])), seq(0.9 * ctr, d[2L])]
  xtl[, seq(1, 0.33 * ncol(xtr))] <- 0
  xtr[, seq(0.4 * ncol(xtr), ncol(xtr))] <- 0

  rtl <- get_mark(xtl, "row")
  rtr <- get_mark(xtr, "row") ## may be affected by text close to the mark, hence not used
  rt <- as.vector(round(get_mean(unique(rtl))))
  
  if(abs((rb - rt) / (cr - cl) - (270 - 13) / (190 - 20)) > 0.02)
    warning("PNG does not seem to be correctly scaled")

  ## extract subimage within markings
  x[rt:rb, cl:cr]
}

## set up row x col regressors with gray values of pixels
digit_regressors <- function(x, nrow = 7, ncol = 5)
{
  d <- dim(x)
  rw <- round(d[1L] * (0L:nrow/nrow))
  cl <- round(d[2L] * (0L:ncol/ncol))
  ix <- as.matrix(expand.grid(1:nrow, 1:ncol))
  rw1 <- rw[ix[, 1L]] + 1L
  rw2 <- rw[ix[, 1L] + 1L]
  cl1 <- cl[ix[, 2L]] + 1L
  cl2 <- cl[ix[, 2L] + 1L]
  rval <- sapply(1:nrow(ix), function(i) round(mean(x[rw1[i]:rw2[i], cl1[i]:cl2[i]]), digits = 4L))
  rval <- as.data.frame(t(rval))
  names(rval) <- paste("x", 1:ncol(rval), sep = "")
  rval$width <- round(d[2L]/d[1L], digits = 4L)
  return(rval)
}

## classify digits 
read_nops_digits <- function(x, type = c("type", "id", "scrambling"), tesseract = FALSE, adjust = FALSE)
{
  ## adjustment for coordinates (e.g. for reading 2nd string page)
  if(identical(adjust, TRUE)) adjust <- c(0.2065, 0)
  if(identical(adjust, FALSE)) adjust <- c(0, 0)
  
  ## extract image of numbers
  type <- match.arg(type)
  z <- switch(type,
    "type" = shave_box(subimage(x, c(0.3925 - adjust[1L], 0.074 - adjust[2L]), c(0.035, 0.078))),
    "id" = shave_box(subimage(x, c(0.3925 - adjust[1L], 0.275 - adjust[2L]), c(0.035, 0.19))),
    "scrambling" = {
      y <- shave_box(subimage(x, c(0.337 - adjust[1L], 0.545 - adjust[2L]), c(0.035, 0.078)), clip = FALSE)
      y[round(0.7 * nrow(y)):nrow(y), round(0.43 * ncol(y)):round(0.57 * ncol(y))] <- 0
      shave(y)
    })
  n <- switch(type,
    "type" = 3L,
    "id" = 11L,
    "scrambling" = 2L)

  ## split
  le <- NULL
  thresh <- 0
  while(length(le) != (2L * n - 1L) & thresh < 0.2) {
    thresh <- thresh + 0.01
    le <- rle(colMeans(z) < thresh)$lengths
  }
  if(length(le) != (2L * n - 1L)) return(paste(rep("X", n), collapse = ""))
  le <- cumsum(c(0L, le))
  d <- lapply(1L:(length(le)/2L), function(i) z[, (le[2 * i - 1L] + 1L):(le[2 * i]), drop = FALSE])
  
  ## transform to regressors
  d <- do.call("rbind", lapply(d, digit_regressors))

  ## get digits
  y <- ifelse(d$width < 0.5, 1L,
    ifelse(d$x8 < 0.15, 4L,
    ifelse(d$x30 < 0.15, 5L,
    ifelse(d$x1 > 0.55, 7L,
    ifelse((d$x7 + d$x12) > 1.05, 2L,
    ifelse((d$x17 + d$x18 + d$x19) < 0.18, 0L,
    ifelse((d$x4 + d$x11) < 0.5, 3L,
    ifelse((d$x4 + d$x32) < 1.05, 8L, 
    ifelse(d$x31 > 0.4, 9L,
    6L)))))))))
  y <- paste(y, collapse = "")

  if(tesseract) {
    y2 <- tesseract(z)
    if(y != y2) cat(sprintf("(%s != %s)", y2, y))
  }
  return(y)
}

read_nops_answers <- function(x, threshold = c(0.04, 0.42), size = 0.029, n = 45L, adjust = FALSE)
{
  ## adjustment for coordinates (e.g. for reading 2nd string page)
  if(identical(adjust, TRUE)) adjust <- c(0.4243, -0.50025)
  if(identical(adjust, FALSE)) adjust <- c(0, 0)

  ## number of answer fields to read
  if(!(is.numeric(n) && isTRUE(n %in% 1L:45L))) n <- 45L

  ## 1-15
  coord1 <- cbind(0.5532 + rep(0L:2L, each = 25L) * 0.147 + rep(0L:4L, each = 5L) * 0.027,
    0.04125 + rep(0L:4L, 15L) * 0.047)
  ## 16-30
  coord2 <- coord1 + cbind(rep(0, 5 * 15), 0.374)
  ## 31-45
  coord3 <- coord2 + cbind(rep(0, 5 * 15), 0.374 * 60/64)
  coord <- rbind(coord1, coord2, coord3)

  ## ## zap numbers next to the boxes
  ## subimage(x, c(0.7542,         0.0095), c(0.42, 0.019)) <- 0L
  ## subimage(x, c(0.7542, 0.373 + 0.0095), c(0.42, 0.019)) <- 0L
  ## subimage(x, c(0.7542, 0.723 + 0.0095), c(0.42, 0.019)) <- 0L

  y <- matrix(sapply(1:(n * 5L), function(i)
    has_mark(subimage(x, coord[i,] - adjust, size), threshold = threshold)), ncol = 5L, byrow = TRUE)
  rval <- paste(apply(y, 1, paste, collapse = ""), collapse = " ")
  if(n < 45L) rval <- paste(rval, paste(rep.int("00000", 45L - n), collapse = " "))
  return(rval)
}

read_nops_registration <- function(x, threshold = c(0.04, 0.42), size = 0.029)
{
  coord <- cbind(0.166 + rep(0L:9L, each = 7L) * 0.027,
    0.681 + rep(0L:6L, 10L) * 0.047)

  y <- try(matrix(sapply(1:nrow(coord), function(i)
    has_mark(subimage(x, coord[i,], size), threshold = threshold, fuzzy = TRUE)), ncol = 7L, byrow = TRUE),
    silent = TRUE)
  if(inherits(y, "try-error")) return("0000000")
  if(!all(apply(y, 2, function(z) any(z > 0)))) return("0000000")
  paste(apply(y, 2, which.max) - 1, collapse = "")
}

read_nops_backup <- function(x, threshold = 0.15, size = 0.01)
  format(as.numeric(mean(subimage(x, c(0.381, 0.574), size)) > threshold[1L]))

## crude tesseract interface
tesseract <- function(x, digits = TRUE) {
  writeBin(png::writePNG(1 - x), ".tesseract-temp-image.png")
  system("tesseract -psm 6 .tesseract-temp-image.png .tesseract-temp-text",
    ignore.stderr = TRUE)
  rval <- readLines(".tesseract-temp-text.txt")
  file.remove(c(".tesseract-temp-image.png", ".tesseract-temp-text.txt"))

  if(digits) {
    rval <- gsub(" ", "", rval, fixed = TRUE)
    rval <- gsub(",", "", rval, fixed = TRUE)
    rval <- gsub(".", "", rval, fixed = TRUE)
    rval <- gsub("_", "", rval, fixed = TRUE)
    rval <- gsub("x", "", rval, fixed = TRUE)
    rval <- gsub("\342", "", rval, fixed = TRUE)
    rval <- gsub("\200", "", rval, fixed = TRUE)
    rval <- gsub("\230", "", rval, fixed = TRUE)
    rval <- gsub("O", "0", rval, fixed = TRUE)
    rval <- gsub("C", "0", rval, fixed = TRUE)
    rval <- gsub("D", "0", rval, fixed = TRUE)
    rval <- gsub("Q", "0", rval, fixed = TRUE)
    rval <- gsub("o", "0", rval, fixed = TRUE)
    rval <- gsub("c", "0", rval, fixed = TRUE)
    rval <- gsub("I", "1", rval, fixed = TRUE)
    rval <- rval[nchar(rval) > 0]

    if(length(rval) > 1L) paste("ERROR:", rval[1L], sep = "")
  }

  return(rval)
}

## simple plotting function
imageplot <- function(x, ...) {
  d <- dim(x)
  xcoord <- t(which(x > 0.5, arr.ind = TRUE))
  xcoord <- t(xcoord/d)
  par(mar = rep(1, 4))
  plot(0, 0, type = "n", xlab = "", ylab = "", axes = FALSE, xlim = c(0, 1), ylim = c(0, 1), ...)
  if(prod(dim(xcoord)) > 0L) rect(xcoord[,2L] - 1/d[2L], 1 - (xcoord[,1L] - 1/d[1L]),
    xcoord[,2L], 1 - xcoord[,1L], col = "black", border = "transparent") 
}
