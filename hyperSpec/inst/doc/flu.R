### R code from vignette source 'flu.Rnw'

###################################################
### code chunk number 1: startup
###################################################
#system ("mkdir fig")
options(SweaveHooks=list(fig=function() {
  par(mar = c (4.1, 4.1, 1, .6))

  trellis.pars <- trellis.par.get ("layout.heights")
  trellis.pars [grep ("padding", names (trellis.pars))] <- 0
  trellis.par.set(layout.heights = trellis.pars)
  
  trellis.pars <- trellis.par.get ("layout.widths")
  trellis.pars [grep ("padding", names (trellis.pars))] <- 0
  trellis.par.set(layout.widths = trellis.pars)
}))
options ("width" = 100, "digits" = 5)
library (hyperSpec)

# redefine lattice functions so that the result is printed without external print command
setMethod ("plot",
           signature (x = "hyperSpec", y = "character"),
           function (x, y, ...){
             tmp <- hyperSpec:::.plot (x, y, ...)
             if (is (tmp, "trellis"))
               print (tmp)
             invisible (tmp)
           })


plotmap <- function (...) print (hyperSpec:::plotmap (...))

setMethod ("levelplot", signature (x = "hyperSpec", data = "missing"),
   function (x, data, ...) {
	   l <- hyperSpec:::.levelplot (x = formula (spc ~ x * y), data = x, ...)
		print (l)
	}
)

setMethod ("levelplot", signature (x = "formula", data = "hyperSpec"), 
   function (x, data, ...) print (hyperSpec:::.levelplot (x, data, ...))
)

plotc <- function (...){
   call <- match.call () 
   call [[1]] <- hyperSpec:::plotc 
   print (eval (call))
}

ploterrormsg <- function (fn, pkg) {
  plot (0, 0, type = "n", axes = FALSE, bty = "n", xlab = "", ylab = "")
  text (0, 0, paste ("Function", fn, "not available:\npackage", pkg, "needed."))
}
griderrormsg <- function (fn, pkg) {
  require (grid)
  grid.text (label = paste ("Function", fn, "not available:\npackage", pkg, "needed."))
  NA
}
texterrormsg <- function (fn, pkg) {
  cat ("Function", fn, "not available:\npackage", pkg, "needed.\n")
}

nice.paste <- function (...){
  fnames <- c (...)
  
  if (length (fnames) == 2L)
    fnames <- paste (fnames, collapse = " and ")
  if (length (fnames) > 1L){
    fnames [length (fnames)] <- paste ("and", tail (fnames, 1))
    fnames <- paste (fnames, collapse = ", ")
  }

  fnames
}

check.req.pkg <- function (pkg = stop ("pkg needed"), 
                           texterrors = NULL, ploterrors = NULL, griderrors = NULL,
                           hynstext = NULL, hynsplot = NULL, hynsgrid = NULL, 
                           donothing = NULL, special = NULL, v = TRUE){
  if (v) cat ("\\item[\\Rpackage{", pkg, "}:] ", sep = "")
  
  dummies <- list ()
  
  if (pkg.exists (pkg)){
    if (v) cat ("available\n")
  } else {
    for (fn in as.character (texterrors))
      dummies <- c (dummies, bquote (.(fn) <- function (...) texterrormsg (.(fn), .(pkg))))
    for (fn in as.character (ploterrors))
      dummies <- c (dummies, bquote (.(fn) <- function (...) ploterrormsg (.(fn), .(pkg))))
    for (fn in as.character (griderrors))
      dummies <- c (dummies, bquote (.(fn) <- function (...) griderrormsg (.(fn), .(pkg))))

    for (fn in as.character (hynstext))
      assignInNamespace (x = fn, 
                         value = eval (bquote (function (...) texterrormsg (.(fn), .(pkg)))), 
                         ns = "hyperSpec")
    for (fn in as.character (hynsplot))
      assignInNamespace (x = fn, 
                         value = eval (bquote (function (...) ploterrormsg (.(fn), .(pkg)))), 
                         ns = "hyperSpec")
    for (fn in as.character (hynsgrid))
      assignInNamespace (x = fn, 
                         value = eval (bquote (function (...) griderrormsg (.(fn), .(pkg)))), 
                         ns = "hyperSpec")

    fnames <- nice.paste (texterrors, ploterrors, griderrors, hynstext, hynsplot, hynsgrid, names (special))
    if (v && length (fnames) > 0L) cat (fnames, "replaced.")
    
    for (fn in as.character (donothing))
      dummies <- c (dummies, bquote (.(fn) <- function (...) invisible (NULL)))
    
    fnames <- nice.paste (donothing)
    if (v && length (fnames) > 0L) cat (fnames, "missing.")
    
    if (v) cat ("\n")
  }
  
  invisible (dummies)
}

plotvoronoi <- function (...) print (hyperSpec:::plotvoronoi (...))

# set standardized color palettes 
seq.palette <- colorRampPalette (c ("white", "dark green"), space = "Lab")

YG.palette <- function (n = 20) rgb (colorRamp (c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", 
                                             "#41AB5D", "#238B45", "#006D2C", "#00441B"), space = "Lab") 
                                # was: brewer.pal (9, "Greens")
                                (seq (1/3, 1, length.out = n)^2), maxColorValue = 255)

										  
div.palette <- colorRampPalette (c("#00008B", "#351C96", "#5235A2", "#6A4CAE", "#8164BA", "#967CC5", 
                                   "#AC95D1", "#C1AFDC", "#D5C9E8", "#E0E3E3", "#F8F8B0", "#F7E6C2", 
											  "#EFCFC6", "#E6B7AB", "#DCA091", "#D08977", "#C4725E", "#B75B46",
											  "#A9432F", "#9A2919", "#8B0000"), space = "Lab")

pkgSuggests <- function (...)
  strsplit (packageDescription (..., fields="Suggests"), ",\\s*")[[1]]

pkg.exists <- function (pkg = stop ("package name needed"), lib.loc = NULL){
  dir <- sapply (pkg, function (p) system.file (package = p, lib.loc = lib.loc))
  nzchar (dir) > 0L 
}
  
is.basepkg <- function (pkg){
  pkg.exists (pkg) && grepl ("^base$", packageDescription (pkg, fields = "Priority"))
}

pkg.or.base <- function (pkg){
  pkg [sapply (pkg, is.basepkg)] <- "base"
  
  pkg
}

citation.or.file <- function (pkg, svd.cit = sprintf ("%s.CITATION", pkg)){
  if (pkg.exists (pkg))
    citation (pkg)
  else if (file.exists (svd.cit))
    readCitationFile (file = svd.cit)
  else
    NULL
}

make.cite.keys <- function (pkg, entries){
  pkg <- pkg.or.base (pkg)

  if (! pkg.exists (pkg))
    return (pkg)
  
  if (missing (entries))
    entries <- citation.or.file (pkg)
  
  keys <- sapply (unclass (entries), attr, "key")
  
  noname <- which (sapply (keys, is.null))

  if (length (keys) == 1L && noname == 1L) {
    keys <- pkg
  } else {
    for (i in noname)
      keys [[i]] <- paste (pkg, i, sep = ".")
  }

  keys <- make.unique (unlist (keys))
  
  keys
}
  
citation.with.key <- function (pkg = "base"){
  pkg <- pkg.or.base (pkg)

  tmp <- citation.or.file (pkg)
  
  keys <- make.cite.keys (pkg, tmp)

  for (entry in seq_along (tmp))
    tmp [entry]$"key" <- keys [[entry]]

  tmp
}

cite.pkg <- function (p, entries, citefun = "cite"){
  paste ("\\\\", citefun, "{", paste (make.cite.keys (p, entries), collapse = ", "), "}", sep = "")
}

make.bib <- function (..., file = NULL) {
  pkg <- c (...)

  if (length (pkg) == 0L) {
    pkg <- loadedNamespaces()
 
    pkg <- unique (pkg.or.base (pkg))
  }
  
  l <- lapply (pkg, citation.with.key)
  l <- do.call ("c", l [! sapply (l, is.null)])

  if (!is.null (file))
    if (is.null (l))
      cat (NULL, file = file)           # touches file
    else
      cat (toBibtex (l), file = file, sep = "\n")
  
  invisible (l)
}



###################################################
### code chunk number 2: mailme
###################################################
cat ("\\newcommand{\\mailme}{\\href{mailto:", 
     packageDescription ("hyperSpec")$Maintainer, 
	  "}{\\texttt{", 
	  packageDescription ("hyperSpec")$Maintainer,
	  "}}}\n", 
	  sep = "")


###################################################
### code chunk number 3: listfunctions
###################################################
texListFun <- function (pattern){
  funs <- ls (envir = getNamespace ("hyperSpec"), pattern = pattern)
  funs <- paste ("\\\\Rfunction{", funs, "}", sep ="")
  nice.paste (funs)
}


###################################################
### code chunk number 4: cleanup (eval = FALSE)
###################################################
## sessionInfo ()
## rm (list = ls ())
## library (tools)


###################################################
### code chunk number 5: read.txt.PE
###################################################
source ("scan.txt.PerkinElmer.R")
flu <- scan.txt.PerkinElmer ("rawdata/flu?.txt", skip = 54)


###################################################
### code chunk number 6: rawspc
###################################################
flu


###################################################
### code chunk number 7: rawfig
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (flu)


###################################################
### code chunk number 8: newdata
###################################################
flu$c <- seq (from = 0.05, to = 0.30, by = 0.05)
labels (flu, "c") <- "c / (mg / l)"
flu
save (flu, file = 'flu.rda')


###################################################
### code chunk number 9: newc
###################################################
flu$c


###################################################
### code chunk number 10: delcol
###################################################
flu$file <- NULL


###################################################
### code chunk number 11: calplot1
###################################################
getOption("SweaveHooks")[["fig"]]()
plotc (flu[,,450])


###################################################
### code chunk number 12: cutspc
###################################################
flu <- flu [,,450]
labels (flu, "spc") <- expression (I ["450 nm"] / a.u.)


###################################################
### code chunk number 13: calplot2
###################################################
getOption("SweaveHooks")[["fig"]]()
plotc (flu, xlim = range (0, flu$c), ylim = range (0, flu$spc))


###################################################
### code chunk number 14: abbrev
###################################################
flu[[]]
flu$.
flu$..


###################################################
### code chunk number 15: cal
###################################################
calibration <- lm (c ~ spc, data = flu$.)


###################################################
### code chunk number 16: summarymodel
###################################################
summary (calibration)


###################################################
### code chunk number 17: pred
###################################################
I <- c (125, 400)
conc <- predict (calibration, newdata = list (spc = as.matrix(I)), interval = "prediction", 
                 level = .99)
conc


###################################################
### code chunk number 18: calplot3
###################################################
getOption("SweaveHooks")[["fig"]]()

int <- list (spc = as.matrix(seq (min (flu), max(flu), length.out = 25)))
ci <- predict (calibration, newdata = int, interval = "confidence", level = 0.99)

panel.ci <-  function (x, y, ...,
                       intensity, ci.lwr, ci.upr, ci.col = "#606060") {
   panel.xyplot (x, y, ...)
   panel.lmline (x, y,...)
   panel.lines (ci.lwr, intensity, col = ci.col)
   panel.lines (ci.upr, intensity, col = ci.col)
}

plotc (flu, panel = panel.ci,
       intensity = int$spc, ci.lwr = ci [, 2], ci.upr = ci [, 3])
## # extrapolate to lower intensities
## int <- list (spc = as.matrix(0 : min (flu)))
## ci <- predict (calibration, newdata = int, interval = "confidence", level = 0.99)
## matlines (ci, int$spc, col = c ("red","#606060","#606060"), lty = 3)

# our example
#lines (conc[-1], rep(I, 2), col = "blue")
#points (conc[1], I, col = "blue", pch = 4, cex = 0.5)


###################################################
### code chunk number 19: calplot4.1
###################################################
flu$type <- "data points"


###################################################
### code chunk number 20: calplot4.2
###################################################
tmp <- new ("hyperSpec", spc = as.matrix(seq (min (flu), max(flu), length.out = 25)),
                         wavelength = 450)
ci <-  predict (calibration, newdata = tmp$., interval = "confidence", level = 0.99)
tmp <- tmp [rep (seq (tmp, index = TRUE), 3)]
tmp$c <- as.numeric (ci)
tmp$type <- rep (colnames (ci), each = 25)

flu <- rbind (flu, tmp)


###################################################
### code chunk number 21: calplot4
###################################################
getOption("SweaveHooks")[["fig"]]()
panel.predict <- function (x, y, ..., 
                 intensity, ci, pred.col = "red", pred.pch = 19, pred.cex = 1) {
   panel.xyplot (x, y, ...)
   mapply (function (i, lwr, upr, ...) {
                 panel.lines (c (lwr, upr), rep (i, 2), ...)
              }, 
           intensity, ci [, 2], ci [, 3], MoreArgs = list (col = pred.col))
   panel.xyplot (ci [, 1], intensity, col = pred.col, pch = pred.pch, cex = pred.cex, type = "p")
}


plotc (flu, groups = type, type = c("l", "p"),
       col = c ("black", "black", "#606060", "#606060"), 
       pch = c (19, NA, NA, NA), cex = 0.5, 
       lty = c (0, 1, 1, 1),
       panel = panel.predict,
       intensity = I,
       ci = conc,
       pred.cex = 0.5)


###################################################
### code chunk number 22: flu.Rnw:240-241
###################################################
sessionInfo ()
rm (list = ls ())
library (tools)


