### R code from vignette source 'plotting.Rnw'
### Encoding: UTF-8

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
### code chunk number 5: check-required
###################################################
msg <- function (...) system (sprintf ("echo '%s'", paste (...)))

required.pkgs <- c("latticeExtra", "deldir", "rgl", "ggplot2")

dummies <- check.req.pkg ("latticeExtra", griderrors = "panel.levelplot.points",
                          hynsgrid = "plotvoronoi")
for (i in seq_along (dummies))
  eval (dummies [[i]])

check.req.pkg ("deldir", hynsgrid = "plotvoronoi")

for (p in required.pkgs [! required.pkgs %in% c("latticeExtra", "deldir")])
  check.req.pkg (p, donothing = "")


###################################################
### code chunk number 6: preproc-chondro
###################################################
chondro.preproc <- chondro - spc.fit.poly.below (chondro)
chondro.preproc <- chondro.preproc / rowMeans (chondro)
chondro.preproc <- chondro.preproc - quantile (chondro, 0.05)

cluster.cols <- c ("dark blue", "orange", "#C02020")
cluster.meansd <- aggregate (chondro.preproc, chondro$clusters, mean_pm_sd)
cluster.means  <- aggregate (chondro.preproc, chondro$clusters, mean)


###################################################
### code chunk number 7: plotspc
###################################################
getOption("SweaveHooks")[["fig"]]()
plotspc (flu)


###################################################
### code chunk number 8: plotspc
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmat (flu)


###################################################
### code chunk number 9: plotflu
###################################################
getOption("SweaveHooks")[["fig"]]()
plotc (flu)


###################################################
### code chunk number 10: levelplot
###################################################
getOption("SweaveHooks")[["fig"]]()
levelplot (spc ~ x * y, chondro, aspect = "iso")


###################################################
### code chunk number 11: plotmap
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmap (chondro)


###################################################
### code chunk number 12: voronoi
###################################################
getOption("SweaveHooks")[["fig"]]()
plotvoronoi (sample (chondro, 300), clusters ~ x * y)


###################################################
### code chunk number 13: plotspcflu
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (flu, "spc")


###################################################
### code chunk number 14: plotchomean
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (chondro.preproc, "spcmeansd")


###################################################
### code chunk number 15: plotchoprctl
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (chondro.preproc, "spcprctile")


###################################################
### code chunk number 16: plotchoprctl5
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (chondro.preproc, "spcprctl5")


###################################################
### code chunk number 17: plotflu2
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (flu, "c")


###################################################
### code chunk number 18: plotts
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (laser [,, 405], "ts")


###################################################
### code chunk number 19: plotdepth
###################################################
getOption("SweaveHooks")[["fig"]]()
depth.profile <- new ("hyperSpec",
    spc = as.matrix (rnorm (20) + 1:20),
    data = data.frame (z = 1 : 20),
    labels = list (spc = "I / a.u.", 
       z = expression (`/` (z, mu*m)),
       .wavelength = expression (lambda)))
plot (depth.profile, "depth")


###################################################
### code chunk number 20: plotmat
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (laser, "mat")


###################################################
### code chunk number 21: plotting.Rnw:231-232
###################################################
plotmat (laser)


###################################################
### code chunk number 22: plotting.Rnw:235-236
###################################################
levelplot (spc ~ .wavelength * .row, laser)


###################################################
### code chunk number 23: plotmapcho2
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (chondro, "map")


###################################################
### code chunk number 24: plotvoronoi
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (sample (chondro, 300), "voronoi")


###################################################
### code chunk number 25: wavelength
###################################################
getOption("SweaveHooks")[["fig"]]()
plotspc (paracetamol [,, 700 ~ 1200])


###################################################
### code chunk number 26: wavelength-2
###################################################
getOption("SweaveHooks")[["fig"]]()
plotspc (paracetamol, 
         wl.range = c (300 ~ 1800, 2800 ~ max), 
         xoffset = 750)


###################################################
### code chunk number 27: abscissa
###################################################
getOption("SweaveHooks")[["fig"]]()
plotspc (paracetamol, wl.reverse = TRUE )


###################################################
### code chunk number 28: colours
###################################################
getOption("SweaveHooks")[["fig"]]()
plotspc (flu, col = matlab.dark.palette (6))


###################################################
### code chunk number 29: dots
###################################################
getOption("SweaveHooks")[["fig"]]()
plotspc (paracetamol [,, 2800 ~ 3200], 
         lines.args = list (pch = 20, type = "p"))


###################################################
### code chunk number 30: mass
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (barbiturates [[1]], lines.args = list (type = "h"))


###################################################
### code chunk number 31: add
###################################################
getOption("SweaveHooks")[["fig"]]()
plotspc (chondro [ 30,,])
plotspc (chondro [300,,], add = TRUE, col = "blue")


###################################################
### code chunk number 32: sd
###################################################
getOption("SweaveHooks")[["fig"]]()
plotspc (chondro.preproc, func = sd)


###################################################
### code chunk number 33: diffline
###################################################
getOption("SweaveHooks")[["fig"]]()
plotspc (paracetamol, 
         zeroline = list (col = "red"))


###################################################
### code chunk number 34: add-line
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (laser, "spcmeansd")
abline (v = c(405.0063, 405.1121, 405.2885, 405.3591), 
        col = c("black", "blue", "red", "darkgreen"))


###################################################
### code chunk number 35: stacked1
###################################################
getOption("SweaveHooks")[["fig"]]()
plotspc (cluster.means, 
         col = cluster.cols,
         stacked = TRUE)


###################################################
### code chunk number 36: stacked2
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cluster.meansd, 
      stacked = ".aggregate",
      fill = ".aggregate",
      col = cluster.cols)


###################################################
### code chunk number 37: stacked3
###################################################
getOption("SweaveHooks")[["fig"]]()
plotspc (cluster.meansd, 
         yoffset = rep (0:2, each = 3), 
         col = rep (cluster.cols, each = 3))


###################################################
### code chunk number 38: stacked4
###################################################
getOption("SweaveHooks")[["fig"]]()
yoffsets <- apply (cluster.means [[]], 2, diff)
yoffsets <- - apply (yoffsets, 1, min) 
plot (cluster.means, yoffset = c (0, cumsum (yoffsets)), 
      col = cluster.cols)


###################################################
### code chunk number 39: stacked5
###################################################
getOption("SweaveHooks")[["fig"]]()
yoffset <- apply (chondro.preproc, 2, quantile, c(0.05, 0.95))
yoffset <- range (yoffset)
plot(chondro.preproc[1], 
     plot.args = list (ylim = c (0, 2) * yoffset),
     lines.args = list( type = "n"))
yoffset <- (0:1) * diff (yoffset)
for (i in 1 : 3){
  plot(chondro.preproc, "spcprctl5", yoffset = yoffset [i],
       col = "gray", add = TRUE)
  plot (chondro.preproc [i], yoffset = yoffset [i], 
        col = matlab.dark.palette (3) [i], add = TRUE, 
        lines.args = list (lwd = 2))
}


###################################################
### code chunk number 40: lin-cal-1
###################################################
getOption("SweaveHooks")[["fig"]]()
plotc (flu [,, 450])


###################################################
### code chunk number 41: lin-cal-3
###################################################
getOption("SweaveHooks")[["fig"]]()
plotc (flu, func = range, groups = .wavelength)


###################################################
### code chunk number 42: plotc2
###################################################
getOption("SweaveHooks")[["fig"]]()
plotc (flu [,, c (405, 445)], spc ~ c | .wavelength, 
       cex = .3, scales = list (alternating = c(1, 1)))


###################################################
### code chunk number 43: plotc3
###################################################
getOption("SweaveHooks")[["fig"]]()
 plotc (flu [,, c (405, 445)], groups = .wavelength)


###################################################
### code chunk number 44: lin-cal-4
###################################################
getOption("SweaveHooks")[["fig"]]()
plotc (flu [,, 450], 
       ylab = expression (I ["450 nm"] / a.u.),
       xlim = range (0, flu$c + .01), 
       ylim = range (0, flu$spc + 10),
       pch = 4)


###################################################
### code chunk number 45: lincal-panel
###################################################
panelcalibration <- function (x, y, ..., clim = range (x), level = .95) {
  panel.xyplot (x, y, ...)
  lm <- lm (y ~ x)
  panel.abline (coef (lm), ...)
  cx <- seq (clim [1], clim [2], length.out = 50)
  cy <- predict (lm, data.frame (x = cx), 
                 interval = "confidence", 
                 level = level) 
  panel.lines (cx, cy [,2], col = "gray")
  panel.lines (cx, cy [,3], col = "gray")
}


###################################################
### code chunk number 46: lin-cal-5
###################################################
getOption("SweaveHooks")[["fig"]]()
plotc (flu [,,405], panel = panelcalibration,
       pch = 4, clim = c (0, 0.35), level = .99)


###################################################
### code chunk number 47: plotc4
###################################################
getOption("SweaveHooks")[["fig"]]()
 plotc (laser [,, c(405.0063, 405.1121, 405.2885, 405.3591)], 
        spc ~ t, 
        groups = .wavelength, 
        type = "b", 
        col = c ("black", "blue", "red", "darkgreen"))


###################################################
### code chunk number 48: levelplot
###################################################
getOption("SweaveHooks")[["fig"]]()
 levelplot (spc ~ x * y, chondro)


###################################################
### code chunk number 49: levelplot-factor
###################################################
getOption("SweaveHooks")[["fig"]]()
 levelplot (clusters ~ x * y, chondro)


###################################################
### code chunk number 50: plotting.Rnw:508-509
###################################################
plot (laser, "mat", col = heat.colors (20))


###################################################
### code chunk number 51: plotmat1
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmat (laser, col = heat.colors (20))


###################################################
### code chunk number 52: plotmat1a
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmat (laser, y = "t")


###################################################
### code chunk number 53: plotting.Rnw:522-523
###################################################
plotmat (laser, y = laser$t, ylab = labels (laser, "t"))


###################################################
### code chunk number 54: plotmat2
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmat (flu, col = matlab.dark.palette (20))
plotmat (flu, col = "white", 
         contour = TRUE, add = TRUE)


###################################################
### code chunk number 55: plotmap-barb
###################################################
getOption("SweaveHooks")[["fig"]]()
require ("latticeExtra")
barb <- do.call (collapse, barbiturates[1:50])
barb <- orderwl (barb)
levelplot (spc ~ .wavelength * z, barb, 
           panel = panel.levelplot.points,
           cex = .33, col.symbol = NA,
           col.regions = matlab.palette)


###################################################
### code chunk number 56: plotmap-chondro
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmap (chondro)


###################################################
### code chunk number 57: plotmap-yx
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmap (chondro, spc ~ y * x)


###################################################
### code chunk number 58: plotmap-clu
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmap (chondro, clusters ~ x * y)


###################################################
### code chunk number 59: plotmap-col
###################################################
getOption("SweaveHooks")[["fig"]]()
print (plotmap (chondro, clusters ~ x * y,
                col.regions = cluster.cols))


###################################################
### code chunk number 60: plotmap-wave
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmap (chondro.preproc [, , c(728, 782, 1098, 
                                1240, 1482, 1577)],
         col.regions = matlab.palette)


###################################################
### code chunk number 61: plotmap-pca
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmap (chondro, 
         spc ~ y * x | x > 5, 
         col.regions = matlab.palette(20))


###################################################
### code chunk number 62: plotmap-pca2
###################################################
getOption("SweaveHooks")[["fig"]]()
pca <- prcomp (~ spc, data = chondro.preproc$.)
scores <- decomposition (chondro, pca$x, 
                         label.wavelength = "PC", 
                         label.spc = "score /  a.u.")
plotmap (scores [,,1:2], 
         spc ~ y * x | as.factor(.wavelength), 
         func = NULL,
         col.regions = matlab.palette(20))


###################################################
### code chunk number 63: plotmap-pca3
###################################################
getOption("SweaveHooks")[["fig"]]()
levelplot (spc ~ y * x | as.factor(.wavelength), 
           scores [,,1:2], 
           aspect = "iso",
           col.regions = matlab.palette(20))


###################################################
### code chunk number 64: voronoi-2
###################################################
getOption("SweaveHooks")[["fig"]]()
plotvoronoi (sample (chondro, 300), clusters ~ x * y, 
             col.regions = matlab.palette(20))


###################################################
### code chunk number 65: missing
###################################################
getOption("SweaveHooks")[["fig"]]()
mark.missing <- function (x, y, z, ...){
  panel.levelplot (x, y, z, ...)

  miss <- expand.grid (x = unique (x), y = unique (y))
  miss <- merge (miss, data.frame (x, y, TRUE), 
                 all.x = TRUE)
  miss <- miss [is.na (miss[, 3]),]
  panel.xyplot (miss [, 1], miss [, 2], pch = 4, ...)
}

plotmap (sample (chondro, 865), 
         col.regions = matlab.palette(20),
         col = "black",
         panel = mark.missing)


###################################################
### code chunk number 66: uneven-prep
###################################################
uneven <- chondro
uneven$x <- uneven$x + round (rnorm (nrow (uneven), sd = 0.05), digits = 1)
uneven$y <- uneven$y + round (rnorm (nrow (uneven), sd = 0.05), digits = 1)


###################################################
### code chunk number 67: uneven-I
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmap (uneven)


###################################################
### code chunk number 68: uneven-II
###################################################
getOption("SweaveHooks")[["fig"]]()
plotvoronoi (uneven)


###################################################
### code chunk number 69: uneven-III
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmap (uneven, panel = panel.levelplot.points, 
				 cex = 0.75, col.symbol=NA)


###################################################
### code chunk number 70: uneven-IV
###################################################
getOption("SweaveHooks")[["fig"]]()
rx <- makeraster (uneven$x, start = -11.55, d = 1, tol = 0.3)
uneven$x <- rx$x
ry <- makeraster (uneven$y, start = -4.77, d = 1, tol = 0.3)
uneven$y <- ry$x
plotmap (uneven)


###################################################
### code chunk number 71: rgl-plot (eval = FALSE)
###################################################
## laser <- laser [,,404.8 ~ 405.6] / 10000
## laser$t <- laser$t / 3600
## cols <- rep (matlab.palette (nrow (laser)), nwl (laser))
## surface3d (y = wl (laser), x = laser$t, 
##            z = laser$spc, col =  cols)
## aspect3d (c (1, 1, 0.25))
## axes3d (c ('x+-', 'y--', 'z--'))
## axes3d ('y--', nticks = 25, labels= FALSE)
## mtext3d ("t / h", 'x+-', line = 2.5)
## mtext3d ("lambda / nm", 'y--', line = 2.5)
## mtext3d ("I / a.u.", edge = 'z--', line = 2.5)


###################################################
### code chunk number 72: rgl-do
###################################################
if (exists (".rgl") && .rgl){        # use extra argument to turn on rgl only on local machine as it doesn't work on r-forge.
  if (require (rgl)){
    open3d (windowRect=c(20,20,600, 350))  # this is needed only for automatically 
                                        # producing the snapshot
laser <- laser [,,404.8 ~ 405.6] / 10000
laser$t <- laser$t / 3600
cols <- rep (matlab.palette (nrow (laser)), nwl (laser))
surface3d (y = wl (laser), x = laser$t, 
           z = laser$spc, col =  cols)
aspect3d (c (1, 1, 0.25))
axes3d (c ('x+-', 'y--', 'z--'))
axes3d ('y--', nticks = 25, labels= FALSE)
mtext3d ("t / h", 'x+-', line = 2.5)
mtext3d ("lambda / nm", 'y--', line = 2.5)
mtext3d ("I / a.u.", edge = 'z--', line = 2.5)
     par3d (userMatrix = matrix (c (-0.52,  0.4, -0.75, 0, 
                                    -0.85, -0.28, 0.44, 0, 
                                    -0.04,  0.87, 0.49, 0, 
                                    -0.75,  0.75,     0, 1), ncol = 4L),
            scale = c (2.75, 5, 0.175),
            windowRect = c(20L, 50L, 520L, 330L),
            zoom = 0.75)
    rgl.snapshot ("fig-3D.png", fmt="png", top=TRUE )
    rgl.quit ()
  } else {
    png ("fig-3D.png")
    ploterrormsg ("", "rgl")
    dev.off ()
  }
}


###################################################
### code chunk number 73: ggplotspc (eval = FALSE)
###################################################
## qplotspc (flu) + aes (colour = c)


###################################################
### code chunk number 74: ggplotmap (eval = FALSE)
###################################################
## qplotmap (chondro) + 
##   scale_fill_gradientn ("spc", colours = matlab.palette ()) 


###################################################
### code chunk number 75: ggplotmeansd (eval = FALSE)
###################################################
## qplotspc (mean (chondro)) +
## geom_ribbon (aes (ymin = mean + sd, 
##                   ymax = mean - sd, 
##                   y = 0, group = NA), 
##              alpha = 0.25, 
##              data = as.t.df (mean_sd (chondro)))


###################################################
### code chunk number 76: ggplotspccut (eval = FALSE)
###################################################
## qplotspc (paracetamol / 1e4, 
##           wl.range = c( min ~ 1800, 2800 ~ max)) +
## 	  scale_x_continuous (breaks = seq (0, 3200, 400)) 


###################################################
### code chunk number 77: ggplot2-do
###################################################
if (require (ggplot2)){
qplotspc (flu) + aes (colour = c)
  ggsave("plotting-fig--ggplotspc.pdf", width = 4, height = 2.6)
qplotmap (chondro) + 
  scale_fill_gradientn ("spc", colours = matlab.palette ()) 
  ggsave("plotting-fig--ggplotmap.pdf", width = 4, height = 2.6)
qplotspc (mean (chondro)) +
geom_ribbon (aes (ymin = mean + sd, 
                  ymax = mean - sd, 
                  y = 0, group = NA), 
             alpha = 0.25, 
             data = as.t.df (mean_sd (chondro)))
  ggsave("plotting-fig--ggplotmeansd.pdf", width = 4, height = 2.6)
qplotspc (paracetamol / 1e4, 
          wl.range = c( min ~ 1800, 2800 ~ max)) +
	  scale_x_continuous (breaks = seq (0, 3200, 400)) 
  ggsave("plotting-fig--ggplotspccut.pdf", width = 4, height = 2.6)
} else {
  for (f in c ("ggplotspc", "ggplotmap", "ggplotmeansd", "ggplotspccut")){
    pdf (sprintf ("plotting-fig--%s.pdf", f), width = 4, height = 2.6) 
    ploterrormsg ("", "ggplot2")
    dev.off ()
  }
}


###################################################
### code chunk number 78: plotting.Rnw:838-839 (eval = FALSE)
###################################################
## spc.identify (plotspc (paracetamol, wl.range = c (600 ~ 1800, 2800 ~ 3200), xoffset = 800))


###################################################
### code chunk number 79: plotting.Rnw:847-848 (eval = FALSE)
###################################################
## map.identify (chondro)


###################################################
### code chunk number 80: plotting.Rnw:855-856 (eval = FALSE)
###################################################
## map.sel.poly (chondro)


###################################################
### code chunk number 81: plotting.Rnw:869-872 (eval = FALSE)
###################################################
## plot (laser, "mat")
## trellis.focus ()
## grid.locator ()


###################################################
### code chunk number 82: bib
###################################################
make.bib (c("latticeExtra", "rgl", "ggplot2", "playwith",  "plotrix", "deldir", "tripack"), 
          file = "plotting-pkg.bib")
print (as.matrix(Sys.info()))
sessionInfo ()
rm (list = ls ())
library (tools)


