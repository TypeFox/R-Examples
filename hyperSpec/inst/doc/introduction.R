### R code from vignette source 'introduction.Rnw'

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
check.req.pkg ("pls", special = list (msc = function (x) {texterrormsg ("msc", "pls"); x}))
check.req.pkg ("baseline", 
               special = list (
                 baseline = function (x) {texterrormsg ("baseline", "baseline"); x},
                 getCorrected = function (x) {texterrormsg ("getCorrected", "baseline"); x}
                 ))
check.req.pkg ("ggplot2", donothing = "")
check.req.pkg ("compiler", donothing = "")
check.req.pkg ("inline", donothing = "")


###################################################
### code chunk number 6: introduction.Rnw:142-144 (eval = FALSE)
###################################################
## chk.hy (object)                         
## validObject (object)                    


###################################################
### code chunk number 7: introduction.Rnw:159-160 (eval = FALSE)
###################################################
## sweep (flu, 2, mean, `-`)


###################################################
### code chunk number 8: introduction.Rnw:165-166
###################################################
`+` (3, 5)


###################################################
### code chunk number 9: introduction.Rnw:175-176 (eval = FALSE)
###################################################
## wl (flu) <- new.wavelength.values


###################################################
### code chunk number 10: init
###################################################
library ("hyperSpec")


###################################################
### code chunk number 11: checkCompleteOptionTable
###################################################
stopifnot (all (names (hy.getOptions(TRUE)) %in% c ("debuglevel", "gc", "file.remove.emptyspc", 
                  "file.keep.name", "tolerance")))


###################################################
### code chunk number 12: print
###################################################
chondro
summary (chondro)


###################################################
### code chunk number 13: nwl
###################################################
nrow (chondro)
nwl (chondro)
ncol (chondro)
dim (chondro)


###################################################
### code chunk number 14: names
###################################################
colnames (chondro)


###################################################
### code chunk number 15: introduction.Rnw:312-313 (eval = FALSE)
###################################################
## spc <- new ("hyperSpec", spc = spectra.matrix, wavelength = wavelength.vector, data = extra.data)


###################################################
### code chunk number 16: introduction.Rnw:334-336
###################################################
pcov <- pooled.cov (chondro, chondro$clusters)
rnd <- rmmvnorm (rep (10, 3), mean = pcov$mean, sigma = pcov$COV)


###################################################
### code chunk number 17: simspc
###################################################
getOption("SweaveHooks")[["fig"]]()
cluster.cols <- c ("dark blue", "orange", "#C02020")
plot (rnd, col = cluster.cols [rnd$.group])


###################################################
### code chunk number 18: lda
###################################################
require ("MASS")
rnd <- rmmvnorm (rep (200, 3), mean = pcov$mean, sigma = pcov$COV)
lda <- lda (clusters ~ spc, rnd)
 
pred.chondro <- predict (lda, chondro)
pred.sim     <- predict (lda)


###################################################
### code chunk number 19: simlda
###################################################
getOption("SweaveHooks")[["fig"]]()
colors <- c("#00008040", "#FFA50040", "#C0202040")
plot (pred.chondro$x, col = colors [chondro$clusters], pch = 3)
points (pred.sim$x, col = colors [rnd$clusters], pch = 20, cex = 0.5)


###################################################
### code chunk number 20: selspc
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (flu, col = "gray")
plot (flu [1 : 3], add = TRUE)


###################################################
### code chunk number 21: delspc
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (flu, col = "gray")
plot (flu [-3], add = TRUE)


###################################################
### code chunk number 22: selspc2
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (flu, col = "gray")
plot (flu [flu$c > 0.2], add = TRUE)


###################################################
### code chunk number 23: sample
###################################################
sample (chondro, 3)


###################################################
### code chunk number 24: isample
###################################################
isample (chondro, 3)


###################################################
### code chunk number 25: seq
###################################################
seq (chondro, length.out = 3, index = TRUE)
seq (chondro, by = 100)


###################################################
### code chunk number 26: data
###################################################
colnames (chondro)
chondro [[1 : 3, 1]]
chondro [[1 : 3, -4]]
chondro [[1 : 3, "x"]]
chondro [[1 : 3, c (TRUE, FALSE)]]      # note the recycling!


###################################################
### code chunk number 27: data2
###################################################
flu$c


###################################################
### code chunk number 28: data3
###################################################
flu$n <- list (1 : 6, label = "sample no.")


###################################################
### code chunk number 29: data2
###################################################
indexmatrix <- matrix (c (1 : 3, 1 : 3), ncol = 2)
indexmatrix
chondro [[indexmatrix, wl.index = TRUE]]
diag (chondro [[1 : 3, , min ~ min + 2i]])


###################################################
### code chunk number 30: data2
###################################################
indexmatrix <- matrix (c (1 : 3, 1 : 3), ncol = 2)
indexmatrix
chondro [[indexmatrix, wl.index = TRUE]]
diag (chondro [[1 : 3, , min ~ min + 2i]])


###################################################
### code chunk number 31: wl2ivec
###################################################
wl2i (flu, 405 : 410)


###################################################
### code chunk number 32: wl2ivec2
###################################################
wl2i (flu, 405 ~ 410)


###################################################
### code chunk number 33: wl2ivec3
###################################################
wl2i (chondro, 1000 : 1010)


###################################################
### code chunk number 34: wl2ivec4
###################################################
wl2i (chondro, 1000 ~ 1010)


###################################################
### code chunk number 35: wl2i.minmax
###################################################
wl2i (flu, min ~ 410)


###################################################
### code chunk number 36: wl2i.im
###################################################
wl2i (flu, 450 - 2i ~ 450 + 2i)
wl2i (flu, max - 2i ~ max)


###################################################
### code chunk number 37: wl2i.list
###################################################
wl2i (flu, c (min ~ 406.5, max - 2i ~ max))


###################################################
### code chunk number 38: introduction.Rnw:560-561
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (paracetamol [,, 2800 ~ 3200])


###################################################
### code chunk number 39: introduction.Rnw:566-567
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (paracetamol [,, 2800 : 3200, wl.index = TRUE])


###################################################
### code chunk number 40: introduction.Rnw:576-577
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (paracetamol [,, -(500 : 1000), wl.index = TRUE])


###################################################
### code chunk number 41: introduction.Rnw:585-586
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (paracetamol [,, c (min ~ 1750, 2800 ~ max)])


###################################################
### code chunk number 42: merged
###################################################
laser
wavelengths <- wl (laser)
frequencies <- 2.998e8 / wavelengths / 1000
wl (laser) <- frequencies
labels (laser, ".wavelength") <- "f / THz"
laser
rm (laser)


###################################################
### code chunk number 43: introduction.Rnw:616-617
###################################################
wl (laser, "f / THz") <- frequencies


###################################################
### code chunk number 44: introduction.Rnw:620-621
###################################################
wl (laser) <- list (wl = frequencies, label = "f / THz")


###################################################
### code chunk number 45: orderwl
###################################################
barb <- collapse (barbiturates [1 : 3])
wl (barb)
barb <- orderwl (barb)
wl (barb)


###################################################
### code chunk number 46: introduction.Rnw:638-642
###################################################
flu <- flu [,,400 ~ 407] # make a small and handy version of the flu data set
as.data.frame (flu)
colnames (as.data.frame (flu))
as.data.frame (flu) $ spc


###################################################
### code chunk number 47: introduction.Rnw:648-650
###################################################
flu$.
flu$..


###################################################
### code chunk number 48: introduction.Rnw:653-654
###################################################
flu [[, c ("c", "spc")]]


###################################################
### code chunk number 49: introduction.Rnw:663-664
###################################################
as.t.df (apply (flu, 2, mean_pm_sd))


###################################################
### code chunk number 50: introduction.Rnw:669-670
###################################################
head (as.long.df (flu), 20)


###################################################
### code chunk number 51: introduction.Rnw:677-679
###################################################
flu [[]]
class (flu [[]])


###################################################
### code chunk number 52: introduction.Rnw:682-683
###################################################
flu [[1:3,, 406 ~ 407]]


###################################################
### code chunk number 53: introduction.Rnw:686-687
###################################################
flu [[1:3, c ("file", "spc"), 406 ~ 407]]


###################################################
### code chunk number 54: introduction.Rnw:690-691
###################################################
rm (flu)


###################################################
### code chunk number 55: cbind
###################################################
dim (flu)
dim (cbind (flu, flu))
dim (rbind (flu, flu))


###################################################
### code chunk number 56: collapse
###################################################
barb <- collapse (barbiturates)
wl (barb) [1 : 25]


###################################################
### code chunk number 57: collapse-orderwl
###################################################
barb <- orderwl (barb)
barb [[1:3, , min ~ min + 10i]]


###################################################
### code chunk number 58: merge-sample
###################################################
chondro.low <- sample (chondro [,, 600 ~ 1200], 700)
nrow (chondro.low)
chondro.high <- sample (chondro [,, 1400 ~ 1800], 700)
nrow (chondro.high)


###################################################
### code chunk number 59: merge
###################################################
chondro.merged <- merge (chondro.low, chondro.high)
nrow (chondro.merged)


###################################################
### code chunk number 60: introduction.Rnw:767-769
###################################################
chondro.merged <- merge (chondro.low, chondro.high, all = TRUE)
nrow (chondro.merged)


###################################################
### code chunk number 61: missing
###################################################
getOption("SweaveHooks")[["fig"]]()
print (levelplot (spc ~ x * y | as.factor (paste (.wavelength, "  1/cm")), 
                  chondro.merged [,,c(1000, 1650)], 
                  aspect = "iso", col.regions = matlab.palette ()))


###################################################
### code chunk number 62: introduction.Rnw:782-786
###################################################
png ("introduction-fig--merged.png", width = 500, height = 425, res=100)
plot (chondro.merged [1 : 100], "mat")
dev.off()
rm (chondro)


###################################################
### code chunk number 63: introduction.Rnw:795-797
###################################################
merged <- merge (chondro [1:7,, 610 ~ 620], chondro [5:10,, 615 ~ 625], all = TRUE)
merged$.


###################################################
### code chunk number 64: approxfun
###################################################
approxfun <- function (y, wl, new.wl){
  approx (wl, y, new.wl, method = "constant",
          ties = function (x) mean (x, na.rm = TRUE)
          )$y
}


###################################################
### code chunk number 65: introduction.Rnw:820-824
###################################################
merged <- apply (merged, 1, approxfun, 
                 wl = wl (merged), new.wl = unique (wl (merged)), 
                 new.wavelength = "new.wl")
merged$.


###################################################
### code chunk number 66: cut.wl
###################################################
flu [,, min ~ 408.5]
flu [[,, c (min ~ min + 2i, max - 2i ~ max)]]


###################################################
### code chunk number 67: introduction.Rnw:879-881
###################################################
tmp <- chondro
wl (tmp) <- wl (tmp) - 10


###################################################
### code chunk number 68: shift-wl
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (chondro [135])
plot (tmp [135,,], add = TRUE, col = "red")


###################################################
### code chunk number 69: fun-interpolate
###################################################
interpolate <- function (spc, shift, wl){
  spline (wl + shift, spc, xout = wl, method = "natural")$y
}


###################################################
### code chunk number 70: introduction.Rnw:904-905
###################################################
tmp <- apply (chondro, 1, interpolate, shift = -10, wl = wl (chondro))


###################################################
### code chunk number 71: shift-interp
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (chondro [135])
plot (tmp [135], add = TRUE, col = "red")


###################################################
### code chunk number 72: shift-untsch
###################################################
getOption("SweaveHooks")[["fig"]]()
tmp <- chondro [135,, 990 ~ 1010]
plot (tmp, lines.args = list (type = "b", pch = 19, cex = 0.5))
wl (tmp) <- wl(tmp) - 0.5
plot (tmp, lines.args = list (type = "b", pch = 19, cex = 0.5), add = TRUE, col = "red")
tmp <- chondro [135]
tmp <- apply (tmp, 1, function (x, wl, shift)
              spline (wl + shift, x, xout = wl)$y,
              wl = wl (tmp), shift = -0.5)
plot (tmp, lines.args = list (type = "b", pch = 19, cex = 0.5), add = TRUE, col = "blue")


###################################################
### code chunk number 73: introduction.Rnw:929-934
###################################################
shifts <- rnorm (nrow (chondro))
tmp <- chondro [[]]
for (i in seq_len (nrow (chondro)))
  tmp [i, ] <- interpolate (tmp [i, ], shifts [i], wl = wl (chondro))
chondro [[]] <- tmp


###################################################
### code chunk number 74: introduction.Rnw:956-972
###################################################

find.max <- function (y, x){
  pos <- which.max (y) + (-1:1)
  X <- x [pos] - x [pos [2]]
  Y <- y [pos] - y [pos [2]]
  
  X <- cbind (1, X, X^2)
  coef <- qr.solve (X, Y)

  - coef [2] / coef [3] / 2 + x [pos [2]]
}

bandpos <- apply (chondro [[,, 990 ~ 1020]], 1, find.max,  wl (chondro [,, 990 ~ 1020]))
refpos <- find.max (colMeans (chondro[[,, 990 ~ 1020]]),  wl (chondro [,, 990 ~ 1020]))

shift1 <- refpos - bandpos


###################################################
### code chunk number 75: introduction.Rnw:977-979
###################################################
chondro <- chondro - spc.fit.poly.below (chondro [,,min+3i ~ max - 3i], chondro)
chondro <- sweep (chondro, 1, rowMeans (chondro [[]], na.rm = TRUE), "/")


###################################################
### code chunk number 76: introduction.Rnw:983-994
###################################################
targetfn <- function (shift, wl, spc, targetspc){
  error <- spline (wl + shift, spc, xout = wl)$y - targetspc
  sum (error^2)
}

shift2 <- numeric (nrow (chondro))
tmp <- chondro [[]]
target <- colMeans (chondro [[]])
for (i in 1 : nrow (chondro))
  shift2 [i] <- unlist (optimize (targetfn, interval = c (-5, 5), wl = chondro@wavelength, 
                                 spc = tmp[i,], targetspc = target)$minimum)


###################################################
### code chunk number 77: shift-fit
###################################################
getOption("SweaveHooks")[["fig"]]()
df <- data.frame (shift = c (shifts, shifts + shift1, shifts + shift2), 
                method = rep (c ("original", "find maximum", "interpolation"), 
                  each = nrow (chondro)))
plot (histogram (~ shift | method, data = df, breaks = do.breaks(range (df$shift), 25),
           layout = c (3,1)))


###################################################
### code chunk number 78: introduction.Rnw:1023-1027 (eval = FALSE)
###################################################
## ir.spc <- chondro / 1500 ## fake IR data
## high.int <- apply (ir.spc > 1, 1, any) # any point above 1 is bad
## low.int <- apply (ir.spc, 1, max) < 0.1 # the maximum should be at least 0.1
## ir.spc <- ir.spc [! high.int & ! low.int] 


###################################################
### code chunk number 79: introduction.Rnw:1031-1040
###################################################
mean_sd_filter <- function (x, n  = 5) {
  x <- x - mean (x)
  s <- n * sd (x) 
  (x <= s) & (x > -s)
}

OK <- apply (chondro [[]], 2, mean_sd_filter, n = 4) # logical matrix

spc.OK <- chondro [apply (OK, 1, all)] 


###################################################
### code chunk number 80: filter
###################################################
getOption("SweaveHooks")[["fig"]]()

plot (chondro [! apply (OK, 1, all)])
i <- which (! OK, arr.ind = TRUE)
points (wl (chondro) [i [,2]], chondro[[!OK]], pch = 19, col = "red", cex = 0.5)


###################################################
### code chunk number 81: introduction.Rnw:1059-1062
###################################################
spc <- chondro [1 : 3,, min ~ min + 15i]
spc [[cbind (1:3, sample (nwl (spc), 3)), wl.index = TRUE]] <- 0
spc [[]]


###################################################
### code chunk number 82: introduction.Rnw:1066-1068
###################################################
spc [[spc < 1e-4]] <- NA
spc [[]]


###################################################
### code chunk number 83: introduction.Rnw:1077-1079
###################################################
spc.corrected <- spc.NA.linapprox (spc)
spc.corrected [[]]


###################################################
### code chunk number 84: bad
###################################################
getOption("SweaveHooks")[["fig"]]()
spc [[is.na (spc)]] <- 0
plot (spc)
spc [[spc < 1e-4]] <- NA
plot (spc.NA.linapprox (spc), add = TRUE, col = "blue", lines.args = list (type = "b", pch = 19, cex = 0.5))


###################################################
### code chunk number 85: fig-loess
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (paracetamol, wl.range = c (300 ~ 1800, 2800 ~ max), xoffset = 850)
p <- spc.loess (paracetamol, c(seq (300, 1800, 2), seq (2850, 3150, 2)))
plot (p, wl.range = c (300 ~ 1800, 2800 ~ max), xoffset = 850, col = "red", add = TRUE)
b <- spc.bin (paracetamol, 4)
plot (b, wl.range = c (300 ~ 1800, 2800 ~ max), xoffset = 850, 
      lines.args = list (pch = 20, cex = .3, type = "p"), col = "blue", add = TRUE)


###################################################
### code chunk number 86: fig-loess-kl
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (paracetamol [, , 1600 ~ 1670])
plot (p [, , 1600 ~ 1670], col = "red", add = TRUE)
plot (b [, , 1600 ~ 1670], col = "blue", add = TRUE)


###################################################
### code chunk number 87: ofs
###################################################
offsets <- apply (chondro, 1, min)
chondro.offset.corrected <- sweep (chondro, 1, offsets, "-")


###################################################
### code chunk number 88: ofs2
###################################################
chondro.offset.corrected <- sweep (chondro, 1, min, "-")


###################################################
### code chunk number 89: bl
###################################################
bl <- spc.fit.poly.below (chondro)
chondro <- chondro - bl


###################################################
### code chunk number 90: do-bl
###################################################
corrected <- hyperSpec::chondro [1] # start with the unchanged data set

require ("baseline")
bl <- baseline (corrected [[]], method = "modpolyfit", degree = 4)
corrected [[]] <- getCorrected (bl)


###################################################
### code chunk number 91: introduction.Rnw:1198-1202
###################################################
getOption("SweaveHooks")[["fig"]]()
baseline <- corrected 
baseline [[]] <- getBaseline (bl)
plot (hyperSpec::chondro [1], plot.args = list (ylim = range (hyperSpec::chondro [1], 0)))
plot (baseline, add = TRUE, col = "red")


###################################################
### code chunk number 92: introduction.Rnw:1206-1207
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (corrected, plot.args = list (ylim = range (hyperSpec::chondro [1], 0)))


###################################################
### code chunk number 93: introduction.Rnw:1215-1216
###################################################
rm (bl, chondro)


###################################################
### code chunk number 94: normalize1
###################################################
chondro <- sweep (chondro, 1, mean, "/")


###################################################
### code chunk number 95: norm
###################################################
factors <- 1 / apply (chondro [, , 1600 ~ 1700], 1, mean)
chondro <- sweep (chondro, 1, factors, "*")


###################################################
### code chunk number 96: centre-flu
###################################################
getOption("SweaveHooks")[["fig"]]()
flu.centered <- scale (flu, scale = FALSE)
plot (flu.centered)


###################################################
### code chunk number 97: perc
###################################################
getOption("SweaveHooks")[["fig"]]()
chondro <- scale (chondro, center = quantile (chondro, 0.05), scale = FALSE)
plot (chondro, "spcprctl5")


###################################################
### code chunk number 98: msc (eval = FALSE)
###################################################
## require (pls)
## chondro.msc <- chondro
## chondro.msc [[]] <- msc (chondro [[]])


###################################################
### code chunk number 99: label (eval = FALSE)
###################################################
## labels (absorbance.spectra)$spc <- "A"


###################################################
### code chunk number 100: pca
###################################################
pca <- prcomp (~ spc, data = chondro$., center = FALSE)


###################################################
### code chunk number 101: pca-auto
###################################################
pca <- prcomp (~ spc, data = chondro, center = FALSE)


###################################################
### code chunk number 102: decomp
###################################################
scores <- decomposition (chondro, pca$x, label.wavelength = "PC", 
                         label.spc = "score / a.u.")
scores


###################################################
### code chunk number 103: loadings
###################################################
loadings <- decomposition (chondro, t(pca$rotation), scores = FALSE, 
                           label.spc = "loading I / a.u.")
loadings


###################################################
### code chunk number 104: retain.col
###################################################
loadings <- decomposition (chondro, t(pca$rotation), scores = FALSE, 
                           retain.columns = TRUE, label.spc = "loading I / a.u.")
loadings[1]$..


###################################################
### code chunk number 105: retain
###################################################
chondro$measurement <- 1
loadings <- decomposition (chondro, t(pca$rotation), scores = FALSE, 
                           label.spc = "loading I / a.u.")
loadings[1]$..


###################################################
### code chunk number 106: pca-load
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (loadings [1:3], stacked = TRUE)


###################################################
### code chunk number 107: pca-score
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmap (scores [,,3], col.regions = div.palette (20))


###################################################
### code chunk number 108: pca-smooth
###################################################
smoothed <- scores [,, 1:10] %*% loadings [1:10]


###################################################
### code chunk number 109: ggplot (eval = FALSE)
###################################################
## require (ggplot2)
## ggplot (as.long.df (chondro [1]), aes (x = .wavelength, y = spc)) + geom_line ()                           


###################################################
### code chunk number 110: ggplot-do
###################################################
if (require (ggplot2)){
require (ggplot2)
ggplot (as.long.df (chondro [1]), aes (x = .wavelength, y = spc)) + geom_line ()                           
  ggsave (file = "introduction-fig--ggplot-do.pdf", width = 8, height = 4)
} else {
  pdf ("introduction-fig--ggplot-do.pdf", width = 8, height = 4)
  ploterrormsg (NULL, "ggplot2")
  dev.off()
}


###################################################
### code chunk number 111: hca
###################################################
dist <- pearson.dist (chondro [[]])


###################################################
### code chunk number 112: hca-asmatrix
###################################################
dist <- pearson.dist (chondro) 
dendrogram <- hclust (dist, method = "ward")


###################################################
### code chunk number 113: dend
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (dendrogram)


###################################################
### code chunk number 114: dendcut
###################################################
chondro$clusters <- as.factor (cutree (dendrogram, k = 3))


###################################################
### code chunk number 115: clustname
###################################################
levels (chondro$clusters) <- c ("matrix", "lacuna", "cell")


###################################################
### code chunk number 116: clustmap
###################################################
getOption("SweaveHooks")[["fig"]]()
plotmap (chondro, clusters ~ x * y, col.regions = cluster.cols)


###################################################
### code chunk number 117: clustmean
###################################################
getOption("SweaveHooks")[["fig"]]()
means <- aggregate (chondro, by = chondro$clusters, mean_pm_sd)
plot (means, col = cluster.cols, stacked = ".aggregate", fill = ".aggregate")


###################################################
### code chunk number 118: split
###################################################
clusters <- split (chondro, chondro$clusters)
clusters


###################################################
### code chunk number 119: speed1
###################################################
tmp <- chondro [1 : 50]
shifts <- rnorm (nrow (tmp))
system.time ({
  for (i in seq_len (nrow (tmp)))
    tmp [[i]] <- interpolate (tmp [[i]], shifts [i], wl = wl (tmp))
})


###################################################
### code chunk number 120: speed3
###################################################
tmp <- chondro [1 : 50]
system.time ({
  tmp.matrix <- tmp [[]]
  wl <- wl (tmp)
  for (i in seq_len (nrow (tmp)))
    tmp.matrix [i, ] <- interpolate (tmp.matrix [i, ], shifts [i], wl = wl)
  tmp [[]] <- tmp.matrix
})


###################################################
### code chunk number 121: tab-fn
###################################################
make.fn.table <- function (){
load ("functions.RData")
functions <- subset (functions, !internal)
functions$group <- functions$group[,drop=TRUE]

TeX.escape <- function (x){
#  x <- gsub ("^\\\\([^\\\\])", "\\\\\\\\\\1", x)
#  x <- gsub ("[^\\\\]\\\\$", "\\1\\\\\\\\", x)
  x <- gsub ("([^\\\\]|^)\\$", "\\1\\\\$", x)
  x <- gsub ("([^\\\\]|^)_", "\\1\\\\_", x)
  x <- gsub ("([^\\\\]|^)%", "\\1\\\\%", x)
  x
}

for (g in levels (functions$group)){
  cat ("\\multicolumn{2}{l}{\\emph{",g, "}}\\\\\n", sep = "")
  df <- t (functions [functions$group == g, c ("name", "description")])
  cat (paste (paste ("\\verb+", df[1,], "+", sep = ""), df[2,], sep = " & ", collapse ="\\\\\n"),"\\\\\n")
}
}
make.fn.table()


###################################################
### code chunk number 122: introduction.Rnw:1635-1640
###################################################
make.bib (c ("baseline", "compiler", "Rcpp", "inline"), file = "introduction-pkg.bib")
print (as.matrix(Sys.info()))
sessionInfo ()
rm (list = ls ())
library (tools)
#for (f in Sys.glob ("fig/*.pdf")) 
#  compactPDF ('$@', qpdf = '', gs_quality = 'screen', gs_extras = '-dDownsampleColorImages=false') 


