### R code from vignette source 'laser.Rnw'
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
check.req.pkg ("rgl", donothing = "")


###################################################
### code chunk number 6: rawspc
###################################################
getOption("SweaveHooks")[["fig"]]()
laser <- scan.txt.Renishaw ("rawdata/laser.txt.gz", data = "ts")
plot (laser, "spcprctl5") 


###################################################
### code chunk number 7: cut
###################################################
getOption("SweaveHooks")[["fig"]]()
laser <- laser [,,-75~0]
plot (laser, "spcprctl5") 


###################################################
### code chunk number 8: wlspc1
###################################################
wl (laser) <- wl (laser) + 50


###################################################
### code chunk number 9: wlcalc
###################################################
wl (laser) <- list (
   wl = 1e7 / (1/405e-7 - wl (laser)),
   label = expression (lambda / nm)
   )


###################################################
### code chunk number 10: wlspc
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (laser, "spcprctl5") 


###################################################
### code chunk number 11: save
###################################################
save (laser, file = "laser.rda") 
laser


###################################################
### code chunk number 12: locator (eval = FALSE)
###################################################
## wls <- locator()$x


###################################################
### code chunk number 13: laser.Rnw:119-120
###################################################
wls <-  c(405.0063, 405.1121, 405.2885, 405.3591)


###################################################
### code chunk number 14: markspc
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (laser, "spcmeansd")
cols <- c("black", "blue", "red", "darkgreen")
abline (v = wls, col = cols )  


###################################################
### code chunk number 15: ts
###################################################
getOption("SweaveHooks")[["fig"]]()
plotc (laser [,, wls], spc ~ t, groups = .wavelength, type = "b",  cex = 0.3, col = cols)


###################################################
### code chunk number 16: tsextra
###################################################
getOption("SweaveHooks")[["fig"]]()
plotc (laser [,, wls], spc ~ t | .wavelength, type = "b", cex = 0.3, col = "black")


###################################################
### code chunk number 17: plotmatr
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (laser, "mat", contour = TRUE, col = "#00000060")


###################################################
### code chunk number 18: plotmatt
###################################################
getOption("SweaveHooks")[["fig"]]()
levelplot (spc ~ .wavelength * t, laser, contour = TRUE, col ="#00000080")


###################################################
### code chunk number 19: libraryrgl
###################################################
require (rgl)


###################################################
### code chunk number 20: rgl-plot (eval = FALSE)
###################################################
## laser <- laser [,,404.8 ~ 405.6] / 10000
## laser$t <- laser$t / 3600
## cols <- rep (matlab.palette (nrow (laser)), nwl (laser))
## 
## surface3d(y = wl(laser), x = laser$t, z = laser$spc, col =  cols)
## surface3d(y = wl(laser), x = laser$t, z = laser$spc + .1 * min (laser), 
##           col =  "black", alpha = .2,front = "lines", line_antialias = TRUE)
## 
## aspect3d (c(1, 1, 0.25))
## 
## axes3d (c ('x+-', 'y--', 'z--'))
## axes3d ('y--', nticks = 25, labels= FALSE)
## mtext3d ("t / h",    'x+-', line = 2.5)
## mtext3d ("lambda / nm",   'y--', line = 2.5)
## mtext3d ("I / a.u.", 'z--', line = 2.5)


###################################################
### code chunk number 21: rgl-do
###################################################
if (exists (".rgl") && .rgl){        # use extra argument to turn on rgl only on local machine as it doesn't work on r-forge.
  if (require (rgl)){
    open3d (windowRect=c(20,20,600, 350))  # this is needed only for automatically 
                                        # producing the snapshot
laser <- laser [,,404.8 ~ 405.6] / 10000
laser$t <- laser$t / 3600
cols <- rep (matlab.palette (nrow (laser)), nwl (laser))

surface3d(y = wl(laser), x = laser$t, z = laser$spc, col =  cols)
surface3d(y = wl(laser), x = laser$t, z = laser$spc + .1 * min (laser), 
          col =  "black", alpha = .2,front = "lines", line_antialias = TRUE)

aspect3d (c(1, 1, 0.25))

axes3d (c ('x+-', 'y--', 'z--'))
axes3d ('y--', nticks = 25, labels= FALSE)
mtext3d ("t / h",    'x+-', line = 2.5)
mtext3d ("lambda / nm",   'y--', line = 2.5)
mtext3d ("I / a.u.", 'z--', line = 2.5)
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
### code chunk number 22: makebib
###################################################
make.bib (file = "laser-pkg.bib")


###################################################
### code chunk number 23: laser.Rnw:241-243
###################################################
make.bib (file = "laser-pkg.bib") 
sessionInfo ()
rm (list = ls ())
library (tools)


