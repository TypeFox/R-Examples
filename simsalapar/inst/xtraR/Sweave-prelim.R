## Copyright (C) 2012, 2014 Marius Hofert and Martin Maechler
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

if(FALSE) ## "load me" via
source(system.file("xtraR", "Sweave-prelim.R", package = "simsalapar", mustWork=TRUE))


## Custom graphics device (for cropping .pdf):
pdfCrop <- function(name, width, height, ...) {
    f <- paste0(name, ".pdf")
    grDevices::pdf(f, width=width, height=height, onefile=FALSE)
    assign(".pdfCrop.name", f, envir=globalenv())
}
pdfCrop.off <- function() { # used automagically
    grDevices::dev.off() # closing the pdf device
    f <- get(".pdfCrop.name", envir=globalenv())
    system(paste("pdfcrop --pdftexcmd pdftex", f, f, "1>/dev/null 2>&1"),
           intern=FALSE) # crop the file (relies on PATH)
}

if(!exists("JSS")) JSS <- FALSE # just in case
if(JSS) { # JSS style, but ugly
options(width=70, useFancyQuotes=FALSE, prompt="R> ", continue="+  ")
} else {
options(width=70, useFancyQuotes=FALSE, prompt="> ", continue="  ")
}
## as we use pdfCrop(): unneeded:
## options(SweaveHooks=list(fig=function() par(mar=c(4, 4, 0.4, 0.7))))

## patchDVI:::patchSynctex setup
.TexRoot <- "parallel.tex"

## "Global"  system.time() saving etc:

if(getRversion() < "3.2") dir.exists <- function(x)
  is.character(x) && file.exists(x) && file.info(path.expand(x))$isdir
## for back compatibility:
.dir.exists <- function(x) { .Deprecated("dir.exists") ; dir.exists(x) }

##' Create directory if needed
mkDir <- function(name) if(!dir.exists(name)) dir.create(name, recursive=TRUE)

##' Create "unique" file name - dependent on current time
##' @examples mkFname("safe", time="min", ext="rds")
mkFname <- function(stem, ext, sep="_", time = c("none", "hr", "min", "sec")) {
    stopifnot(is.character(stem))
    time <- match.arg(time)
    tim <- Sys.time()
    timCh <- switch(time, ##-> ?format.POSIXct
                    "none" = format(as.Date(tim)),
                    "hr"  = format(tim, "%Y-%m-%d_%H"),
                    "min" = format(tim, "%Y-%m-%d_%H:%M"),
                    "sec" = format(tim, "%Y-%m-%d_%H:%M:%S"))
    fn <- paste(stem, timCh, sep=sep)
    if(nzchar(ext)) paste(fn, ext, sep=".") else fn
}

##' Create "unique" file name inside note-specific directory
mkD.Fname <- function(dir, stem, ext="", first = Sys.info()[["nodename"]]) {
    stopifnot(is.character(dir), is.character(stem), is.character(first))
    mkDir(dir <- file.path(dir, first))
    file.path(dir, mkFname(stem, ext=ext))
}
