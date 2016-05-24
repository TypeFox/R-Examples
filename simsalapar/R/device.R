## Copyright (C) 2012 Marius Hofert and Martin Maechler
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 2 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


##' @title Cropping and Font Embedding PDF Device
##' @param file output file name (including extension .pdf)
##' @param crop crop command
##'        - NULL: crop with a (Unix-)suitable default commands;
##'        - "...": own crop command;
##'        - "": no cropping
##' @param embedFonts font embedding command
##'        - NULL: font embedding with a (Unix-)suitable default command;
##'        - "...": own font embedding command;
##'        - "": no font embedding
##' @param ... additional arguments passed to dev.off()
##' @return invisible
##' @author Marius Hofert
##' @note see Murrell, Ripley (2006),
##'       - http://stackoverflow.com/questions/4909214/issue-displaying-pdf-figures-created-with-r-on-ios-devices,
##'       - sfsmisc's pdf.end()
dev.off.pdf <- function(file="Rplots.pdf", crop=NULL, embedFonts="", ...)
{
    ## close device
    r <- dev.off(...)
    ## non-Unix case
    if(.Platform$OS.type != "unix") { # non-Unix
        iNcrop <- is.null(crop)
        iNeF <- is.null(embedFonts)
        if(iNcrop || iNeF)
            warning("'crop = NULL' and 'embedFonts = NULL' are only suitable for Unix")
        if(iNcrop) crop <- "" # continue without cropping
        if(iNeF) embedFonts <- "" # continue without font embedding
    }
    ## cropping
    f <- file.path(getwd(), file)
    if(is.null(crop)) { # crop with default command
        system(paste("pdfcrop --pdftexcmd pdftex", f, f, "1>/dev/null 2>&1"))
    } else if(nzchar(crop)) { # crop != "" crop with provided command
        system(crop)
    }
    ## font embedding
    if(is.null(embedFonts)) { # embed fonts with default command
        embedFonts(f, options="-dSubsetFonts=true -dEmbedAllFonts=true -dPDFSETTINGS=/printer -dUseCIEColor")
    } else if(nzchar(embedFonts)) { # embedFonts != "" embed fonts with provided command
        system(embedFonts)
    }
    ## return invisibly
    invisible(r)
}
