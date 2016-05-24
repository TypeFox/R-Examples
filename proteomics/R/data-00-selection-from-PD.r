# Copyright (C) 2012-2014 Thomas W. D. Möbius (kontakt@thomasmoebius.de)
#
#     This program is free software: you can redistribute it and/or
#     modify it under the terms of the GNU General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see
#     <http://www.gnu.org/licenses/>.

################################
# Methods for data parsing     #
#   when the data come from PD #
################################

listModifications <- function(mods) {
    a <- as.character(levels(mods))
    b <- strsplit(a, ";? ")
    return(extractModifications(unique(unlist(b))))
}

extractModifications <- function (modstring) {
    mod  <- unlist(strsplit(as.character(modstring), ";? "))
    pos.from <- regexpr('[0-9]+\\([a-zA-Z]', mod)
    pos.till <- pos.from +  attributes(pos.from)$match.length
    position <- as.numeric(substr(mod, pos.from, pos.till - 3))
    mod.from <- regexpr('\\([a-zA-Z0-9]+\\)', mod)
    mod.till <- mod.from + attributes(mod.from)$match.length
    modifica <- tolower(substr(mod, mod.from + 1, mod.till - 2))
    code <- mapvalues(modifica, unique(modifica), substr(unique(modifica), 1,1))
    mat <- data.frame(modstring, pos=position, mod=modifica, code=code)
    mat <- mat[order(mat$pos, decreasing=T),]
    return(mat)
}

attachModifications <- function (peptide, modmat) {
    peptide <- toupper(as.character(peptide))
    modmat <- na.omit(modmat[order(modmat$pos, decreasing=T),])
    if (dim(modmat)[1] > 0) {
        for (i in 1:dim(modmat)[1]) {
            peptide <- sub(sprintf('(?<=.{%d})', modmat$pos[i]-1), modmat$code[i], peptide, perl=TRUE)
        }
    }
    return(peptide)
}

toPeptide <- function(sec, modstring) {
    return(attachModifications(sec, extractModifications(modstring)))
}

#' Data parsing -- from Proteom Discover v1.4
#'
#' Has been tested with PD v1.4
#'
#' This is a rather neat function that allows to get data from an export form
#' the software Proteom Discoverer into R and parsed into a reasonable data
#' frame such one can work with it. It will also add a few statistics and
#' create unique identifiers for all identified peptides. You may argue that
#' this functionality alone is worth the import of the whole package.
#'
#' @param dwide raw data from a PD export.
#' @param ch the column names which hold the reporter ion intensities.
#' @param ref the colmun name which holds the reporter ion intesities of the
#' reference channel.
#' @examples
#' \dontrun{
#' bio1 <- read.csv("my-proteome-discoverer-v1.4-export-experiment-1.csv")
#' bio2 <- read.csv("my-proteome-discoverer-v1.4-export-experiment-2.csv")
#' run1 <- droplevels(bio1[bio1$Quan.Usage == "Used",])
#' run2 <- droplevels(bio2[bio2$Quan.Usage == "Used",])
#' channels <- c("X113", "X114", "X115", "X116", "X117", "X118", "X119", "X121")
#' reference <- c("X121")
#'
#' run1 <- meetSelection(run1, channels, reference)
#' run2 <- meetSelection(run2, channels, reference)
#'
#' run1$experiment <- factor(1, levels=1:2, labels=c("iTRAQ-1", "iTRAQ-2"))
#' run2$experiment <- factor(2, levels=1:2, labels=c("iTRAQ-1", "iTRAQ-2"))
#' runs <- rbind(run1, run2)
#' }
#' @export
meetSelection <- function (dwide, ch, ref) {
    if(!(is.na(ref)  || (ref %in% ch))) {ch <- union(ch, ref)}
    dwide <- dwide[c(  ch
                 , "Intensity"
                 , "Sequence"
                 , "Modifications"
                 , "Protein.Group.Accessions"
                 , "QuanResultID"
                 , "RT..min."
                 , "Charge"
                 , "X.M..ppm."
                 , "Isolation.Interference...."
                 , "X..Missed.Cleavages"
                 , "Spectrum.File")]
    dwide$Sequence <- factor(toupper(as.character(dwide$Sequence)))
    dwide$QuanResultID <- factor(dwide$QuanResultID)
    dwide <- rename(dwide, replace=c(  Intensity="intensity"
                                     , Sequence="sequence"
                                     , Protein.Group.Accessions="protein"
                                     , QuanResultID="id"
                                     , RT..min.="retention"
                                     , Modifications="mods"
                                     , Charge="charge"
                                     , X.M..ppm.="ppm"
                                     , Isolation.Interference....="interference"
                                     , X..Missed.Cleavages="missed.cleavages"
                                     , Spectrum.File="injection"
                                     ))
    dwide <- droplevels(dwide)
    dwide$peptide <- mapply(toPeptide, dwide$sequence, dwide$mods)
    dwide$peptide <- as.factor(dwide$peptide)
    dwide$charge  <- factor(dwide$charge, ordered=TRUE)
    dwide$mods <- NULL

    dwide$missing <- apply(dwide[ch], 1, function(i) any(is.na(i)))

    attr(dwide, "channelnames") <- ch
    attr(dwide, "channels") <- match(ch, names(dwide))
    attr(dwide, "referencename") <- ref
    attr(dwide, "reference") <- match(ref, names(dwide))
    return(dwide)
}
