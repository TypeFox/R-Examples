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

##################################################################################
## Merging the summerised iTRAQ results from multiple files into one data frame ##
##################################################################################

barcoding <- function (id, dlong, sampledesign) {
    samp <- droplevels(sampledesign[sampledesign$experiment == id,])
    dlong <- droplevels(dlong)
    dlong$barcode <- mapvalues(dlong$variable, from=paste(samp$channel), to=paste(samp$barcode))
    return(dlong)
}

#' Merging multiple experiments
#'
#' At the end each channel in each iTRAQ experiment can be uniquly identified
#' by a barcode.  If two channels of different experiments correspond to the
#' same subject, the same barcode may be used and a method of combinig these
#' measurments be chosen.
#'
#' @param files data frame of file names and corresponding ids.
#' @param path leading to the files
#' @param sampledesign data frame of ids, channelnames and corresponding barcodes.
#'
#' @export
mergeFrames <- function (files, path, sampledesign) {
    frames <- function (id) {
        dl <- readRDS(paste0(path,files$name[files$id == id],".rds"))
        return(barcoding(id, dl, sampledesign))
    }

    return(do.call(rbind,lapply(files$id, frames)))
}
