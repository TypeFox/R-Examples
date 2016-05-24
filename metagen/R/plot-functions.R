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

#######################
### Colour Palettes ###
#######################

#' Colour palettes for colour blind people
#'
#' The palette with grey.
#'
#' This palette is directly taken from
#'
#'     http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
#'
#' Hence, I don't take any credit for this.
#' @examples
#' scale_fill_discrete <- function(...) scale_fill_manual(...,
#'   values=cbgPalette)
#' scale_colour_discrete <- function(...) scale_fill_manual(...,
#'   values=cbgPalette)
#' @export
cbgPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73"
                , "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#' Colour palettes for colour blind people
#'
#' The palette with black.
#'
#' This palette is directly taken from
#'
#'     http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
#'
#' Hence, I don't take any credit for this.
#' @examples
#' scale_fill_discrete <- function(...) scale_fill_manual(...,
#'   values=cbbPalette)
#' scale_colour_discrete <- function(...) scale_fill_manual(...,
#'   values=cbbPalette)
#' @export
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73"
                , "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##############################
### Rendering plot to disk ###
##############################

#' Render plot: To PDF
#'
#' Renders obj into a pdf-file of name: path++name.  Neat feature is
#' that the default size in A4.  Simply use the `scale` parameter to
#' adjust the size of the plot to a fraction of a page.
#'
#' @param name Should be self explanatory.
#' @param plotObj Should be self explanatory.
#' @param path Should be self explanatory.
#' @param scale Should be self explanatory.
#' @param height Should be self explanatory.
#' @param width Should be self explanatory.
#' @export
render <- function(  name, plotObj, path, scale=1
                   , height=11.6, width=8.2) {
    height <- height * scale
    message("Creating pdf-figure... ", name)
    pdf(paste(path, name, ".pdf", sep=""), height=height, width=width)
    print(plotObj)
    dev.off()
}

#' Render plot: To SVG
#'
#' Renders obj into a svg-file of name: path++name.  Neat feature is
#' that the default size in A4.  Simply use the `scale` parameter to
#' adjust the size of the plot to a fraction of a page.
#'
#'
#' @param name Should be self explanatory.
#' @param plotObj Should be self explanatory.
#' @param path Should be self explanatory.
#' @param scale Should be self explanatory.
#' @param height Should be self explanatory.
#' @param width Should be self explanatory.
#' @export
renderSVG <- function(  name, plotObj, path, scale=1
                      , height=11.6, width=8.2) {
    height <- height * scale
    message("Creating svg-figure... ", name)
    svg(paste(path, name, ".svg", sep=""), height=height, width=width)
    print(plotObj)
    dev.off()
}
