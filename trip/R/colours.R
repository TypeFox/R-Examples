# $Id: oc.theme.R 74 2013-03-21 15:28:34Z sluque $



#' SeaWiFS ocean colour colours
#'
#'
#' Generate ocean colour colours, using the SeaWiFS scheme
#'
#'
#' This is a high-contrast palette, log-scaled originally for ocean
#' chlorophyll.
#'
#' @aliases oc.theme oc.colors
#' @param x Number of colours to generate as part of a theme
#' @param n Number of colours to generate
#' @return A set of colours or a theme object.
#' @seealso
#'
#' Similar functions in sp \code{\link[sp]{sp.theme}},
#' \code{\link[sp]{bpy.colors}}
#' @keywords color
#' @examples
#' \dontrun{
#'  oc.colors(10)
#'  library(lattice)
#'  trellis.par.set(oc.theme())
#'  utils::example("trip-methods", package="trip",
#'                 ask=FALSE, echo=FALSE)
#'  tg <- tripGrid(tr)
#'  spplot(tg)
#' }
#' @export oc.theme
oc.theme <- function(x=50) list(regions=list(col=oc.colors(x)))

##' @rdname trip-internal
.oc.col <- c("#FFFFFF", "#A00AFF", "#A013FF", "#A01DFF", "#A027FF",
             "#A031FF", "#A03BFF", "#A045FF", "#A050FF", "#A05AFF",
             "#9C5FFF", "#9560FF", "#8F60FF", "#8860FF", "#8160FF",
             "#7A60FF", "#7460FF", "#6D60FF", "#6660FF", "#5F60FF",
             "#5560FF", "#4B60FF", "#4160FF", "#3760FF", "#2D60FF",
             "#2360FF", "#1960FF", "#0F60FF", "#0560FF", "#0061FF",
             "#0066FF", "#006BFF", "#0071FF", "#0076FF", "#007BFF",
             "#0080FF", "#0085FF", "#008AFF", "#008FFF", "#0094FF",
             "#0099FF", "#009EFF", "#00A3FF", "#00A8FF", "#00ADFF",
             "#00B2FF", "#00B7FF", "#00BCFF", "#03BFFF", "#08C0FF",
             "#0FC0FF", "#16C0FF", "#1DC0FF", "#23C0FF", "#2AC0FF",
             "#31C0FF", "#38C0FF", "#3EC0FF", "#45C0FF", "#4CC0FF",
             "#53C0FF", "#59C0FF", "#60C0FF", "#67C0FF", "#6EC0FF",
             "#74C0FF", "#7BC0FF", "#7CC0FC", "#78C1F7", "#73C3F2",
             "#6EC5EE", "#69C6E9", "#64C8E4", "#5FCADF", "#5ACBDA",
             "#55CDD5", "#50CFD0", "#4BD0CB", "#46D2C6", "#41D4C1",
             "#3CD6BC", "#36D7B6", "#31D9B1", "#2CDBAC", "#27DCA7",
             "#22DEA2", "#1FDE9C", "#1DDD94", "#1BDB8B", "#19D983",
             "#18D87B", "#16D672", "#14D46A", "#13D361", "#11D159",
             "#0FCF50", "#0ECE49", "#0CCC40", "#0ACA38", "#08C82F",
             "#07C727", "#05C51E", "#03C316", "#02C20D", "#00C005",
             "#01C100", "#05C300", "#0AC700", "#0FCA00", "#14CD00",
             "#19D100", "#1ED400", "#24D800", "#29DB00", "#2EDE00",
             "#33E100", "#38E500", "#3DE800", "#42EB00", "#47EF00",
             "#4CF200", "#51F500", "#56F800", "#5BFC00", "#60FE00",
             "#63FF00", "#67FF00", "#6AFF00", "#6DFF00", "#71FF00",
             "#74FF00", "#77FF00", "#7BFF00", "#7EFF00", "#81FF00",
             "#85FF00", "#88FF00", "#8CFF00", "#8FFF00", "#92FF00",
             "#96FF00", "#99FF00", "#9CFF00", "#A0FF00", "#A4FF00",
             "#A8FF00", "#ACFF00", "#B1FF00", "#B5FF00", "#B9FF00",
             "#BDFF00", "#C2FF00", "#C6FF00", "#C9FF00", "#CCFF00",
             "#CFFF00", "#D2FF00", "#D5FF00", "#D8FF00", "#DBFF00",
             "#DEFF00", "#E1FF00", "#E4FF00", "#E7FF00", "#EAFF00",
             "#ECFF00", "#EFFF00", "#F2FF00", "#F5FF00", "#F8FF00",
             "#FAFF00", "#FDFE00", "#FFFC00", "#FFF900", "#FFF600",
             "#FFF200", "#FFEF00", "#FFEC00", "#FFE900", "#FFE500",
             "#FFE200", "#FFDF00", "#FFDB00", "#FFD800", "#FFD500",
             "#FFD100", "#FFCE00", "#FFCA00", "#FFC700", "#FFC400",
             "#FFC000", "#FFBD00", "#FFBA00", "#FFB700", "#FFB300",
             "#FFB000", "#FFAC00", "#FFA900", "#FFA600", "#FFA200",
             "#FF9F00", "#FF9C00", "#FF9800", "#FF9500", "#FF9100",
             "#FF8E00", "#FF8B00", "#FF8700", "#FF8400", "#FF8100",
             "#FF7D00", "#FF7A00", "#FF7600", "#FF7300", "#FF7000",
             "#FF6C00", "#FF6900", "#FF6600", "#FF6200", "#FF5F00",
             "#FF5A00", "#FF5500", "#FF5000", "#FF4B00", "#FF4600",
             "#FF4100", "#FF3B00", "#FF3600", "#FF3100", "#FF2C00",
             "#FF2700", "#FF2200", "#FF1D00", "#FF1800", "#FF1300",
             "#FF0E00", "#FF0900", "#FF0400", "#FE0000", "#F90000",
             "#F40000", "#EF0000", "#EA0000", "#E50000", "#E00000",
             "#DB0000", "#D60000", "#D10000", "#CC0000", "#C70000",
             "#C20000", "#BD0000", "#B80000", "#B30000", "#AE0000",
             "#AA0000")

##' @rdname oc.theme
##' @export
oc.colors <- grDevices::colorRampPalette(.oc.col)



###_ + Emacs local variables
## Local variables:
## allout-layout: (+ : 0)
## End:
