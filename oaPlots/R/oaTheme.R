
#' Apply OA ggplot2 theme
#' @param p ggplot2 plot object
#' @param useOAColors boolean which indicates wether or not to use the oaColors 
#' package to provide a color scheme. Default: TRUE
#' @param expand specify wether or not to expand the axis valid options are: 
#' (both, x, y, none) Default: both
#' @param bgColor specify a different background color 
#' (useful for plotting colors with alpha values) Default: gray(0.9)
#' @return ggplot2 plot object
#' @author Willem Ligtenberg
#' @importFrom oaColors oaPalette
#' @importFrom grid unit
#' @importFrom ggplot2 theme_update
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 rel
#' @export
oaTheme <- function(p, useOAColors = TRUE, expand = "both", bgColor = gray(0.9)){
  #Apply basic theme
  theme_update(panel.grid.minor = element_blank())
  theme_update(axis.ticks = element_blank())
  theme_update(panel.grid.major = element_line(size = 0.6, color = "white"))
  theme_update(panel.margin = unit(1, "lines")) # changes space between panels
  theme_update(axis.text.x = element_text(hjust = 0, vjust = 1, color = gray(0.5), size = rel(0.85))) # nice to have them a bit higher
  theme_update(axis.text.y = element_text(hjust = 1, vjust = 0, color = gray(0.5), size = rel(0.85)))
  theme_update(axis.title.x = element_text(vjust = -0.5, color = gray(0.3), size = rel(1)))
  theme_update(axis.title.y = element_text(angle = 90, vjust = 0.3, color = gray(0.3), size = rel(1)))
  theme_update(axis.ticks.length = unit(0.2, "lines"))
  theme_update(axis.ticks.margin = unit(0, "lines"))
  theme_update(plot.title = element_text(color = gray(0.2), size = rel(1.5), vjust = 1.5))
  theme_update(legend.title = element_text(color = gray(0.2), size = rel(1)))
  theme_update(legend.text = element_text(color = gray(0.2), size = rel(0.9)))
  theme_update(strip.text = element_text(color = gray(0.2), size = rel(0.9)))
  #theme_update(strip.background = element_rect())
  #theme_update(plot.margin = unit(c(1, 1, 1, 1), "lines"))
  theme_update(panel.border = element_rect(fill = NA, size = 0, color = "white")) # seems like there is a bug, use white for now
  theme_update(panel.background = element_rect(fill = bgColor, size = 0))
  
  if(expand %in% c("both", "x")){
    p <- p + scale_x_continuous(expand=c(0.1, 0.1)) # expand the x-axis
  }
  if(expand %in% c("both", "y")){
    p <- p + scale_y_continuous(expand=c(0.1, 0.1)) # expand the x-axis
  }
  
  #Apply oaColors
  #For each layer get the colour/fill variable
  colorVars <- c()
  fillVars <- c()
  
  for(layer in p$layers){
    colorVars <- c(colorVars, as.character(layer$mapping[["colour"]]))
    fillVars <- c(fillVars, as.character(layer$mapping[["fill"]]))
  }
  uniqueColorVars <- unique(colorVars)
  uniqueFillVars <- unique(fillVars)
  
  #Process color
  totalColorLevels <- 0
  for(colorVar in uniqueColorVars){
    if(class(p$data[, colorVar]) %in% c("factor", "character")){
      totalColorLevels <- totalColorLevels + length(unique(p$data[, colorVar]))
    }
  }
  
  if(totalColorLevels > 0 & totalColorLevels < 10 & useOAColors){
    p <- p + scale_colour_manual(values=unname(oaPalette(totalColorLevels)))
  }
  
  #Process fill
  totalFillLevels <- 0
  for(fillVar in uniqueFillVars){
    if(class(p$data[, fillVar]) %in% c("factor", "character")){
      totalFillLevels <- totalFillLevels + length(unique(p$data[, fillVar]))
    }
  }
  
  if(totalFillLevels > 0 & totalFillLevels < 10 & useOAColors){
    p <- p + scale_fill_manual(values=unname(oaPalette(totalFillLevels, alpha = 0.09)))
  }
  return(p)
}
