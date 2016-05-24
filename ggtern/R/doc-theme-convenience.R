#' Theme Convenience Functions
#' 
#' @description
#' \code{ggtern} has made available a number of convenience functions for rapid tweaking of the various theme elements, 
#' for a full list of the available theme elements which can be manually modified, see \link[=themeelements]{HERE}.
#' 
#' @section Convenience Functions:
#' Convenience functions that ship with \code{ggtern}, to assist in the rapid modification of key theme elements:
#' \itemize{
#'   \item \code{\link[=theme_showtitles]{Show/Hide Axis Titles}}
#'   \item \code{\link[=theme_showarrows]{Show/Hide Arrows}} 
#'   \item \code{\link[=theme_showgrid]{Show/Hide Grids}}
#'   \item \code{\link[=theme_showprimary]{Show/Hide Primary/Secondary Ticks}}
#'   \item \code{\link[=theme_showlabels]{Show/Hide Axis Ticklabels}}
#'   \item \code{\link[=theme_clockwise]{Clockwise/Anticlockwise Axis Precession}}  
#'   \item \code{\link[=theme_ticksoutside]{Ticks Inside or Outside of the Main Plot Area}}
#'   \item \code{\link[=atomic_percent]{Atomic or Weight Percent Arrow Label Suffix.}}
#'   \item \code{\link[=theme_rotate]{Rotate the plot by X degrees or radians}}
#' }
#' 
#' @section Manual Modification:
#' For manual modification on a per-element basis:
#' \itemize{
#'   \item \code{\link[=theme_elements]{Ternary Theme Elements}} 
#' }
#' 
#' @section Default Themes:
#' Default (complete) themes which ship with \code{ggtern}:
#' \itemize{
#'   \item \code{\link[=theme_complete]{Complete Themes}} 
#' }
#' 
#' @aliases theme_convenience convenience convenience_functions
#' @name convenience_functions
#' @rdname convenience_functions
#' @examples
#' 
#' #Load data and create the base plot.
#' plot <- ggtern(data=data.frame(x=1,y=1,z=1),aes(x,y,z)) + geom_point() + theme_bw()
#' plot
#' 
#' #Show or Hide Arrows
#' last_plot() + theme_showarrows()
#' #last_plot() + theme_hidearrows()
#' 
#' #Major/Minor Grids?
#' last_plot() + theme_nogrid_minor()
#' #last_plot() + theme_nogrid_major()
#' #last_plot() + theme_showgrid()
#' #last_plot() + theme_nogrid()
#' 
#' #Clockwise/Anticlockwise Precession
#' last_plot() + theme_clockwise()
#' #last_plot() + theme_anticlockwise()
#' 
#' #Ticks Inside or Outside
#' last_plot() + theme_ticksoutside()
#' #last_plot() + theme_ticksinside()
#' 
#' #Show/Hide BOTH Primary and Secondary Ticks
#' last_plot() + theme_showticks()
#' #last_plot() + theme_hideticks()
#' 
#' #Show/Hide EITHER Primary OR Secondary Ticks.
#' last_plot() + theme_showprimary() + theme_hidesecondary()
#' #last_plot() + theme_hideprimary() + theme_showsecondary()
#' 
#' #Atomic / Weight Percent
#' last_plot() + theme_showarrows() + atomic_percent()
#' #last_plot() + theme_showarrows() + weight_percent()
#' #last_plot() + theme_showarrows() + custom_percent("Atomic Percent")
#' 
#' #Rotation
#' last_plot() + theme_rotate(30)
NULL