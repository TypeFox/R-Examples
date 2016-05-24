#' New Theme Elements
#' 
#' \code{ggtern} creates many new theme elements and inheritances, the following is an outline:
#'
#' Theme elements can inherit properties from other theme elements.
#' For example, \code{axis.title.x} inherits from \code{axis.title}, 
#' which in turn inherits from \code{text}. All text elements inherit
#' directly or indirectly from \code{text}; all lines inherit from
#' \code{line}, and all rectangular objects inherit from \code{rect}.
#'
#' Modifying the newly created items requires the same procedures as introduced in the ggplot2 \code{\link[ggplot2]{theme}} documentation.
#' Some convenience functions have been also newly created, proceed to \code{\link{convenience}} for additional information.
#' 
#' @section New/Additional Inheritance Structures:
#' \Sexpr[results=rd,stage=build]{ggtern:::rd_theme()}
#' **  \strong{NB:} \code{tern.panel.background}, whilst the ternary area is 'triangular' per-se, \code{\link{element_rect}} has been used, 
#' as it actually holds NO information regarding the geometry (width, height), only fill, color, 
#' size and linetype border (ie the style of how it will be rendered).
#' @aliases themeelements elements newelements theme-tern newthemes newtheme theme 
#' tern.panel.background tern.plot.background tern.axis.clockwise tern.axis 
#' tern.axis.hshift tern.axis.vshift tern.axis.line tern.axis.line.T tern.axis.line.L 
#' tern.axis.line.R tern.axis.text tern.axis.text.T tern.axis.text.L tern.axis.text.R 
#' tern.axis.arrow.show tern.axis.arrow.sep tern.axis.arrow.start tern.axis.arrow.finish 
#' tern.axis.arrow tern.axis.arrow.T tern.axis.arrow.L tern.axis.arrow.R tern.axis.arrow.text 
#' tern.axis.arrow.text.T tern.axis.arrow.text.L tern.axis.arrow.text.R tern.axis.title.show
#' tern.axis.title tern.axis.title.T tern.axis.title.L tern.axis.title.R tern.axis.ticks.length.major 
#' tern.axis.ticks.length.minor tern.axis.ticks.outside tern.axis.ticks.primary.show tern.axis.ticks.secondary.show 
#' tern.axis.ticks tern.axis.ticks.major tern.axis.ticks.major.T tern.axis.ticks.major.L tern.axis.ticks.major.R 
#' tern.axis.ticks.minor tern.axis.ticks.minor.T tern.axis.ticks.minor.L tern.axis.ticks.minor.R tern.panel.grid.major.show 
#' tern.panel.grid tern.panel.grid.major tern.panel.grid.major.T tern.panel.grid.major.L tern.panel.grid.major.R 
#' tern.panel.grid.minor.show tern.panel.grid.minor 
#' tern.panel.grid.minor.T tern.panel.grid.minor.L tern.panel.grid.minor.R ternary.options
#' panel.margin.tern tern.panel.expand tern.panel.rotate tern.panel.grid.ontop tern.panel.border.ontop
#' tern.axis.text.show
#' @name theme_elements
NULL


#SEARCH FOR THE ORIGINAL FUNCTIONS
#ggint$.element_tree     <- find_global_tern(".element_tree")
ggint$.element_tree.orig <- ggint$.element_tree #To determine the new set relative to the existing.
.el_def                  <- ggint$el_def

##TERNARY PANEL
ggint$.element_tree$tern.panel.background          = .el_def("element_rect", "panel.background",     description="Background of Ternary Clipping Area**")
ggint$.element_tree$tern.plot.background           = .el_def("element_rect", "plot.background",      description="Background of Ternary Plot Area**")

##AXIS ARROWS
ggint$.element_tree$tern.axis                      = .el_def("element_line", "line",                 description="Base line for ggtern object") #
ggint$.element_tree$tern.axis.hshift               = .el_def("numeric",                              description="Amount to nudge the plot horizontally") #
ggint$.element_tree$tern.axis.vshift               = .el_def("numeric",                              description="Amount to nudge the plot vertically") #
ggint$.element_tree$tern.axis.clockwise            = .el_def("logical",                              description="Clockwise or Anticlockwise Precession")

ggint$.element_tree$tern.axis.line                 = .el_def("element_line", "tern.axis",            description="Base Line") #
ggint$.element_tree$tern.axis.line.T               = .el_def("element_line", "tern.axis.line",       description="Line for TOP Axis") #
ggint$.element_tree$tern.axis.line.L               = .el_def("element_line", "tern.axis.line",       description="Line for LHS Axis") #
ggint$.element_tree$tern.axis.line.R               = .el_def("element_line", "tern.axis.line",       description="Line for RHS Axis") #
ggint$.element_tree$tern.axis.line.ontop           = .el_def("logical",                              description="Bring Axis Borders on Top of Everything")

#Axis Titles
ggint$.element_tree$tern.axis.title                = .el_def("element_text", "tern.axis.text",       description="Base Apex Title") #
ggint$.element_tree$tern.axis.title.T              = .el_def("element_text", "tern.axis.title",      description="Apex Title for TOP Axis") #
ggint$.element_tree$tern.axis.title.L              = .el_def("element_text", "tern.axis.title",      description="Apex Title for LHS Axis") #
ggint$.element_tree$tern.axis.title.R              = .el_def("element_text", "tern.axis.title",      description="Apex Title for RHS Axis") #
ggint$.element_tree$tern.axis.title.show           = .el_def("logical",                              description="Apex Titles Show or Hide")

#Axis Text
ggint$.element_tree$tern.axis.text                 = .el_def("element_text", "text",                 description="Base Text") #
ggint$.element_tree$tern.axis.text.T               = .el_def("element_text", "tern.axis.text",       description="Text for TOP Axis") #
ggint$.element_tree$tern.axis.text.L               = .el_def("element_text", "tern.axis.text",       description="Text for LHS Axis") #
ggint$.element_tree$tern.axis.text.R               = .el_def("element_text", "tern.axis.text",       description="Text for RHS Axis") #
ggint$.element_tree$tern.axis.text.show            = .el_def("logical",                              description="Axis Labels Show or Hide")

#Arrow
ggint$.element_tree$tern.axis.arrow                = .el_def("element_line", "tern.axis",            description="Base Arrow Line") #
ggint$.element_tree$tern.axis.arrow.T              = .el_def("element_line", "tern.axis.arrow",      description="Arrow Line for TOP Axis") #
ggint$.element_tree$tern.axis.arrow.L              = .el_def("element_line", "tern.axis.arrow",      description="Arrow Line for LHS Axis") #
ggint$.element_tree$tern.axis.arrow.R              = .el_def("element_line", "tern.axis.arrow",      description="Arrow Line for RHS Axis") #
ggint$.element_tree$tern.axis.arrow.text           = .el_def("element_text", "tern.axis.text",       description="Base Arrow Label") #
ggint$.element_tree$tern.axis.arrow.text.T         = .el_def("element_text", "tern.axis.arrow.text", description="Arrow Label on TOP Axis") #
ggint$.element_tree$tern.axis.arrow.text.L         = .el_def("element_text", "tern.axis.arrow.text", description="Arrow Label on LHS Axis") #
ggint$.element_tree$tern.axis.arrow.text.R         = .el_def("element_text", "tern.axis.arrow.text", description="Arrow Label on RHS Axis") #
ggint$.element_tree$tern.axis.arrow.show           = .el_def("logical",                              description="Arrows Show or Hide")
ggint$.element_tree$tern.axis.arrow.sep            = .el_def("numeric",                              description="Arrows Seperation from Axis")
ggint$.element_tree$tern.axis.arrow.start          = .el_def("numeric",                              description="Proportion of Axis when Arrow Starts")
ggint$.element_tree$tern.axis.arrow.finish         = .el_def("numeric",                              description="Proportion of Axis when Arrow Finishes")

#Ticks
ggint$.element_tree$tern.axis.ticks                = .el_def("element_line", "tern.axis",            description="Base Ticks") #
ggint$.element_tree$tern.axis.ticks.major          = .el_def("element_line", "tern.axis.ticks",      description="Base Major Ticks") #
ggint$.element_tree$tern.axis.ticks.major.T        = .el_def("element_line", "tern.axis.ticks.major",description="Base Major Ticks for TOP Axis") #
ggint$.element_tree$tern.axis.ticks.major.L        = .el_def("element_line", "tern.axis.ticks.major",description="Base Major Ticks for LHS Axis") #
ggint$.element_tree$tern.axis.ticks.major.R        = .el_def("element_line", "tern.axis.ticks.major",description="Base Major Ticks for RHS Axis") #
ggint$.element_tree$tern.axis.ticks.length.major   = .el_def("unit",                                 description="Ticks Major Ticklength")
ggint$.element_tree$tern.axis.ticks.length.minor   = .el_def("unit",                                 description="Ticks Minor Ticklength")
ggint$.element_tree$tern.axis.ticks.outside        = .el_def("logical",                              description="Ticks Outside or Inside") #
ggint$.element_tree$tern.axis.ticks.primary.show   = .el_def("logical",                              description="Ticks Show Primary")
ggint$.element_tree$tern.axis.ticks.secondary.show = .el_def("logical",                              description="Ticks Show Secondary")
ggint$.element_tree$tern.axis.ticks.minor          = .el_def("element_line", "tern.axis.ticks",      description="Base Minor Ticks") #
ggint$.element_tree$tern.axis.ticks.minor.T        = .el_def("element_line", "tern.axis.ticks.minor",description="Base Minor Ticks for TOP Axis") #
ggint$.element_tree$tern.axis.ticks.minor.L        = .el_def("element_line", "tern.axis.ticks.minor",description="Base Minor Ticks for LHS Axis") #
ggint$.element_tree$tern.axis.ticks.minor.R        = .el_def("element_line", "tern.axis.ticks.minor",description="Base Minor Ticks for RHS Axis") #

#Panel Grids
ggint$.element_tree$tern.panel.grid                = .el_def("element_line", "tern.axis",            description="Base Gridline") #
ggint$.element_tree$tern.panel.grid.major          = .el_def("element_line", "tern.panel.grid",      description="Base Major Gridline") #
ggint$.element_tree$tern.panel.grid.major.T        = .el_def("element_line", "tern.panel.grid.major",description="Major Gridline for TOP Axis") #
ggint$.element_tree$tern.panel.grid.major.L        = .el_def("element_line", "tern.panel.grid.major",description="Major Gridline for LHS Axis") #
ggint$.element_tree$tern.panel.grid.major.R        = .el_def("element_line", "tern.panel.grid.major",description="Major Gridline for RHS Axis") #
ggint$.element_tree$tern.panel.grid.major.show     = .el_def("logical",                              description="Show or Hide Major Gridline")
ggint$.element_tree$tern.panel.grid.minor          = .el_def("element_line", "tern.panel.grid",      description="Base Minor Gridline") #
ggint$.element_tree$tern.panel.grid.minor.T        = .el_def("element_line", "tern.panel.grid.minor",description="Minor Gridline for TOP Axis") #
ggint$.element_tree$tern.panel.grid.minor.L        = .el_def("element_line", "tern.panel.grid.minor",description="Minor Gridline for LHS Axis") #
ggint$.element_tree$tern.panel.grid.minor.R        = .el_def("element_line", "tern.panel.grid.minor",description="Minor Gridline for RHS Axis") #
ggint$.element_tree$tern.panel.grid.minor.show     = .el_def("logical",                              description="Show or Hide Minor Gridline")
ggint$.element_tree$tern.panel.grid.ontop          = .el_def("logical",                              description="Bring grids on top of everything else")
ggint$.element_tree$tern.panel.expand              = .el_def("numeric",                              description="The amount to expand the ternary plotting panel, in ratio to npc units")
ggint$.element_tree$tern.panel.rotate              = .el_def("numeric",                              description="The amount to rotate the ternary diagram in degrees")
