#' @title Construct default structure legend to plot on a \code{strucplot}
#' display.
#'
#' @description A structure legend is always printed on the console. It can also
#' be optionally added to the trellis plot. This non-exported function
#' constructs a default legend for this option.
#' 
#' @return A text grob that can be included as part of the \code{strucplot}
#' trellis plot.
#' 
#' @param struc The "structure" attribute of a \code{strucplot} object
#' 
#' @param legendLoc One of c("top","right","bottom","left","newpage")
#'  The first 4 specify on which side of the trellis display to place the legend.
#'  The last indicates that the legend will be plotted at the center of a
#'  separate page.
#' 
#' @param lbl A character string title for the legend. 
#'   Default = \code{"PLOT STRUCTURE"} 
#' 
#' @param heading A character vector of length 2 giving the headings for the 
#'   horizontal and vertical conditioning variables portions of the legend.
#'   Default = \code{c("Horizontal", "Vertical")}
#' 
#' @param miss A character string to use when a conditioning variable is absent.
#' Default = \code{"None"}
#' 
#' @param cex.font  Multiplier for text font size. Default = 1
#'
#' @param cex.lbl  Multiplier for label font size. Default = 1.25
#' 
#' @param col Text color for text. Default = \code{"black"}
#' 
#' @param abbrevLength  \code{minlength} argument for \code{\link{abbreviate}} 
#'   function to abbreviate names of the conditioning factors. Default = 0 = no 
#'   abbreviation.
#'  
#' @param ... Additional arguments to be used in a \code{gp} list for
#'  controlling text appearance. See \code{\link[grid]{gpar}} for possibilities.
#'  
#' @seealso \code{\link{print.structured}}
#'  
defaultStrucLegend <- 
  function(struc
   ,legendLoc = c("top","right","bottom","left","newpage")
   ,lbl = "PLOT STRUCTURE"  ## legend label
   ,heading = c("Horizontal", "Vertical") ## Heading for dimensions
   ,miss = "None" ## text if there is no conditioning in a dimension
   ,cex.font = 1 ## multiplier for text font size
   ,cex.lbl = 1.25 ## multiplier for label font size
   ,col = "black" ## text color for text
   ,abbrevLength = 0 ## minlength argument for abbreviate function to abbreviate
  ## names of the conditioning factors. 0 = no abbreviation is the default
   ,... ## Additional arguments for controlling font appearance to be used as
  ## part of gp argument
)
{
  ## Create vector of text labels for horiz and vertical directions
  ## First check the arguments
  if(!(is.character(miss) && length(miss)==1)) stop(
    "miss must be a single character string")
  if(!(is.character(lbl) && length(lbl)==1)) stop(
    "lbl must be a single character string")
  if(!(is.character(heading) && length(heading)==2)) stop(
    "heading must be a length 2 character vector")  
  dots <- list(...)
  if(is.null(dots[["col"]]))dots[["col"]] <- col
  gp <- do.call(gpar,dots)
  ## create text display
  txt <- displayStruc(struc= struc
    ,miss = miss
    ,heading = heading
    ,abbrevLength = abbrevLength)
  
  ## construct grobs
  titl <- textGrob(lbl,gp = gpar(cex = cex.lbl))
  hz <- textGrob(txt[[1]], gp = gpar(cex=cex.font), x = 0, y=1,
                 just = c("left","top"))
  vrt <- textGrob(txt[[2]],gp = gpar(cex=cex.font), x = 0, y=1,
                  just = c("left", "top"))
  parts <- list(titl,hz,vrt)
  grbwid <- lapply(parts,grobWidth)
  grbht <- lapply(parts,grobHeight)
  spcht <- unit(1,"line")
  spcwid <- unit(1,"char")
  if(legendLoc %in% c("right","left","newpage")){
    owid <- do.call(max,grbwid)
    midwid <- do.call(max,grbwid[2:3])
    sidewid <- (owid - midwid )*.5
    pad <- (if(legendLoc == "right")c(.5,0) else c(0,.5))*spcwid
    glay <- grid.layout(nrow = 5, ncol = 5,
              widths = unit.c(pad[1],sidewid, midwid,sidewid,pad[2]),
              heights = unit.c(grbht[[1]],spcht,grbht[[2]],spcht,grbht[[3]])
            )
    frm <- frameGrob(layout = glay, gp = gp )
    frm <- placeGrob(frm, parts[[1]],row = 1,col = 2:4)
    frm <- placeGrob(frm,parts[[2]], row = 3, col = 3)
    frm <- placeGrob(frm,parts[[3]], row = 5, col = 3)
  } else{
    leftwid <- grbwid[[2]] + spcwid
    w <- max(leftwid + grbwid[[3]],grbwid[[1]])
    pad <- (if(legendLoc == "top")c(0,.5) else c(.5,0)) * spcht
    glay <- grid.layout(nrow = 5, ncol = 3,
                    widths = unit.c(grbwid[[2]],spcwid, w - leftwid),
                    heights = unit.c(pad[1],grbht[[1]],spcht, 
                            do.call(max,grbht[2:3]),pad[2])
              )
    frm <- frameGrob(layout = glay, gp = gp)
    frm <- placeGrob(frm, parts[[1]],row = 2)
    frm <- placeGrob(frm, parts[[2]],row = 4, col = 1)
    frm <- placeGrob(frm, parts[[3]],row = 4, col = 3)
  }
  frm
}

#' @describeIn defaultStrucLegend An auxiliary function that constructs the
#' character strings for the conditioning variables.
displayStruc <- function(struc
                         ,miss = "No Conditioning"
                         ,heading = c("Horizontal","Vertical")
                         ,abbrevLength = 0)
  
  ## Helper function for constructing structure legends  
{
  txt <- lapply(struc,function(u){
    if(!length(u)) miss
    else 
      paste(if(abbrevLength)abbreviate(names(u), minlength = abbrevLength) 
            else names(u),sapply(u, paste, collapse = "  "), sep = ":  ")
  })
  sprintf("_%s_\n%s",heading,sapply(txt,paste,collapse = "\n"))
}