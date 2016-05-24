#' Generate prallel coordinate plot with interactive functions
#' 
#' \code{vpcp} generates parallel coordinate plot of specific data columns of data frame \code{data} with interactive functions.
#' 
#' @docType methods
#' @param data data frame for default data set
#' @param columns a vector of variables (either names or indices) to be axes in the plot
#' @param name character for the name of the generated scatter plot
#' @param tag character for the common name of a series of linked plots
#' @param groupColumn a single variable to group (color) by
#' @param scale method used to scale the variables
#' @param alphaLines value of alpha scaler for the lines of the parcoord plot or a column name of the data
#' @param missing method used to handle missing values
#' @importFrom grDevices dev.off pdf
#' @export
#' @examples
#' data(vsfuk2012)
#' vpcp(vsfuk2012, 4:17, "pcp1", "vsfuk2012", scale="uniminmax")
#' vlaunch(vsfuk2012, "main", "vsfuk2012", browse=FALSE)
#' 

vpcp <- function(data, columns, name, tag,
                 groupColumn=NULL, scale="std", alphaLines=0.5, missing="exclude"){

  jspath <- file.path(system.file(package="vdmR"), "exec/vdmr_pcp.js")
  file.copy(jspath, paste(name, ".", tag, ".js", sep=""), overwrite=TRUE)
  
  pdf(file=NULL, width=7, height=5)
  pcp <- GGally::ggparcoord(data, groupColumn=groupColumn, columns=columns, alphaLines=alphaLines, scale=scale, missing=missing) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  pcpgrob <- ggplot2::ggplotGrob(pcp)
  grid::grid.newpage()
  grid::grid.draw(pcpgrob)
  grid::grid.force()
  grid::grid.gedit("GRID.polyline", name="GRID.polyline")
  
  gridSVG::grid.script(paste("var ncol=", length(columns), ";"))
  gridSVG::grid.script(paste("var nrow=", nrow(data), ";"))
  gridSVG::grid.script(file=paste(name, ".", tag,".js", sep=""))
  gridSVG::grid.script(paste("var winname= '", name, "';", sep=""))
  
  svgfn <- paste0(name, ".", tag, ".svg")
  
  gridSVG::grid.export(svgfn, htmlWrapper=TRUE, exportMappings="file",
                       xmldecl="<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")

  invisible(dev.off())
  
}
