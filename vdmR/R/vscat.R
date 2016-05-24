#' Generate scatter plot with interactive functions
#'
#' \code{vscat} generates scatter plot of variable \code{x} and \code{y} of data frame \code{data} with interactive functions.
#'
#' @docType methods
#' @param x,y column name of data frame \code{data} for x-axis and y-axis of scatter plot
#' @param data data frame for default data set
#' @param name character for the name of the generated scatter plot
#' @param tag character for the common name of a series of linked plots
#' @param ... aesthetic mappings to be passed to ggplot2 methods
#' @importFrom grDevices dev.off pdf
#' @importFrom stats asOneSidedFormula
#' @importFrom utils packageVersion
#' @export
#' @examples
#' data(vsfuk2012)
#' vscat(MarriageRate, DivorceRate, vsfuk2012, "scat1", "vsfuk2012", colour=Type)
#' vlaunch(vsfuk2012, "main", "vsfuk2012", browse=FALSE)
#'

vscat <- function(x, y, data, name, tag, ...){

  argnames <- names(as.list(match.call(expand.dots = FALSE)[-1]))
  arguments <- as.list(match.call()[-1])
  aesthetics <- plyr::compact(arguments[allaes])
  aesthetics <- aesthetics[!is.const(aesthetics)]
  aes_names <- names(aesthetics)
  aesthetics <- rename.aes(aesthetics)
  class(aesthetics) <- "uneval"

  params <- arguments[setdiff(names(arguments), c(aes_names,argnames))]
  params <- lapply(params, eval)
  class(params)   <- "uneval"

  jspath <- file.path(system.file(package="vdmR"), "exec/vdmr_scat.js")
  file.copy(jspath, paste(name, ".", tag, ".js", sep=""), overwrite=TRUE)
  pdf(file=NULL, width=7, height=5)

  scat <- ggplot2::ggplot(data, aesthetics)
  if(packageVersion("ggplot2")>'1.0.1'){
    params$na.rm <- TRUE
    scat <- scat + ggplot2::layer(geom="point", position="identity",
                                  stat="identity", params=params)
  } else {
    scat <- scat + ggplot2::layer(geom="point", geom_params=params)
  }


  scatgrob <- ggplot2::ggplotGrob(scat)

  grid::grid.newpage()
  grid::grid.draw(scatgrob)
  grid::grid.force()
  grid::grid.gedit("geom_point.points", name="geom_point.points")
  gridSVG::grid.script(file=paste0(name, ".", tag,".js"))
  gridSVG::grid.script(paste0("var winname= '", name,"';"))
  gridSVG::grid.script(paste0("var x= ", rjson::toJSON(scat$data[,as.character(scat$mapping$x)]), ";"))
  gridSVG::grid.script(paste0("var y= ", rjson::toJSON(scat$data[,as.character(scat$mapping$y)]), ";"))

  svgfn <- paste0(name, ".", tag, ".svg")

  gridSVG::grid.export(svgfn, htmlWrapper=TRUE, exportMappings="file",
                       xmldecl="<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")

  invisible(dev.off())


}
