#' Generate histogram with interactive functions
#'
#' \code{vscat} generates histogram of variable \code{x} of data frame \code{data} with interactive functions.
#'
#' @docType methods
#' @param x column name of data frame \code{data} for drawing histogram
#' @param data data frame for default data set
#' @param name character for the name of the generated histogram
#' @param tag character for the common name of a series of linked plots
#' @param ... aesthetic mappings to be passed to ggplot2 methods
#' @importFrom grDevices dev.off pdf
#' @importFrom stats asOneSidedFormula
#' @importFrom utils packageVersion
#' @export
#' @examples
#' data(vsfuk2012)
#' vhist(FertilityRate, vsfuk2012, "hist1", "vsfuk2012", fill=Type)
#' vlaunch(vsfuk2012, "main", "vsfuk2012", browse=FALSE)
#'

vhist <- function(x, data, name, tag, ...){

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

  jspath <- file.path(system.file(package="vdmR"), "exec/vdmr_hist.js")
  file.copy(jspath, paste(name, ".", tag, ".js", sep=""), overwrite=TRUE)

  pdf(file=NULL, width=7, height=5)

  histqp <- ggplot2::ggplot(data, aesthetics)
  if(packageVersion("ggplot2")>'1.0.1'){
    params$na.rm <- TRUE
    histqp <- histqp + ggplot2::layer(geom="bar", stat="bin",
                                      position="identity", params=params)
  } else {
    histqp <- histqp + ggplot2::layer(geom="histogram", geom_params=params)
  }

  # retrieving data from histogram
  histdata <- histqp$data[,as.character(arguments$x)]

  gghist <- print(histqp)
  histparam <- gghist$data[[1]]
  plotranges <- gghist$panel$ranges[[1]]

  # "unique" is for multiple colored histogram
  xmax <- unique(histparam$xmax)
  xmin <- unique(histparam$xmin)
  hdmtx <- matrix(rep(histdata, length(xmax)), ncol=length(xmax))
  brkmax <- sweep(hdmtx, 2, signif(xmax), "<")
  brkmin <- sweep(hdmtx, 2, signif(xmin), ">=")

  whichcls <- apply((brkmax&brkmin),1,function(x){match(TRUE,x)})

  count <- as.vector(tapply(histparam$count, factor(histparam$xmax), sum))

  grid::grid.force()
  
  grid::downViewport("panel.3-4-3-4")

  dvp <- grid::dataViewport(xscale=plotranges$x.range, yscale=plotranges$y.range)

  grid::grid.rect(unique(histparam$x), 0, xmax-xmin,
                   0, default.units="native", gp=grid::gpar(fill="red", color=NA), name="hlbar", vp=dvp)

  grid::pushViewport(dvp)
  
  grid::grid.gedit("geom_rect.rect", name="geom_rect.rect")

  grid::upViewport(2)

  gridSVG::grid.script(paste("var xmin=",rjson::toJSON(xmin),";"))
  gridSVG::grid.script(paste("var xmax=",rjson::toJSON(xmax),";"))
  gridSVG::grid.script(paste("var count=",rjson::toJSON(count), ";"))
  gridSVG::grid.script(paste("var data=", rjson::toJSON(histdata), ";"))
  gridSVG::grid.script(paste("var whichcls=", rjson::toJSON(whichcls), ";"))
  gridSVG::grid.script(file=paste(name, ".", tag,".js", sep=""))
  gridSVG::grid.script(paste("var winname= '", name, "';", sep=""))

  histgrob <- grid::grid.grab(wrap=TRUE)

  grid::grid.newpage()
  grid::grid.draw(histgrob)

  svgfn <- paste0(name, ".", tag, ".svg")

  e <- try(gridSVG::grid.export(svgfn, htmlWrapper=TRUE, exportMappings="file",
                                xmldecl="<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"),
           silent=TRUE)

  if (class(e) == "try-error") {
    gridSVG::grid.export(svgfn, htmlWrapper=TRUE, exportMappings="file",
                         xmldecl="<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
  }

  invisible(dev.off())

}
