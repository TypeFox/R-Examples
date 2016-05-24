#' Loading Fonts For The 'showtext' Package
#' 
#' This function loads fonts that will be used by the \pkg{showtext} package.
#' 
#' @export
#' @author Yixuan Qiu <\url{http://statr.me/}>
#' @examples load_fonts()
load_fonts = function()
{
    fontfile = system.file("fonts", "wqy-microhei.ttc.zip", package = "showtextdb")
    outdir = tempdir()
    unzip(fontfile, exdir = outdir)
    outfile = file.path(outdir, "wqy-microhei.ttc")

    sysfonts::font.add("wqy-microhei", outfile)
    
    invisible(NULL)
}
