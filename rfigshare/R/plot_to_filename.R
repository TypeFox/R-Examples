#' Convienence function to save a ggplot2 plot, and return its filename.
#' 
#' @param plotobj ggplot2 plot object (should add support for base plots too)
#' @param filename Filename, don't include the file type extension.
#' @param path Path where you want to save the file. 
#' @return A file name, to use in fs_upload
#' @importFrom ggplot2 ggsave
#' @seealso \code{\link{fs_upload}}
#' @examples \dontrun{ 
#' # include in your fs_upload call
#' library(ggplot2)
#' p <- qplot(mpg, wt, data=mtcars)
#' plott <- fs_create(title="my title", description="some description", type="figure")
#' fs_add_categories(plott, "Ecology")
#' fs_upload(plott, plot_to_filename(p, "myfilename", "~"))
#' }
#' @export
plot_to_filename <- function(plotobj, filename, path = ".")
{
  filename2 <- paste0(filename, ".png")
  if(any(class(plotobj)=="ggplot"))
    ggsave(paste0(path, "/", filename2), plotobj)
  else
    message("figure out how to save base plots inside a function")
  #     png(paste0(path, "/", filename2))
  #     plotobj
  #     dev.off()
  return(filename2)
}