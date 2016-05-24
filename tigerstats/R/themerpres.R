#' Lattice Theme or R Presentations
#' 
#' Modifies the current theme for use with lattice graphics in R Presentation dicuments. Increases size of title,
#'  axis lables and axis numbers, thickens sone lines, etc.
#'
#'@usage themerpres()
#' 
#' @return Returns a list to be supplied as the \code{theme} to the \code{lattice} function
#' \code{\link{trellis.par.set}()}.
#'
#'
#' @seealso \code{\link{trellis.par.set}}, \code{\link{show.settings}} 
#' 
#' @rdname themerpres
#'
#' @export
#' @examples
#' trellis.par.set(theme=themerpres())
#' 
#' @keywords graphics 
#' 
themerpres <- function () {
  my.theme <- lattice::trellis.par.get() #Would like to be able to pick up the default theme
  my.theme$add.text[[2]] <- 2.5
  my.theme$box.rectangle[[5]] <- 3
  my.theme$box.umbrella[[4]] <- 3
  my.theme$box.dot[[3]] <- 2
  my.theme$box.dot[[5]] <- 19
  my.theme$plot.line[[4]] <- 3
  my.theme$axis.text[[2]] <- 2.5
  my.theme$par.xlab.text[[2]] <- 2.5
  my.theme$par.ylab.text[[2]] <- 2.5
  my.theme$par.main.text[[2]] <- 3
  return(my.theme)
}

