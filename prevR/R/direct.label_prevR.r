#' Direct label on a ggplot object
#' 
#' Direct label a ggplot2 grouped plot
#' 
#' @param p The ggplot object.
#' @param method Method for direct labeling as described in ?label.positions.
#' @param debug Show debug output?
#' 
#' @return The ggplot object with direct labels added.
#' 
#' @note This function is based on and similar to \code{\link[directlabels]{direct.label.ggplot}}
#' except that legend is not hidden.
#' 
#' @seealso \code{\link[directlabels]{direct.label.ggplot}}\{\pkg{directlabels}\}
#' 
#' @export

direct.label_prevR <- function (p, method = NULL, debug = FALSE) 
# based on http://www.rdocumentation.org/packages/directlabels/functions/direct.label.ggplot
{
  getData <- function(colour.or.fill) {
    for (L in p$layers) {
      m <- p$mapping
      m[names(L$mapping)] <- L$mapping
      colvar <- m[[colour.or.fill]]
      if (!is.null(colvar)) {
        return(list(layer = L, colvar = as.character(colvar)))
      }
    }
  }
  dl.info <- getData("colour")
  if (is.null(dl.info)) {
    dl.info <- getData("fill")
  }
  if (is.null(dl.info)) {
    stop("Need colour or fill aesthetic to infer default direct labels.")
  }
  L <- dl.info$layer
  colvar <- dl.info$colvar
  geom <- L$geom$objname
  if (is.null(method)) 
    method <- default.picker("ggplot")
  data <- if ((!is.null(L$data)) && (length(L$data) > 0)) {
    L$data
  }
  else {
    NULL
  }
  a <- aes_string(label = colvar, colour = colvar)
  a2 <- structure(c(L$mapping, a), class = "uneval")
  dlgeom <- geom_dl(a2, method, stat = L$stat, debug = debug, 
                    data = data)
  dlgeom$stat_params <- L$stat_params
  # leg.info <- legends2hide(p)
  leg.info <- NULL # We want to keep the legend
  guide.args <- as.list(rep("none", length(leg.info$hide)))
  names(guide.args) <- leg.info$hide
  guide.args$colour <- "none"
  guide <- do.call(guides, guide.args)
  p + dlgeom + guide
}
