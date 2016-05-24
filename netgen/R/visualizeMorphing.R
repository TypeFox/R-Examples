#' Fancy visualization of morphing.
#'
#' Takes two instances of equal size and some alpha values. Computes the point matching
#' and morphings for the alpha values and visualizes the transition of points
#' of the first instance towards their matched counterparts of the second instance
#' with two different methods.
#'
#' @template arg_first_network
#' @template arg_second_network
#' @param point.matching [\code{matrix}]\cr
#'   Point matching which shall be used for morphing. If \code{NULL}, an optimal
#'   point matching is generated via function \code{\link{getOptimalPointMatching}}.
#'   Default is \code{NULL}.
#' @param alphas [\code{numeric}]\cr
#'   Vector of coefficients 'alpha' for convex combinations.
#' @param arrows [\code{logical(1)}]\cr
#'   Draw arrows originating in the points of \code{x} and ending in the
#'   points matched in \code{y}. Default is \code{TRUE}.
#' @param in.one.plot [\code{logical(1)}]\cr
#'   Currently the function offers two different types of plot. If \code{in.one.plot}
#'   is \code{TRUE}, which is the default value, the morphing is dipicted in one plot.
#'   This is in particular useful for small instances. If set to \code{FALSE},
#'   a matrix of plots is generated via \code{\link[ggplot2]{facet_grid}}. One
#'   plot for each alpha value in \code{alphas}.
#' @return [\code{\link[ggplot2]{ggplot}}]
#' @seealso \code{\link{morphInstances}}
#' @export
visualizeMorphing = function(x, y,
  point.matching = NULL,
  alphas = c(0.25, 0.5, 0.75), arrows = TRUE, in.one.plot = TRUE) {
  assertClass(x, "Network")
  assertClass(y, "Network")
  assertNumeric(alphas, any.missing = FALSE, lower = 0, upper = 1)
  assertFlag(arrows)
  assertFlag(in.one.plot)

    # we compute the point matching here one time additionaly for later use
  if (is.null(point.matching)) {
    point.matching = getOptimalPointMatching(x$coordinates, y$coordinates)
  }
  assertMatrix(point.matching, mode = "numeric", any.missing = FALSE)

  # do the morphing
  instances = lapply(alphas, function(alpha) {
    morphInstances(x, y, alpha = alpha)
  })

    # put in the original instances as extreme instances
  alphas = c(0, alphas, 1)
  instances = c(list(y), instances, list(x))

  # build data frame of points
  df.points = lapply(1:length(instances), function(i) {
    df = as.data.frame(instances[[i]], include.extra = TRUE)
    # we need to add membership here for ggplot
    # This is only the case if we morph a clustered and a non-clustered instance
    if (is.null(df$membership)) {
      df$membership = 1L
    }
    df$alpha = alphas[i]
    df
  })

  df.points = do.call(rbind, df.points)
  df.points$membership = as.factor(df.points$membership)
  df.points$alpha = as.factor(df.points$alpha)

  lower = x$lower
  upper = x$upper

  #FIXME: outsource both variants to seperate files
  if (in.one.plot) {
    pl = ggplot(data = df.points, mapping = aes_string(x = "x1", y = "x2", colour = "alpha"))

    if (arrows) {
      # build arrows from x to (matched) y coordinates
      df.arrows = cbind(as.data.frame(x$coordinates), as.data.frame(y$coordinates[point.matching[, 2], ]))
      colnames(df.arrows) = c("x1", "x2", "end1", "end2")
      pl = pl + geom_segment(
        data = df.arrows,
        mapping = aes_string(x = "x1", y = "x2", xend = "end1", yend = "end2"),
        arrow = grid::arrow(type = "closed", length = grid::unit(0.1, "inches")),
        colour = "gray"
      )
    }
    pl = pl + geom_point()
    pl = pl + geom_point(data = df.points[(df.points$types == "depot"), ], colour = "black", size = 4)
    pl = pl + geom_point(data = df.points[(df.points$types == "depot"), ], colour = "white", size = 3)
  } else {
    # we want nice 'alpha = value' labels for the facets
    df.points$alpha = paste("alpha == ", df.points$alpha)
    pl = ggplot(data = df.points[(df.points$types != "depot"), ], mapping = aes_string(x = "x1", y = "x2"))
    pl = pl + geom_point(mapping = aes_string(colour = "membership"))
    pl = pl + geom_point(data = df.points[(df.points$types == "depot"), ], colour = "black", size = 4)
    pl = pl + geom_point(data = df.points[(df.points$types == "depot"), ], colour = "white", size = 3)
    # keep in mind labeller to parse expressions!
    pl = pl + facet_grid(. ~ alpha, labeller = label_parsed)
  }
  pl = pl + ggtitle(paste("#Nodes:", getNumberOfNodes(x)))
  pl = decorateGGPlot(pl, lower, upper)
  return(pl)
}
