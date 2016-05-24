#' Show a three- or two-dimensional plot of a prcomp object
#' 
#' Show a three- two-dimensional plot of a prcomp object or a matrix, using
#' different symbols and colors for groups of data
#'
#' The pca3d function shows a three dimensional representation of a PCA object or any
#' other matrix. It uses the rgl package for rendering.
#'
#' pca2d is the 2D counterpart. It creates a regular, two-dimensional plot
#' on the standard graphic device. However, it takes exactly the same
#' options as pca3d, such that it is easy to create 2D variants of the 3D
#' graph.
#'
#' Often, PCA visualisation requires using different symbols and colors for
#' different groups of data. pca3d() and pca2d() aim at creating reasonable
#' defaults, such that a simple call with two parameters -- the pca object
#' and the vector with group assignments of the samples -- is sufficient for
#' a basic diagnosis.
#'
#' @section Biplots:
#'
#'   If option \option{biplot} is TRUE, a biplot showing both the PCA results
#'   (samples) and variables is shown. This corresponds to the
#'   \code{\link{biplot}} function which works for the \code{\link{prcomp}}
#'   class objects. However, a biplot showing all variable loadings will be
#'   unreadable if the data is highly dimensional (for example, gene
#'   expression data). Therefore, the option \option{biplot.vars} specifies
#'   which variables are shown on the biplot.
#'
#'   If \option{biplot.vars} is a vector of length larger than one, it will be
#'   interpreted as a direct selection of the variables to be shown; for
#'   example, for a \code{\link{prcomp}} object \var{pca}, the variable selection will
#'   happen through \code{pca$rotation[biplot.vars,]}.
#'
#'   If \option{biplot.vars} is a single number, then for each of the
#'   components shown, a number of variables equal to \option{biplot.vars}
#'   with the highest absolute loadings will be shown on the biplot.
#'
#' @param pca 
#'     Either a prcomp object or a matrix with at least three columns
#' @param components 
#'   Vector of length 3 (\code{pca3d}) or 2 (\code{pca2d}) containing the components to be shown
#' @param col 
#'   Either a single value or a vector of length equal to number of rows,
#'   containing color definitions for the plot points to be shown
#' @param title 
#'   Window title
#' @param new 
#'   Use TRUE to open a new window
#' @param axes.color 
#'   Axis color
#'   This option has no effect in pca2d.
#' @param bg 
#'   Background color
#' @param palette 
#'     Specifies the color palette when colors are automatically assigned to
#'     the groups. See Details.
#' @param fancy 
#'     set \option{show.labels}, \option{show.shadows},
#'     \option{show.centroids} and \option{show.group.labels} to TRUE.
#' @param radius 
#'   Scaling item for the size of points to be shown.
#'   In pca2d, this corresponds to the cex parameter.
#' @param biplot 
#'     Specify whether to show a biplot (see section \sQuote{biplots} below)
#' @param biplot.vars 
#'     Specify which loading to show on the biplot (see section
#'     \sQuote{biplots} below)
#' @param legend
#'     If NULL, no legend will be drawn. Otherwise the value specifies the
#'     legend position in a form accepted by \code{\link{legend}} and \code{\link{legend3d}}.
#' @param group 
#'   either NULL or a factor of length equal to number of rows. Factor levels
#'   can be used to automatically generate symbols and colors for the points
#'   shown
#' @param shape 
#'   Either a single value or a character vector describing the shapes to be
#'   used when drawing data points. Allowed shapes are: sphere, tetrahaedron
#'   and cube, and may be abbreviated.
#'   In pca2d, the parameter is passed directly on to the pch option of the
#'   points() function.
#' @param show.labels 
#'   TRUE for showing labels (taken from the coordinate matrix or the prcomp
#'   object). Alternatively, a vector with labels to be shown next to the data
#'   points.
#' @param labels.col 
#'   Single value or vector describing the colors of the labels.
#' @param show.scale 
#'   TRUE for showing a numeric scale at the edges of the plot.
#'   This option has no effect in pca2d.
#' @param show.axes 
#'   TRUE to show the axes.
#'   This option has no effect in pca2d.
#' @param show.axe.titles 
#'   If TRUE, show axe titles (PC 1, PC 2 etc.)
#'   This option has no effect in pca2d.
#' @param axe.titles 
#'   A vector with two (pca2d) or three (pca3d) values containing the axe
#'   titles (corresponds to xlab and ylab in regular plot). If missing, but
#'   show.axe.titles is TRUE, axe titles will be generated automatically.
#' @param show.plane 
#'   If TRUE, show a grey horizontal plane at y = 0.
#'   This option has no effect in pca2d.
#' @param show.shadows 
#'   If TRUE, show a "lollipop" representation of the points on the y = 0
#'   plane: a vertical line joining the data point with the plane and a
#'   shadow.
#'   In pca2d, for each sample at (x,y), a grey line is drawn from (x,y) to
#'   (x,0).
#' @param show.centroids 
#'   If TRUE and the group variable is defined, show cluster centroids (using
#'   apropriate group symbols) and lines from each data point to the
#'   corresponding centroid.
#' @param show.group.labels 
#'     Either TRUTH/FALSE or a vector equal to the number of unique values in
#'     the \code{group} parameter.  If set, labels for each of the defined
#'     group will be shown at the group's centroid. If the value of the
#'     parameter is TRUE, then the group names will be taken from the
#'     \code{group} parameter. Otherwise, the values from this parameter will
#'     be used.
#' @param show.shapes 
#'       A TRUTH/FALSE value indicating whether the different symbols (shapes)
#'       for the shown data points should be plotted (default TRUE).
#' @param show.ellipses
#'       A TRUTH/FALSE value indicating whether to show confidence interval ellipses or
#'       ellipsoids around each defined group
#' @param ellipse.ci
#'       The confidence level of a pairwise confidence region for the CI.  The
#'       default is 0.95, for a 95% region.  This is used to control
#'       the size of the ellipse being plotted.  
#' @param ... 
#'     For pca2d, any further argument will be passed on to the points()
#'     function.
#' @return Both pca2d and pca3d return invisibly a data frame which can be used to generate
#' a legend for the figure. The data frame has as many rows as there are
#' groups, and column with the group name, assigned color and assigned shape.
#' @examples
#'   data( metabo )
#'   pca <- prcomp( metabo[,-1], scale.= TRUE )
#'
#'   pca3d( pca, group= metabo[,1] )
#'   pca2d( pca, group= metabo[,1] )
#'
#'   ## a bit more fancy:
#'   ## black background, white axes,
#'   ## centroids
#'   pca3d( pca, group= metabo[,1], 
#'     fancy= TRUE, bg= "black", 
#'     axes.color= "white", new= TRUE )
#' @keywords PCA princomp biplot prcomp
#' @name pca3d-package
#' @import rgl
#' @importFrom ellipse ellipse
#' @import graphics
#' @importFrom utils head
#' @importFrom stats cov
NULL
