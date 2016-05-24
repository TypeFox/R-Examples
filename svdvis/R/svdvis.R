#' Scree plot
#'
#' Creates a scree plot, where percentages of variance explained (PVE) by singular vectors are plotted.
#' Additional customizations can be done by adding \code{ggplot2} functions to the output.
#'
#' @param svd.obj A list, resulted from applying svd to a dataset, with \code{u}, \code{d}, and \code{v} corresponding to left singular vector, singular values, and right singular vectors, respectively. Alternatively, supply a vector of singular values, \code{d}.
#' @param subr An optional positive integer to display PVE corresponding to the first \code{subr} singular values.
#' @param maintitle A figure title (by default, "Scree Plot").
#' @param axis.title.x A title for x axis (by default, "Singular Vectors").
#' @param axis.title.y A title for y axis (by default, "Percent Variance Explained").
#'
#' @return \code{svd.scree} creates and draws a figure, which is a \code{ggplot2} when \code{subr=NULL} or a {gtable} object when \code{subr} is specified.
#'
#' @author Neo Christopher Chung \email{nchchung@gmail.com}
#'
#' @export svd.scree
#' @importFrom gridExtra grid.arrange
#' @import ggplot2
#'
#' @examples
#' set.seed(1234)
#' dat = matrix(rnorm(1000), 100, 10)
#' svd.obj = svd(dat)
#' colnames(svd.obj$v) = paste0("V",1:10)
#' svd.scree(svd.obj)
svd.scree <- function(svd.obj, subr=NULL, maintitle="Scree Plot", axis.title.x="Singular Vectors", axis.title.y="Percent Variance Explained") {
  if(is.list(svd.obj) & all(names(svd.obj) %in% c("u","d","v"))) {
    print("Your input data is treated as a SVD output, with u, d, v corresponding to left singular vector, singular values, and right singular vectors, respectively.")
  } else {
    print("Your input data is treated as a vector of singular values. For example, it should be svd.obj$d from a SVD output.")
    svd.obj = list(d=svd.obj)
  }

  print("Scree Plot")

  if(is.null(names(svd.obj$d))) {
    names(svd.obj$d) = paste0("V",1:ncol(svd.obj$v))
  }

  pve = svd.obj$d^2/sum(svd.obj$d^2) * 100
  pve = data.frame(names=names(svd.obj$d), pve)
  pve$names = factor(pve$names, levels=unique(names(svd.obj$d)))
  g = ggplot(pve, aes(names, pve)) + geom_point() + theme_bw()
  gout = g + ylim(0,NA) +
    labs(title=maintitle, axis.title.x=axis.title.x, axis.title.y=axis.title.y)

  if(!is.null(subr)) {
    gsub = g + coord_cartesian(xlim = c(0.5, subr+.5)) + ylim(min(pve$pve[1:subr])-1, max(pve$pve[1:subr])+1) +
      labs(title=paste("First",subr,"singular values"), axis.title.x=axis.title.x, axis.title.y=axis.title.y)
    gout = grid.arrange(gout, gsub, nrow=1)
  }

  return(gout)
}

#' Visualizing Singular Vectors or Principal Components by Scatterplot Matrices
#'
#' Creates a set of multiple scatter plots from all pairs of selected singular vectors or principal components.
#' Principal components can be plotted by setting \code{weights = "sv"}.
#' Since it largely uses \code{ggpairs} from the \code{GGally} package, optional arguments for \code{ggpairs} can be specified.
#'
#' @param svd.obj A list, resulted from applying svd to a dataset, with \code{u}, \code{d}, and \code{v} corresponding to left singular vector, singular values, and right singular vectors, respectively. Alternatively, supply singular vectors, \code{v}.
#' @param r A positive integer to use only the first \code{r} vectors in visualization. If not specified, all vectors available in \code{svd.obj$v} are visualized.
#' @param group A vector of length \code{n}, specifying groups (e.g., phenotypes or conditions for \code{n} samples).
#' @param weights A vector of length \code{r}. If "sv", singular values contained in \code{svd.obj$d[1:r]} are used.
#' @param alpha A numeric value for transparency.
#' @param axisLabels Set to either "none" (default), "show", or "internal".
#' @param ... Additional arguments to pass onto \code{ggpair}.
#'
#' @return \code{svd.scatter} creates and draws a figure, which is a \code{ggpair} object.
#'
#' @author Neo Christopher Chung \email{nchchung@gmail.com}
#'
#' @export svd.scatter
#' @importFrom GGally ggpairs
#' @import ggplot2
#'
#' @seealso \link{ggpairs}
#' @examples
#' set.seed(1234)
#' dat = matrix(rnorm(1000), 100, 10)
#' svd.obj = svd(dat)
#' colnames(svd.obj$v) = paste0("V",1:10)
#' svd.scatter(svd.obj, r=3, group=c(rep("Group1",5), rep("Group2",5)))
svd.scatter <- function(svd.obj, r=NULL, group=NULL, weights=NULL, alpha=0.7, axisLabels="none", ...) {
  if(is.list(svd.obj) & all(names(svd.obj) %in% c("u","d","v"))) {
    print("Your input data is treated as a SVD output, with u, d, v corresponding to left singular vector, singular values, and right singular vectors, respectively.")
  } else {
    print("Your input data is treated as (typically, right) singular vectors. For example, it should be svd.obj$v  from a SVD output.")
    svd.obj = list(v=svd.obj)
  }

  print("Multiple Scatter Plots")

  if(is.null(dimnames(svd.obj$v))) {
    colnames(svd.obj$v) = paste0("V",1:ncol(svd.obj$v))
    rownames(svd.obj$v) = paste0("ID ",1:nrow(svd.obj$v))
  }

  if(is.null(r)) r=ncol(svd.obj$v)
  if(r > 15) print("It may not be good to visualize too many singular vectors or principal components at one.")
  if(is.null(group)) group=rep(1,nrow(svd.obj$v))
  if(!is.null(weights)) {
    if(weights == "sv") {
      print("Singular values are used as weights.")
      weights = svd.obj$d[1:r]
    } else if (length(weights) != r) {
      stop("The length of weights must equal r.")
    }
  }

  if(!is.null(weights)) {
    mv = as.data.frame(weights * svd.obj$v[, 1:r])
    mv = cbind(mv, group)
  } else {
    mv = as.data.frame(svd.obj$v[, 1:r])
    mv = cbind(mv, group)
  }

  g = ggpairs(mv, columns = 1:(r), color="group", alpha=alpha, axisLabels=axisLabels, diag=list(continuous="density"), upper=list(continuous='blank'), ...)
  return(g)
}

#' Visualizing Singular Vectors or Principal Components by Heatmaps
#'
#' Creates a heatmap from selected singular vectors or principal components.
#' Principal components can be plotted by setting \code{weights = "sv"}.
#' Colors for heatmap can be specified by optional arguments \code{low} and \code{high} colors.
#'
#' @param svd.obj A list, resulted from applying svd to a dataset, with \code{u}, \code{d}, and \code{v} corresponding to left singular vector, singular values, and right singular vectors, respectively. Alternatively, supply singular vectors, \code{v}.
#' @param r A positive integer to use only the first \code{r} vectors in visualization. If not specified, all vectors available in \code{svd.obj$v} are visualized.
#' @param group A vector of length \code{n}, specifying groups (e.g., phenotypes or conditions for \code{n} samples).
#' @param weights A vector of length \code{r}. If "sv", singular values contained in \code{svd.obj$d[1:r]} are used.
#' @param alpha A numeric value for transparency.
#' @param low A hex color code to color the lowest value.
#' @param high A hex color code to color the highest value.
#'
#' @return \code{svd.heatmap} creates and draws a figure, which is a \code{ggplot} object.
#'
#' @author Neo Christopher Chung \email{nchchung@gmail.com}
#'
#' @export svd.heatmap
#' @importFrom reshape2 melt
#' @import ggplot2
#'
#' @examples
#' set.seed(1234)
#' dat = matrix(rnorm(1000), 100, 10)
#' svd.obj = svd(dat)
#' colnames(svd.obj$v) = paste0("V",1:10)
#' rownames(svd.obj$v) = paste0("Sample",1:10)
#' svd.heatmap(svd.obj, r=5)
svd.heatmap <- function(svd.obj, r=NULL, group=NULL, weights=NULL, alpha=0.7, low = "#FFFFFF", high = "#9E0142") {
  if(is.list(svd.obj) & all(names(svd.obj) %in% c("u","d","v"))) {
    print("Your input data is treated as a SVD output, with u, d, v corresponding to left singular vector, singular values, and right singular vectors, respectively.")
  } else {
    print("Your input data is treated as (typically, right) singular vectors. For example, it should be svd.obj$v  from a SVD output.")
    svd.obj = list(v=svd.obj)
  }

  print("SVD Heatmap")

  if(is.null(dimnames(svd.obj$v))) {
    colnames(svd.obj$v) = paste0("V",1:ncol(svd.obj$v))
    rownames(svd.obj$v) = paste0("ID ",1:ncol(svd.obj$v))
  }

  if(is.null(r)) r=ncol(svd.obj$v)
  if(r > 30) print("It may not be good to visualize too many singular vectors or principal components at one.")
  if(is.null(group)) group=rep(1,nrow(svd.obj$v))
  if(!is.null(weights)) {
    if(weights == "sv") {
      print("Singular values are used as weights.")
      weights = svd.obj$d[1:r]
    } else if (length(weights) != r) {
      stop("The length of weights must equal r.")
    }
  }

  if(!is.null(weights)) {
    mv = as.data.frame(weights * svd.obj$v[, 1:r])
    #mv = cbind(mv, group)
  } else {
    mv = as.data.frame(svd.obj$v[, 1:r])
    #mv = cbind(mv, group)
  }

  mv$rnames <- rownames(mv)
  mv = melt(as.data.frame(mv), id.vars="rnames", value.name="value", variable.name="cnames")
  g = ggplot(mv, aes_string(x="rnames", y="cnames")) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    theme(legend.position = "bottom", axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0, colour = "grey50")) +
    geom_tile(aes_string(fill = "value"), colour = "white") +
    scale_fill_gradient("", low = low, high = high)

  return(g)
}

#' Visualizing Singular Vectors or Principal Components by Parallel Coordinates Plots
#'
#' Creates a Parallel Coordinates Plot from selected singular vectors or principal components.
#' Principal components can be plotted by setting \code{weights = "sv"}.
#' Since it largely uses \code{ggparcoord} from the \code{GGally} package, optional arguments for \code{ggparcoord} can be specified.
#'
#' @param svd.obj A list, resulted from applying svd to a dataset, with \code{u}, \code{d}, and \code{v} corresponding to left singular vector, singular values, and right singular vectors, respectively. Alternatively, supply singular vectors, \code{v}.
#' @param r A positive integer to use only the first \code{r} vectors in visualization. If not specified, all vectors available in \code{svd.obj$v} are visualized.
#' @param group A vector of length \code{n}, specifying groups (e.g., phenotypes or conditions for \code{n} samples).
#' @param weights A vector of length \code{r}. If "sv", singular values contained in \code{svd.obj$d[1:r]} are used.
#' @param alpha A numeric value for transparency.
#' @param ... Additional arguments to pass onto \code{ggparcoord}.
#'
#' @return \code{svd.parallel} creates and draws a figure, which is a \code{ggplot} object.
#'
#' @author Neo Christopher Chung \email{nchchung@gmail.com}
#'
#' @export svd.parallel
#' @importFrom GGally ggparcoord
#' @import ggplot2
#'
#' @seealso \link{ggparcoord}
#' @examples
#' set.seed(1234)
#' dat = matrix(rnorm(1000), 100, 10)
#' svd.obj = svd(dat)
#' colnames(svd.obj$v) = paste0("V",1:10)
#' rownames(svd.obj$v) = paste0("Sample",1:10)
#' svd.parallel(svd.obj, r=4)
svd.parallel <- function(svd.obj, r=NULL, weights=NULL, group=NULL, alpha=.7, ...) {
  if(is.list(svd.obj) & all(names(svd.obj) %in% c("u","d","v"))) {
    print("Your input data is treated as a SVD output, with u, d, v corresponding to left singular vector, singular values, and right singular vectors, respectively.")
  } else {
    print("Your input data is treated as (typically, right) singular vectors. For example, it should be svd.obj$v  from a SVD output.")
    svd.obj = list(v=svd.obj)
  }

  print("Parallel Coordinates Plot")

  if(is.null(dimnames(svd.obj$v))) {
    colnames(svd.obj$v) = paste0("V",1:ncol(svd.obj$v))
    rownames(svd.obj$v) = paste0("ID ",1:ncol(svd.obj$v))
  }

  if(is.null(r)) r=ncol(svd.obj$v)
  if(r > 15) print("It may not be good to visualize too many singular vectors or principal components at one.")
  if(is.null(group)) group=rep(1,nrow(svd.obj$v))
  if(!is.null(weights)) {
    if(weights == "sv") {
      print("Singular values are used as weights.")
      weights = svd.obj$d[1:r]
    } else if (length(weights) != r) {
      stop("The length of weights must equal r.")
    }
  }

  if(!is.null(weights)) {
    mv = as.data.frame(weights * svd.obj$v[,1:r])
    mv = cbind(mv, group)
  } else {
    mv = as.data.frame(svd.obj$v[,1:r])
    mv = cbind(mv, group)
  }

  g = ggparcoord(mv, columns=c(1:r), groupColumn=r+1, alphaLines=alpha,  ...)
  return(g)
}

#' Visualizing Singular Vectors or Principal Components by Radial Coordinates Plots
#'
#' Creates a Radial Coordinates Plot from selected singular vectors or principal components.
#' Principal components can be plotted by setting \code{weights = "sv"}.
#' It uses the \code{radviz} function, optional arguments for \code{radviz} can be specified.
#'
#' @param svd.obj A list, resulted from applying svd to a dataset, with \code{u}, \code{d}, and \code{v} corresponding to left singular vector, singular values, and right singular vectors, respectively. Alternatively, supply singular vectors, \code{v}.
#' @param r A positive integer to use only the first \code{r} vectors in visualization. If not specified, all vectors available in \code{svd.obj$v} are visualized.
#' @param group A vector of length \code{n}, specifying groups (e.g., phenotypes or conditions for \code{n} samples).
#' @param weights A vector of length \code{r}. If "sv", singular values contained in \code{svd.obj$d[1:r]} are used.
#' @param alpha A numeric value for transparency.
#' @param ... Additional arguments to pass onto \code{radviz}.
#'
#' @return \code{svd.radial} creates and draws a figure, which is a \code{ggplot} object.
#'
#' @author Neo Christopher Chung \email{nchchung@gmail.com}
#'
#' @export svd.radial
#' @importFrom GGally ggparcoord
#' @import ggplot2
#'
#' @seealso \link{radviz}
#'
#' @examples
#' set.seed(1234)
#' dat = matrix(rnorm(1000), 100, 10)
#' svd.obj = svd(dat)
#' colnames(svd.obj$v) = paste0("V",1:10)
#' rownames(svd.obj$v) = paste0("Sample",1:10)
#' svd.radial(dat, group=c(rep("Group1",5), rep("Group2",5)))
svd.radial <- function(svd.obj, r=NULL, weights=NULL, group=NULL, alpha=1, ...) {
  print("Radial Visualization Plots")
  if(is.list(svd.obj) & all(names(svd.obj) %in% c("u","d","v"))) {
    print("Your input data is treated as a SVD output, with u, d, v corresponding to left singular vector, singular values, and right singular vectors, respectively.")
  } else {
    print("Your input data is treated as (typically, right) singular vectors. For example, it should be svd.obj$v  from a SVD output.")
    svd.obj = list(v=svd.obj)
  }

  if(is.null(dimnames(svd.obj$v))) {
    colnames(svd.obj$v) = paste0("V",1:ncol(svd.obj$v))
    rownames(svd.obj$v) = paste0("ID ",1:nrow(svd.obj$v))
  }

  if(is.null(r)) r=ncol(svd.obj$v)
  if(r > 10) print("It may not be good to visualize too many singular vectors or principal components at one.")
  if(is.null(group)) group=rep(1,nrow(svd.obj$v))
  if(!is.null(weights)) {
    if(weights == "sv") {
      print("Singular values are used as weights.")
      weights = svd.obj$d[1:r]
    } else if (length(weights) != r) {
      stop("The length of weights must equal r.")
    }
  }

  if(!is.null(weights)) {
    mv = as.data.frame(weights * svd.obj$v[,1:r])
  } else {
    mv = as.data.frame(svd.obj$v[,1:r])
  }
  return(radviz(mv, group=group, alpha=alpha,  ...))
}
