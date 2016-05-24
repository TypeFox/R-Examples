## all defined shapes
all.shapes<- c("sphere", "tetrahedron", "cube", "octahedron")

#' Show all permissible shapes for pca3d
#'
#' Show all permissible shapes for pca3d
#'
#' Show all permissible shapes for the functions pca3d and pca2d. The
#' shapes may be abbreviated using (matching is done with
#' \code{\link{pmatch}}.
#' @return A data frame with permissible 3d shapes for plotting and their
#' pch counterparts is
#' returned invisibly.
#' @export
listShapes <- function() {
  pch <- shape2pch(all.shapes)
  df <- data.frame(Shape=all.shapes, pch=pch, stringsAsFactors=FALSE)
  rownames(df) <- NULL
  print(df)
  return(invisible(df))
}

## a default set of colors. Nice.
#' Default palette
#'
#' Default set of colors for the pca3d package. This is a
#' colorblind-friendly palette, following the R cookbook.
#'
#' The default palette contains 21 colors. 
#' @param n Number of colors to return
#' @param transparent character string which will be pasted to each color
#' @param d3 if true, no transparency information will be added to the
#' colors
#' @return A character vector with the color palette
#' @export
defaultPalettePCA3D <- function(n= NULL, transparent= NULL, d3=FALSE) {

  # source: R cookbook, http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
  pal <- "E69F00 56B4E9 009E73 F0E442 0072B2 D55E00 CC79A7 999999 E69F00 56B4E9 009E73 F0E442 0072B2 D55E00 CC79A7"

  pal <- unlist(strsplit(pal, ' '))
  pal <- paste0("#", pal)
  if(!is.null(transparent)) {
    pal <- paste0(pal, transparent)
  }

  if(! is.null(n)) {
    if(n > length(pal)) {
      pal <- rep(pal, ceiling(n / length(pal)))
    } 
    pal <- pal[ 1:n ]
  }

  return(pal)
}
## little helpers
printf <- function(...) print(sprintf(...)) 
catf   <- function(...) cat(sprintf(...)) 

## determine average position in each group
calc.centroids <- function(coords, group) {

  centr.coords <- apply(coords, 2, function(x) tapply(x, group, mean))
  if(length(unique(group)) == 1)
    centr.coords <- matrix(centr.coords, nrow=1)

  rownames(centr.coords) <- levels(group)

  return(centr.coords)
}

## extract components, check dimensionality etc.
get.pca.coords <- function(pca, n, components=1:n) {

  if(length(components) < n) {
    stop(sprintf("Length of components is %d, but should be %d",
      length(components), n))
  }
  
  # if too many components, we take only first n
  if(length(components) > n) {
    warning(sprintf("Using first %d components", n))
    components <- components[1:n]
  }

  ret <- switch(class(pca),
    prcomp=pca$x,
    princomp=pca$scores,
    matrix=pca,
    stop("pca must be either a matrix or a PCA object (prcomp, princomp etc.)"))

  if(ncol(ret) < n || ncol(ret) < max(components)) {
    stop(sprintf("Not enough dimensions: %d, but at least %d are needed",
       ncol(ret), max(components, n)))
  }

  return(ret[,components])
}

## extract coordinates of variable loadings
get.biplot.coords <- function(pca, biplot) {

  if(class(biplot) == "matrix") return(biplot)
  if(class(pca) == "prcomp")    return(pca$rotation)
  if(class(pca) == "princomp")  return(pca$loadings)

  stop("For a biplot, another matrix or a prcomp object as pca is needed")
}

## select variables for plotting biplots
get.biplot.vars <- function(biplot.coords, biplot.vars) {

  if(length(biplot.vars) > 1) 
    return(biplot.vars)

  if(biplot.vars > nrow(biplot.coords)) biplot.vars <- nrow(biplot.coords)

  res <- c()
  for(i in 1:ncol(biplot.coords)) {
    res <- c(res, head(order(abs(biplot.coords[,i]), decreasing= TRUE), biplot.vars))
  }

  return(unique(res))
}

## print out the legend
print.legend <- function(group, group.col, group.shape, print=FALSE) {

  l.g <- max(nchar(c(levels(group), unique(group.col), unique(group.shape))))
  fmt <- sprintf("%% %ds", l.g)
  fmt <- sprintf("%s: %s, %s\n", fmt, fmt, fmt)

  if(print) {
    cat("\nLegend:\n")
    catf("%s\n", paste(rep("-", 3*l.g + 4), collapse= ""))
    catf(fmt, "group", "color", "shape")
    catf("%s\n", paste(rep("-", 3*l.g + 4), collapse= ""))
    for(i in 1:length(levels(group))) {
      catf(fmt, levels(group)[i], group.col[i], group.shape[i])
    }
    cat("\n")
  }

  ret <- data.frame(groups= levels(group), colors= group.col, 
    shapes= group.shape, 
    pch= shape2pch(group.shape),
    stringsAsFactors=FALSE)

  return(invisible(ret))
}




## convert 3d shapes to pch numeric vector for use with 2D plots
shape2pch <- function(sh) {

  if(is.numeric(sh)) return(sh)

  sh <- as.character(sh)
  
  # match shapes against known shapes
  ret <- all.shapes[pmatch(sh, all.shapes, duplicates.ok=TRUE)]

  # map to numeric values
  map <- c( sphere=16, tetrahedron=17, cube=15, octahedron=18)
  ret <- map[ret]
  if(any(is.na(ret))) {
    incorrect <- paste(unique(sh[is.na(ret)]), collapse=" ")
    warning(sprintf("Incorrect shapes: %s", incorrect))
  }
  ret
}


## find shapes by matching abbreviations
matchShapes <- function(sh) {

  if(is.numeric(sh)) {
    return(sh)
  }

  sh <- as.character(sh)
  ret <- all.shapes[pmatch(sh, all.shapes, duplicates.ok=TRUE)]

  if(any(is.na(ret))) {
    incorrect <- paste(unique(sh[is.na(ret)]), collapse=" ")
    stop(sprintf("Incorrect shapes: %s", incorrect))
  }
  ret
}


## automatically assing values (e.g. color) depending on the provided
## argument and a palette of character values to choose from
autocol <- function(n, col=NULL, group, palette, msg="col") {

  g.n <- as.numeric(group)

  if(is.null(col)) 
    return(palette[ (g.n-1) %% length(palette) + 1 ])

  if(length(col) == 1)
    return(rep(col, n))

  if(length(col) == n)
    return(col)

  if(length(col) == length(levels(group)) &&
     setequal(names(col), levels(group))) {

    return(col[as.character(group)])
  }

  stop(sprintf("Invalid %s vector length: %d (must be 1 or %d or a named vector of %d)",
    msg, length(col), n, length(levels(group))))
}

## based on a vector of values, return the first value from a group for a
## mapping
getGroupVals <- function(var, group) {

  group.val        <- var[ match(levels(group), group) ]
  names(group.val) <- levels(group)

  group.val
}

