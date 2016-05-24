#' Spectral Map Plot of Multivariate Data
#' Produces a spectral map plot (biplot) of an object of class \code{mpm}
#' 
#' Spectral maps are special types of biplots with the area of the symbols
#' proportional to some measure, usually the row or column mean value and an
#' identification of row- and column-items.  For large matrices, such as gene
#' expression data, where there is an abundance of rows, this can obscure the
#' plot.  In this case, the argument \code{label.tol} can be used to select the
#' most informative rows, i.e. rows that are most distant from the center of
#' the plot.  Only these row-items are then labeled and represented as circles
#' with their areas proportional to the marginal mean value.  For the
#' column-items it can be useful to apply some grouping specified by
#' \code{col.group}. Examples of groupings are different pathologies, such as
#' specified in \code{Golub.grp}
#' 
#' @param x object of class \code{mpm} a result of a call to \code{mpm}.
#' @param scale optional character string specifying the type of factor scaling
#'   of the biplot. This can be either \kbd{"singul"} (singular value scaling),
#'   \kbd{"eigen"} (eigenvalue scaling), \kbd{"uvr"} (unit row-variance
#'   scaling), \kbd{"uvc"} (unit column-variance scaling). The latter is of
#'   particular value when analyzing large matrices, such as gene expression
#'   data. Singular value scaling \kbd{"singul"} is customary in spectral map
#'   analysis. Defaults to \kbd{"singul"}.
#' @param dim optional principal factors that are plotted along the horizontal
#'   and vertical axis. Defaults to \code{c(1,2)}.
#' @param zoom optional zoom factor for row and column items. Defaults to
#'   \code{c(1,1)}.
#' @param show.row optional character string indicating whether all rows
#'   (\kbd{"all"}) are to be plotted or just the positioned rows
#'   \kbd{"position"}.
#' @param show.col optional character string indicating whether all columns
#'   (\kbd{"all"}) are to be plotted or just the positioned columns
#'   \kbd{"position"}.
#' @param col.group optional vector (character or numeric) indicating the
#'   different groupings of the columns, e.g. \code{Golub.grp}. Defaults to 1.
#' @param colors vector specifying the colors for the annotation of the plot;
#'   the first two elements concern the rows; the third till the last element
#'   concern the columns; the first element will be used to color the unlabeled
#'   rows; the second element for the labeled rows and the remaining elements
#'   to give different colors to different groups of columns.
#' @param col.areas logical value indicating whether columns should be plotted
#'   as squares with areas proportional to their marginal mean and colors
#'   representing the different groups (\code{TRUE}), or with symbols
#'   representing the groupings and identical size (\code{FALSE}). Defaults to
#'   \code{TRUE}.
#' @param col.symbols vector of symbols when \code{col.areas=FALSE} corresponds
#'   to the \code{pch} argument of the function \code{plot}.
#' @param sampleNames Either a logical vector of length one or a character
#'   vector of length equal to the number of samples in the dataset. If a
#'   logical is provided, sample names will be displayed on the plot
#'   (\code{TRUE}; default) or not (\code{FALSE}); if a character vector is
#'   provided, the names provided will be used to label the samples instead of
#'   the default column names.
#' @param rot rotation of plot. Defaults to \code{c(-1,-1)}.
#' @param labels character vector to be used for labeling points on the graph;
#'   if \code{NULL}, the row names of \code{x} are used instead
#' @param label.tol numerical value specifying either the percentile
#'   (\code{label.tol<=1}) of rows or the number of rows (\code{label.tol>1})
#'   most distant from the plot-center (0,0) that are labeled and are plotted
#'   as circles with area proportional to the marginal means of the original
#'   data.
#' @param label.col.tol numerical value specifying either the percentile
#'   (\code{label.col.tol<=1}) of columns or the number of columns (\code{label.col.tol>1})
#'   most distant from the plot-center (0,0) that are labeled and are plotted
#'   as circles with area proportional to the marginal means of the original
#'   data.
#' @param lab.size size of identifying labels for row- and column-items as
#'   \code{cex} parameter of the \code{text} function
#' @param col.size size in mm of the column symbols
#' @param row.size size in mm of the row symbols
#' @param do.smoothScatter use smoothScatter or not instead of plotting
#'   individual points
#' @param do.plot produce a plot or not
#' @param ... further arguments to \code{eqscaleplot} which draws the canvas
#'   for the plot; useful for adding a \code{main} or a custom \code{sub}
#' @return An object of class \code{plot.mpm} that has the following
#'   components: \item{Rows}{a data frame with the X and Y coordinates of the
#'   rows and an indication \code{Select} of whether the row was selected
#'   according to \code{label.tol}} \item{Columns}{a data frame with the X and
#'   Y coordinates of the columns}
#' @note \code{value} is returned invisibly, but is available for further use
#'   when an explicit assignment is made
#' @author Luc Wouters
#' @seealso \code{\link{mpm}}, \code{\link{summary.mpm}}
#' @references Wouters, L., Goehlmann, H., Bijnens, L., Kass, S.U.,
#'   Molenberghs, G., Lewi, P.J. (2003). Graphical exploration of gene
#'   expression data: a comparative study of three multivariate methods.
#'   \emph{Biometrics} \bold{59}, 1131-1140.
#' @keywords multivariate hplot
#' @method plot mpm
#' @S3method plot mpm
#' @examples
#' 
#'   # Weighted spectral map analysis
#'   data(Golub) # Gene expression data of leukemia patients
#'   data(Golub.grp) # Pathological classes coded as 1, 2, 3
#'   r.sma <- mpm(Golub[,1:39], row.weight = "mean", col.weight = "mean")
#'   # Spectral map biplot with result
#'   r <- plot(r.sma, label.tol = 20, scale = "uvc",
#'             col.group = (Golub.grp)[1:38], zoom = c(1,1.2), col.size = 5)
#'   Golub[r$Rows$Select, 1] # 20 most extreme genes
#' 
#' @export
plot.mpm <- function(
    x, # mpm object
    scale = c("singul", "eigen", "uvr", "uvc"), # Scaling type
    dim = c(1, 2), # Principal factors to plot.
    zoom = rep(1, 2), # Zoom factor for factor scores and loadings
    show.row = c("all", "position"),
    show.col = c("all", "position"),
    col.group = rep(1, length(x$col.names)),
    colors = c("orange1", "red", rainbow(length(unique(col.group)), start=2/6, end=4/6)),
    col.areas = TRUE,
    col.symbols = c(1, rep(2, length(unique(col.group)))),
    sampleNames = TRUE,
    rot = rep(-1, length(dim)), # Mirror all axes
    labels = NULL, # character vector of labels (to allow labels to differ from row.names)
    label.tol = 1,
    label.col.tol = 1, 
    lab.size = 0.725,
    col.size = 10,
    row.size = 10,
    do.smoothScatter = FALSE, # Plot individual points or density maps
    do.plot = TRUE, # This routine can also be used to calculate the 
                    # coordinates, without plotting
    # legend = FALSE,
    ...) 
{

  ### error checking and argument matching  
  if (missing(x))
    stop ("Argument 'x' is missing, with no default.")
  if (!inherits(x, "mpm")) 
    stop("Use only with 'mpm' objects.")
  if (!require("MASS"))
    stop("MASS package could not be loaded.") 
  if (do.smoothScatter && label.tol == 1){
    # if we require all points in the plot to be labelled, no points are left for the density map
    warning("All points selected for labelling, continuing without density map.\n")
    do.smoothScatter <- FALSE
  }
  if(do.smoothScatter && !do.plot){
    warning("Density map plotting requested but plotting not selected. Continuing with plot.\n")
    do.plot <- TRUE # If smoothScatter is set, then also do a plot
  }
  
  # labeling of samples on the plot
  if (is.logical(sampleNames)){
    if (length(sampleNames) > 1){
      stop("'sampleNames' should either be a logical of length one or a character of length equal to the number of observations")
    } else {
      sampleLabels <- if (sampleNames) x$col.names else rep("", length(x$col.names))  
    }
  } else {
    if (length(sampleNames) != length(x$col.names)){
      stop("If 'sampleNames' is a character vector, the length should be equal to the number of observations")
    } else {
      sampleLabels <- sampleNames
    }
  }
  
  
  scale <- match.arg(scale)
  show.row <- match.arg(show.row)
  show.col <- match.arg(show.col)
  if (is.data.frame(col.group))
    col.group <- unlist(col.group)
  if (length(col.group) != length(x$col.names))
    stop("Length of 'col.group' not equal to number of columns in data.")
  col.group <- as.numeric(as.factor(col.group))

  ### extract data from x  
  Wn <- x$Wn        # row weights
  Wp <- x$Wp        # column weights
  Z <- x$TData      # raw data after standardisation
  d <- x$SVD$d      # singular values
  U <- x$SVD$u      # left singular vectors (columns of U)
  V <- t(x$SVD$vt)  # right singular vectors (columns of V)
  
  #
  # Scaling, scores and loadings
  #
  ### scaling: alpha and beta values for different scaling options
  fact.scale <- switch(scale,
      singul = rep(0.5, 2),
      eigen = rep(1, 2),
      uvc = c(1,0),
      uvr = c(0,1))
  
  # Calculate weighted factor scores (for rows)
  S <- sweep(crossprod(t(as.matrix(sweep(Z, 2, sqrt(Wp), "*"))), V[, dim]), 2, 
             d[dim]^(fact.scale[1]-1), "*")
  # (RV: Does the next line do anything? Normally, the eigenvalues sum to 1)
  S <- S * sqrt(sum(x$eigen)) / sqrt(sum(x$eigen)^fact.scale[1])
  # Set columns of S that correspond to small singular values to 0 (dimensions close to singular)
  S[, d[dim] < 1E-6] <- 0
  
  # Calculate weighted factor loadings (for columns)
  L <- sweep(crossprod(as.matrix(sweep(Z, 1, sqrt(Wn), "*")), U[,dim]), 2, d[dim]^(fact.scale[2]-1), "*")
  # (RV: same remark as above for S)
  L <- L * sqrt(sum(x$eigen)) / sqrt(sum(x$eigen)^fact.scale[2])
  L[, d[dim] < 1E-6] <- 0
  #
  # Rotation
  # Default: flips both the X and Y axes
  # (RV: Generalised from ncol=2 to ncol=dim(S)[2] to account for more than 2 dimensions)
  S <- S * matrix(rot, ncol = ncol(S), nrow = nrow(S), byrow = TRUE)
  L <- L * matrix(rot, ncol = ncol(L), nrow = nrow(L), byrow = TRUE)
  
  #
  # zooming
  # (RV: it should not be allowed to zoom S and L differently if they are to be
  # projected on the same plane in the plot)
  ss <- S * zoom[1]
  ll <- L * zoom[2]
  
  #
  # Calculate which rows are far from the origin (most specific)
  #
  DS <- S[, 1]^2 + S[, 2]^2  # distance from center
  DL <- L[, 1]^2 + L[, 2]^2  # distance from center
  
  # Threshold for row labels
  if (label.tol > 1) # label.tol most distant rows are plotted as circles and labelled
    thres <- sort(DS, decreasing = TRUE)[min(length(DS), floor(label.tol))]
  else # label.tol percent most distant rows are plotted as circles and labelled
    thres <- if (label.tol == 0) Inf else quantile(DS, probs = 1-label.tol)
    
  # Threshold for row labels
  if (label.col.tol > 1) # label.tol most distant rows are plotted as circles and labelled
    colThres <- sort(DL, decreasing = TRUE)[min(length(DL), floor(label.col.tol))]
  else # label.tol percent most distant rows are plotted as circles and labelled
    colThres <- if (label.col.tol == 0) Inf else quantile(DL, probs = 1-label.col.tol)
  
  #
  # Positioning
  #
  is <- switch(show.row,
      all = rep(TRUE, nrow(S)), 
      position = x$pos.row)
  
  il <- switch(show.col,
      all = rep(TRUE, nrow(L)),
      position = x$pos.column)
  
  isel <- is & (DS >= thres)
  
  # Only draw a plot if the projection is 2D and the user requested a plot (default)
  if (do.plot && (length(dim) == 2)){
    
    #
    # Compute range of plot
    #
    xrange <- range(ll[!x$pos.column, 1], ss[!x$pos.row, 1], 0)
    yrange <- range(ll[!x$pos.column, 2], ss[!x$pos.row, 2], 0)
    xrange <- ifelse(rep(diff(xrange) == 0, 2),
        c(-1, 1), # Default range [-1,1] if calculated range is 0
        xrange + diff(xrange) * c(-0.1, 0.1)) # Else expand range by 10%
    yrange <- ifelse(rep(diff(yrange) == 0, 2),
        c(-1, 1),
        yrange + diff(yrange) * c(-0.1, 0.1))
                                        # RV: changed diff(xrange) to diff(yrange):
    # the eqscplot funtion (ratio=1) takes care of equally scaled axes, so make sure
    # we have the correct range for each axis here.
    
    #
    # Clipping
    # (replace coordinates outside the plotted range by the min or max coordinates)
    ss[, 1] <- pmin(pmax(ss[,1], xrange[1]), xrange[2])
    ss[, 2] <- pmin(pmax(ss[,2], yrange[1]), yrange[2])
    ll[, 1] <- pmin(pmax(ll[,1], xrange[1]), xrange[2])
    ll[, 2] <- pmin(pmax(ll[,2], yrange[1]), yrange[2])
    
    #
    # Set-up plot
    #
    opar <- par(pty = "m") # preserve configuration
    # Create a window with the maximal plotting region
    # Equal scale plot function from MASS library
    dotList <- list(...)
    if (is.null(dotList$sub)){ 
      dotList$sub <- paste("Closure = ", x$closure, ", Center = ", x$center, ", Norm. = ", 
          x$normal, ", Scale = ", scale, ", RW = ", x$row.weight, ", CW = ", x$col.weight, sep="")
      if (is.null(dotList$cex.sub))
        dotList$cex.sub <- 0.85 
    } 
    dotList$cex.sub <- if (is.null(dotList$cex.sub)) 0.85 else dotList$cex.sub 
    if (is.null(dotList$xlab))
      dotList$xlab <- paste("PC", dim[1], " ", 100 * round(x$contrib[dim[1]], 2), "%", sep = "")
    if (is.null(dotList$ylab))
      dotList$ylab <- paste("PC", dim[2], " ", 100 * round(x$contrib[dim[2]], 2), "%", sep = "")
    
    dotList$x <- xrange
    dotList$y <- yrange
    dotList$ratio <- 1
    dotList$tol <- 0
    dotList$type <- "n"
    dotList$axes <- FALSE
    dotList$cex.lab = 0.85
    
    do.call("eqscplot", dotList)

    # Scales
    # (RV: the drop parameter doesn't have an effect in the following lines)
    usr <- par("usr") # Retrieve extremes of coordinates in the plotting region
    sx <- diff(usr[c(1,2)]) / (25.4 * par("pin")[1, drop = TRUE]) # scale in mm
    sy <- diff(usr[c(3,4)]) / (25.4 * par("pin")[2, drop = TRUE])
    
    #
    # Plot rows close to 0,0 as unlabelled dots or as density maps
    #
    # Select rows to be plotted
    i <- is & (DS < thres)
    if (!do.smoothScatter) # plot inner points as unlabelled dots
      points(ss[i,1], ss[i,2], col = colors[1], cex = 0.825, lwd = 2)
    else # plot as density maps
      smoothScatter(x = ss[i,1], y = ss[i,2], nbin = 256, nrpoints = 0,
          add = TRUE, colramp = colorRampPalette(c("white", "burlywood"))) # add image to current eqscplot axes, instead of overwritting
    # Alternative to adding to the axes set up by eqscplot, use the following
    # parameters to set up smoothScatter's own axes.
    # xlim = xrange, ylim = yrange, # make sure we can add distant points later
    # asp=1, # equal scale axes
    # xaxt = "n", yaxt = "n", # do not plot the axes scales
    # xlab = paste("PC", dim[1], " ", 100 * round(x$contrib[dim[1]], 2), "%", sep = ""),
    # ylab = paste("PC", dim[2], " ", 100 * round(x$contrib[dim[2]], 2), "%", sep = ""),
    

    ### plot distant rows as circles with areas proportional to x$Rm
    sqs <- 0.5 * sx * pmax(0.02, row.size * sqrt(abs((x$Rm)/(max(abs(x$Rm))))))
    yoffset <- sy * (2 + sqs / sx)
    
    
    if (sum(isel) > 0){ # if there is at least 1 point to plot
      symbols(ss[isel, 1], ss[isel, 2], circle = sqs[isel], 
        inches = FALSE, lwd = 3, add = TRUE, fg = colors[2])
      if (is.null(labels)) labels <- x$row.names
      text(ss[isel, 1], ss[isel, 2] - yoffset[isel], adj = c(0.5, 1), 
        cex = lab.size, labels = labels[isel], # x$row.names[isel]
        col = colors[2])
    }
    
    
    
    ### plot columns with indication of column-grouping
    iGroup <- unique(col.group[il]) # unique groups in columns to be plotted
    for (i in 1:length(iGroup)){
      iColGroup <- col.group == iGroup[i] 
      ii <- il & iColGroup # Select columns in group i selected for plotting
      iiExtremes <- ii & (DL >= colThres) # in order to only label most extreme ones
      
      if (sum(ii) > 0){ # at least one point to plot
        
        if (col.areas){ # use squares with size (ie. area) proportional to x$Cm
          # RV: Use same scale formula as for rows
          sqs <- 0.5 * sx * pmax(0.02, col.size * sqrt(abs((x$Cm[ii])/(max(abs(x$Cm))))))
          yoffset <- sy * (5 + sqs / (2 * sx)) # ii based
          
          
          
          symbols(ll[ii, 1], ll[ii, 2],
              square = sqs, inches = FALSE, lwd = 3, add = TRUE, fg = colors[2+iGroup[i]])
          if (sum(iiExtremes) > 0){
            
            sqsExtremes <- 0.5 * sx * pmax(0.02, col.size * sqrt(abs((x$Cm[iiExtremes])/(max(abs(x$Cm))))))
            yoffsetExtremes <- sy * (5 + sqsExtremes / (2 * sx)) # ii based
            
            text(ll[iiExtremes,1], ll[iiExtremes,2] + yoffsetExtremes,
                adj = c(0.5, 1), cex=lab.size, labels=sampleLabels[iiExtremes], col=colors[2+iGroup[i]]) # x$col.names[ii]
          }
          
        } else { # Use different symbols, ignore size
          
          yoffset <- sy * (5 + col.size / 5)
          points(ll[ii, 1], ll[ii, 2],
              pch = col.symbols[iGroup[i]], col = colors[2+iGroup[i]],
              cex = col.size / (25.4 * par("csi")), lwd = 3)
          if (sum(iiExtremes) > 0)
            text(ll[iiExtremes, 1], ll[iiExtremes, 2] + yoffset[iiExtremes],
                adj = c(0.5, 1), cex = lab.size, labels = sampleLabels[iiExtremes], col = colors[2+iGroup[i]]) # x$col.names[ii]
          
        }
      }
    }
    
    ### finish plot, put cross on 0,0, box, and legend 
    lines(c(-2.5, 2.5) * sx, c(0,0), lwd = 3)
    lines(c(0,0), c(-2.5, 2.5) * sy, lwd = 3)
    box()
   
    par(opar) # Restore plotting configuration
  }
  
  #
  # Return value: list of coordinates and indication of most distant (most specific) points
  #
  cols <- as.data.frame(cbind(ss, isel))
  
  if (length(dim) == 1){
    dimnames(cols) <- list(x$row.names, c("X", "Select"))
    dimnames(ll) <- list(x$col.names, c("X"))
  } else if (length(dim) == 2){
    dimnames(cols) <- list(x$row.names, c("X", "Y", "Select"))
    dimnames(ll) <- list(x$col.names, c("X", "Y"))
  } else if (length(dim) == 3){
    dimnames(cols) <- list(x$row.names, c("X", "Y", "Z", "Select"))
    dimnames(ll) <- list(x$col.names, c("X", "Y", "Z"))
  } else { # More than 3 dimensions requested
    dimnames(cols) <- list(x$row.names, c(paste("Prf", dim, sep=""), "Select"))
    dimnames(ll) <- list(x$col.names, paste("Pcf", dim, sep=""))
  }
  r <- list(Rows = cols, Columns = as.data.frame(ll))
  class(r) <- "plot.mpm"
  invisible(r)
}
