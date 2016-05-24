#' Plot Method for Class 'int.lsbclust'
#' 
#' Two-dimensional plot method for object of class 'int.lsbclust' as output by \code{\link{int.lsbclust}}.
#' 
#' @aliases plot.int.lsbclust
#' @param x An object of class \code{int.lsbclust}.
#' @param which A vector indicating which item segments to plot.
#' @param plot.type Character string giving the type of plots to produce: either \code{"biplots"}
#' for the biplots approximating the cluster means, \code{"means"} for level plots of the cluster means 
#' themselves or \code{"estimates"} for level plots of the low-rank approximations of the cluster means 
#' (as represented in the biplots).
#' @param segments A logical vector with two elements,  indicating whether the rows and columns should 
#' be plotted as line segments or not.
#' @param biplot.axes A logical indicating whether to plot calibrated biplot axes for the line
#' segments indicated in \code{segments} or not.
#' @param nmarkers Either a single integer giving the number of desired markers per biplot axis 
#' for all axes, or a named list. This is passed as the argument \code{n} to \code{\link{pretty}}. See
#' \code{Details} for information on the list option.
# @param axes.col The colour used for the biplot axes.
#' @param alpha Numeric value in [0, 1] which determines how the singular values are distributed
#' between rows and columns. It will trigger a recomputation of the updates if it does not correspond
#' to the value used when fitting the model.
#' @param check.alpha Logical indicating whether to look for a better alpha. This is only used when
#' \code{alpha = NULL} is used.
#' @param fix.alpha Logical indicating whether to fix alpha across all clusters or not 
#' when \code{fixed == "none"}.
#' @param probs Argument passed to \code{\link{quantile}} to determine the alpha value. The 
#' corresponding quantile of the distances of all points in the biplots to the origin will be
#' used to determine alpha in case check.alpha = TRUE.
#' @param arrange Logical indicating whether to arrange the plots side-by-side
#' via \code{\link{grid.arrange}} or not.
#' @param fix.limits Logical indicating whether biplot x- and y-limits must be fixed across clusters
#' or not. Note that this is automatically set to \code{TRUE} when \code{fixed == "rows"} or 
#' \code{fixed == "columns"}. When limits are fixed, the axis calibrations are also turned off.
#' @param limit.exp A numeric expansion factor applied multiplicatively to the plot limits, but only
#' when \code{fixed} equals \code{"rows"} or \code{"columns"}.
#' @param lambda.scale Logical indicating whether to apply lambda scaling to the coordinates or not.
#' If true, the scaling is done such that the average squared distance to the origin is equal
#' for the row and column coordinates.
#' @param procrustes.rotation Logical indicating whether to do Procrustes rotations so that the
#' location of the axes indicated as segments (see argument \code{segments}) are similar
#' across configurations.
#' @param fix.lambda Logical indicating whether to fix lambda across all clusters or not.
#' @param labs.grey Logical indicating whether to apply greying to the text labels are well.
# @param labs.min.grey A numeric value between zero and one giving the lightest grey to be
# used for the labels. The value 0.10 corresponds to the default colour scale used for the 
# points and lines for the rows and columns.
#' @param label.0 Logical indicating whether to label the origin or not.
#' @param tick.length The required tick length as a \code{\link{unit}} object. It defaults to a 
#' propoprtion of the width of the plot region (through lazy evaluation).
#' @param axis.col The colour of the biplot axes.
#' @param label.size The size of the labels for the markers on the biplot axes.
#' @param axis.size Line size for biplot axes.
#' @param axis.title.size Size of biplot axis titles.
#' @param draw.axis A list with up to two components which must be named \code{"rows"} and 
#' \code{"columns"}. Each element contains a vector indicating which biplot axes should be drawn. 
#' The vectors can be character vectors containing the names of the axes to be drawn, numeric
#' vectors containing indices indicating which axes to draw, or logical vectors indicating which 
#' biplot axes to draw. In case of the default value \code{NULL}, the elements of \code{segments}
#' are used for the \code{"rows"} and \code{"columns"} entries.
#' @param points.col A named list containing the colours to use for plotting the sets of points. The 
#' elements \code{"rows"} and \code{"columns"} contain vectors giving the colours for the points. Single
#' element vectors are recycled across the different points, otherwise the vectors must be of the 
#' appropriate length.
#' @param offset.tick.labels A numeric value giving the offset factor of the biplot axis marker labels
#' from their respective tick marks. Higher (lower) values lead to labels being further from 
#' (nearer to) their respective tick marks.
#' @param offset.axis.title A names list of (up to) two numeric values giving the fixed length offset of the 
#' biplot axis title label from the end of the axis segment. The two elements must have names 
#' \code{"rows"} and code{"columns"}.
#' @param axis.arrow An \code{\link{arrow}} object to be used for the endpoints of biplot axis 
#' segment lines. This is passed to \code{\link{geom_segment}}.
#' @param \dots Additional arguments passed to \code{\link{theme}}.
#' @details
#' In case \code{nmarkers} is a list, it can have up to two elements. These are required to be named 
#' \code{"rows"} and/or \code{"columns"}, otherwise an error will be thrown. The elements of the list
#' contains either single numeric values each or numeric vectors of the appropriate lengths
#' indicating the \code{n} argument passed to \code{\link{pretty}}.
#' @keywords hplot
#' @method plot int.lsbclust
#' @export
plot.int.lsbclust <- function(x, which = seq_len(nclust), plot.type = c("biplots", "means", "estimates"), 
                              segments = NULL,  biplot.axes = TRUE, nmarkers = 5, alpha = NULL, 
                              check.alpha = TRUE, fix.alpha = FALSE, probs = 0, arrange = FALSE, 
                              fix.limits = TRUE, limit.exp = 1.05, lambda.scale  = TRUE, 
                              procrustes.rotation = x$fixed == "none", 
                              fix.lambda = FALSE, labs.grey = TRUE, #labs.min.grey = 0.1,
                              label.0 = FALSE, tick.length = 0.0075 * diff(lims), axis.col = "grey60", 
                              label.size = 3, axis.size = 0.25, axis.title.size = 4, 
                              draw.axis = NULL, points.col = list(rows = "red", columns = "blue2"), 
                              offset.tick.labels = 3.5, 
                              offset.axis.title = list(rows = 0.015 * max(nchar(rnms)), 
                                                       columns = 0.015 * max(nchar(cnms))),
                              axis.arrow = grid::arrow(angle = 20, length = grid::unit(0.0175, "npc")), ...) {
  # axes.col = "grey50"
  ## Avoid globalVariable issues
  y <- Fit <- Fill <- Row <- Column <- Value <- Label <- xupper <- xlower <- yupper <- ylower <- xlabs <- ylabs <- NULL
  
  ## Determine ndim and number of clusters
  nC <- length(x$C)
  nD <- length(x$D)
  nclust <- max(nC, nD)
  J <- nrow(x$C[[1]])
  K <- nrow(x$D[[1]])
  nvec <- table(x$cluster)
  N <- sum(nvec)
  
  ## Check 'segments' and set to c(FALSE, TRUE) unless nC == 1
  if (is.null(segments))  {
    segments <- switch(x$fixed, "rows" = c(TRUE, FALSE), "columns" = c(FALSE, TRUE), 
                       "none" = c(FALSE, TRUE))
  } else {
    if(length(segments) != 2) stop("Argument 'segments' must be of length 2.")
  }
  
  ## Check ndim
  ndim <- ncol(x$C[[1]])
  if(ndim != 2) stop("Only plot for 2 dimensions are currently implemented.")
  
  ## Check alpha and recalculate Cs and Ds if not the same
  if (!is.null(alpha) && check.alpha) {
    check.alpha <- FALSE
    message("check.alpha set to FALSE since explicit alpha value was supplied.")
  }
  if (is.null(alpha)) alpha <- x$alpha
  if (alpha < 0 || alpha > 1) stop("Alpha must be between 0 and 1.")
  
  ## Check type
  plot.type <- tolower(plot.type)
  plot.type <- match.arg(arg = plot.type, several.ok = FALSE)
  
  ## Set up list containing plot grobs
  plots <- vector(length = max(nC, nD), mode = "list")
  
  ## Do the biplots if required
  if ("biplots" %in% plot.type) {
    ## Get names for rows and columns
    rnms <- rownames(x$C[[1]])
    cnms <- rownames(x$D[[1]])
    
    ## Get indicator for segments to switch on, get total number of segment axes to be drawn
    segind <- paste(as.numeric(segments), collapse = "")
    seglen <- switch(segind, "00" = 0, "01" = K, "10" = J, "11" = J + K)
    
    ## Handle nmarkers
    if (is.list(nmarkers)) {
      
      ## Check names
      names(nmarkers) <- tolower(names(nmarkers))
      if (!all(names(nmarkers) %in% c("rows", "columns"))) 
        stop("Argument 'nmarkers' have incorrect element names.")
      
      ## Check type
      if (!all(sapply(nmarkers, is.numeric))) stop("Argument 'nmarkers' contains an element of incorrect type.")
      
      ## Handle NULL
      if (is.null(nmarkers[["rows"]])) nmarkers[["rows"]] <- rep(0, J)
      if (is.null(nmarkers[["columns"]])) nmarkers[["columns"]] <- rep(0, K)
      
      ## Expand if any of the lengths are 1
      if (length(nmarkers[["rows"]]) == 1) nmarkers[["rows"]] <- rep(nmarkers[["rows"]], J)
      if (length(nmarkers[["columns"]]) == 1) nmarkers[["columns"]] <- rep(nmarkers[["columns"]], K)
      
      ## Check lengths
      if (any(sapply(nmarkers, length) != c(J, K))) stop("Argument 'nmarkers' contains an element of incorrect length.")
    }
    if (is.numeric(nmarkers)) nmarkers <- list(rows = rep(nmarkers, J), columns = rep(nmarkers, K))
    
    ## Check draw.axis type and make sure it is a named list with logical vectors as elements
    if (is.null(draw.axis)) 
      draw.axis <- list(rows = rep(segments[1], J), columns = rep(segments[2], K))
    else {
      if (is.list(draw.axis)) {
        
        ## Check names
        names(draw.axis) <- tolower(names(draw.axis))
        if (!all(names(draw.axis) %in% c("rows", "columns"))) 
          stop("Argument 'draw.axis' have incorrect element names: the names 'rows' and/or 'columns' must be used.")
        
        ## Check row type if not logical already
        if (is.numeric(draw.axis$rows)) {
          tmp <- rep(FALSE, J)
          tmp[draw.axis$rows] <- TRUE
          draw.axis$rows <- tmp
        }
        if (is.character(draw.axis$rows)) {
          draw.axis$rows <- rnms %in% draw.axis$rows
        }
        if (is.null(draw.axis$rows)) {
          draw.axis$rows <- rep(FALSE, J)
        }
        
        ## Check column type if not logical already
        if (is.numeric(draw.axis$columns)) {
          tmp <- rep(FALSE, K)
          tmp[draw.axis$columns] <- TRUE
          draw.axis$columns <- tmp
        }
        if (is.character(draw.axis$columns)) {
          draw.axis$columns <- cnms %in% draw.axis$columns
        }
        if (is.null(draw.axis$columns)) {
          draw.axis$columns <- rep(FALSE, K)
        }
      }
      if (is.logical(draw.axis) && length(draw.axis) == 1) 
        draw.axis <- list(rows = rep(draw.axis, J), columns = rep(draw.axis, K))
    }
    
    ## Check draw.axis length
    if (!all(sapply(draw.axis[c("rows", "columns")], length) %in% c(J, K))) stop("Argument 'draw.axis' is of incorrect length.")
    
    ## Check offset.axis.title
    names(offset.axis.title) <- tolower(names(offset.axis.title))
    if (!all(names(offset.axis.title) %in% c("rows", "columns"))) 
      stop("Argument 'offset.axis.title' have incorrect element names.")
    
    ## Handle points.col
    ## Check names
    names(points.col) <- tolower(names(points.col))
    if (!all(names(points.col) %in% c("rows", "columns"))) 
      stop("Argument 'points.col' have incorrect element names.")
    
    ## Expand if any of the lengths are 1
    if (length(points.col[["rows"]]) == 1) points.col[["rows"]] <- rep(points.col[["rows"]], J)
    if (length(points.col[["columns"]]) == 1) points.col[["columns"]] <- rep(points.col[["columns"]], K)
    
    ## Check lengths
    if (!segments[1] && length(points.col$rows) != J)
      stop("Element 'rows' of argument 'points.col' is of incorrect length.")
    if (!segments[2] && length(points.col$columns) != K)
      stop("Element 'columns' of argument 'points.col' is of incorrect length.")
    
    ## Recalculate Cs and Ds if necessary 
    if(alpha != x$alpha || check.alpha){
      
      ## For fixed = "rows", recalculating function
      reRows <- function(alpha, index = NULL) {
        svdX <- x$svd[[1]]
        
        ## Determine updates from SVD
        Cs <- list(svdX$u %*% diag(svdX$d[1:ndim]^alpha))
        Dstar <- svdX$v %*% diag(svdX$d[1:ndim]^(1 - alpha))
        
        ## Rescale Dstar and divide up in Du's
        Dstar.sc <- Dstar * 1/sqrt(nvec) %x% matrix(1, K, ndim)
        Ds <- split.data.frame(Dstar.sc, rep(1:nclust, each = K))
        return(list(Cs = Cs, Ds = Ds, Dstar = Dstar))        
      }
      
      ## For fixed = "columns", recalculating function
      
      reCols <- function(alpha, index = NULL) {
        svdX <- x$svd[[1]]
        
        ## Determine updates from SVD
        Cstar <- svdX$u %*% diag(svdX$d[1:ndim]^alpha)
        Ds <- list(svdX$v %*% diag(svdX$d[1:ndim]^(1 - alpha)))
        
        ## Rescale Cstar and divide up in Cu's
        Cstar.sc <- Cstar * 1/sqrt(nvec) %x% matrix(1, J, ndim)
        Cs <- split.data.frame(Cstar.sc, rep(1:nclust, each = J))
        return(list(Cs = Cs, Ds = Ds, Cstar = Cstar))  
      }
      
      ## For fixed = "none", recalculating function
      reNone <- function(alpha, index = NULL) {
        if(is.null(index)) {
          
          ## Same alpha for all clusters
          svdX <- x$svd
          
          ## Calculate Cu's and Du's
          Cs <- lapply(svdX, function(x) x$u %*% diag(x$d[1:ndim]^alpha))
          Ds <- lapply(svdX, function(x) x$v %*% diag(x$d[1:ndim]^(1 - alpha)))
          return(list(Cs = Cs, Ds = Ds))  
        } else {
          
          ## Different alpha for each cluster
          svdX <- x$svd[[index]]
          
          ## Calculate Cu's and Du's
          C <- svdX$u %*% diag(svdX$d[1:ndim]^alpha)
          D <- svdX$v %*% diag(svdX$d[1:ndim]^(1 - alpha))
          return(list(Cs = C, Ds = D))  
        }
      }
      
      ## Switch for recalculating functions
      reFun <- switch(x$fixed, rows = reRows, columns = reCols, none = reNone)
      
      ## Optimization for alpha
      if(check.alpha) {
        
        ## Optimization function
        optAlpha <- function(alpha, index) {
          
          updCD <- reFun(alpha, index)
          
          ## Calculate distances
          if(is.null(index)){
            rowdist <- sqrt(sapply(updCD$Cs, function(x) rowSums(x^2)))
            coldist <- sqrt(sapply(updCD$Ds, function(x) rowSums(x^2)))
          } else {
            rowdist <- sqrt(rowSums(updCD$Cs^2))
            coldist <- sqrt(rowSums(updCD$Cs^2))
          }
          
          return(quantile(c(rowdist, coldist), probs = probs))
        }
        if(x$fixed == "none" && !fix.alpha){
          opt <- lapply(seq_len(nclust), function(x) optimize(f = optAlpha, interval = c(0, 1), 
                                                              index = x, maximum = TRUE))
          alpha <- sapply(opt, "[[", "maximum")
          message("alpha chosen as ", paste(round(alpha, digits = min(4, getOption("digits"))), collapse = " "))
          
          ## Replace old Cs and Ds
          updCD <- lapply(seq_len(nclust), function(x) reFun(alpha[x], index = x))
          x$C <- lapply(updCD, "[[", "Cs")
          x$D <- lapply(updCD, "[[", "Ds")
          
        } else {
          opt <- optimize(f = optAlpha, interval = c(0, 1), index = NULL, maximum = TRUE)
          alpha <- opt$maximum
          message("alpha chosen as ", paste(round(alpha, digits = min(4, getOption("digits"))), collapse = " "))
          
          ## Replace old Cs and Ds
          updCD <- reFun(alpha, index = NULL)
          x$C <- updCD$Cs
          x$D <- updCD$Ds
          
          ## Replace Cstar / Dstar
          if(x$fixed == "rows") x$Dstar <- updCD$Dstar
          if(x$fixed == "columns") x$Cstar <- updCD$Cstar
        }
      } else {
        ## If only alpha should only be recalculated
        ## Replace old Cs and Ds
        updCD <- reFun(alpha, index = NULL)
        x$C <- updCD$Cs
        x$D <- updCD$Ds
        
        ## Replace Cstar / Dstar
        if(x$fixed == "rows") x$Dstar <- updCD$Dstar
        if(x$fixed == "columns") x$Cstar <- updCD$Cstar
      }
    }
    
    ## Do lambda-scaling if required
    if(lambda.scale) {
      
      ## ... when fixed == "none"
      if(x$fixed == "none") {
        
        ## Squared Euclidean distances to the origin
        ssC <- sapply(x$C, function(x) sum(x^2))
        ssD <- sapply(x$D, function(x) sum(x^2))
        
        ## Use sum of distances across all clusters if fix.lambda is TRUE
        if(fix.lambda) {
          lambda <- (J * sum(ssD) / (K * sum(ssC)))^0.25
        } else {
          lambda <- (J * ssD / (K * ssC))^0.25
        }
        
        ## Recalculate C's and D's
        x$C <- Map("*", x$C, as.list(lambda))
        x$D <- Map("/", x$D, as.list(lambda))
        message("lambda chosen as ", paste(round(lambda, digits = min(4, getOption("digits"))), collapse = " "))
      }
      
      ## ... when fixed == "columns"
      if(x$fixed == "columns") {
        ## Use sum of squared distances within cluster for C's (instead of Cstar)
        ssC <- sapply(x$C, function(x) sum(x^2))
        ssD <- sum(x$D[[1]]^2)
        lambda <- ((J * nclust * ssD) / (K * sum(ssC)))^0.25
        
        ## Recalculate C's and D's
        x$C <- Map("*", x$C, as.list(lambda))
        x$D[[1]] <- x$D[[1]]/lambda
        message("lambda chosen as ", round(lambda, digits = min(4, getOption("digits"))))
      }
      
      ## ... when fixed == "rows"
      if(x$fixed == "rows") {
        ## Use sum of squared distances within cluster for D's (instead of Dstar)
        ssC <- sum(x$C[[1]]^2)
        ssD <- sapply(x$D, function(x) sum(x^2))
        lambda <- ((J * sum(ssD)) / (K * nclust * ssC))^0.25
        
        ## Recalculate C's and D's
        x$C[[1]] <- x$C[[1]]*lambda
        x$D <- Map("/", x$D, list(lambda))
        message("lambda chosen as ", round(lambda, digits = min(4, getOption("digits"))))
      }
    }
    
    ## If nC != nD, and nC or nD == 1, expand the shorter so that they have equal length
    if(nC != nD){
      if(min(nC, nD) != 1) stop("Cannot interpret the number of C and/or D matrices.")
      if(nC == 1) x$C <- replicate(nD, x$C)
      if(nD == 1) x$D <- replicate(nC, x$D)
    }
    
    ## Do generalized Procrustes rotation if required
    if (procrustes.rotation) {
      if (segments[1] && !segments[2]) {
        
        ## Align C's and counterrotate the D's
        proc <- genproc(configs = x$C)
        x$C <- proc$rotated
        x$D <- Map("%*%", x$D, proc$rotations)
        message("Generalized Procrustes rotations applied to the rows.")
      }
      
      if (!segments[1] && segments[2]) {
        
        ## Align D's and counterrotate the C's
        proc <- genproc(configs = x$D)
        x$D <- proc$rotated
        x$C <- Map("%*%", x$C, proc$rotations)
        message("Generalized Procrustes rotations applied to the columns.")
      }
      
      if (all(segments)) 
        message("Generalized Procrustes rotation cannot be applied with all(segments) (simultaneously for the rows and colums).")
    }
    
    ## Data frame for Cs
    dfC <- as.data.frame(do.call(rbind, x$C), row.names = "")
    colnames(dfC) <- c("x", "y")
    dfC$Cluster <- rep(1:max(nC, nD), each = J)
    dfC$Fit <- c(x$rfit)
    dfC$Label <- rep(rnms, nclust)
    
    #   ## Add grey-scale colours for text
    #   if (labs.grey) dfC$Colour <- grey((1 - labs.min.grey) * (1 - dfC$Fit))
    #   else dfC$Colour <- rep("#000000", nrow(dfC))
    
    ## Data frame for Ds
    dfD <- as.data.frame(do.call(rbind, x$D), row.names = "")
    colnames(dfD) <- c("x", "y")
    dfD$Cluster <- rep(1:max(nC, nD), each = K)
    dfD$Fit <- c(x$cfit)
    dfD$Label <- rep(cnms, nclust)
    
    #   ## Add grey-scale colours for text
    #   if (labs.grey) dfD$Colour <- grey((1 - labs.min.grey) * (1 - dfD$Fit))
    #   else dfD$Colour <- rep("#000000", nrow(dfD))
    
    ## Joint limits for when one of nC or nD == 1
    lims <- limit.exp * range(dfC[, c("x", "y")], dfD[, c("x", "y")])
    
    ## Function to calculate relevant info per axis
    ## Arguments: x and y coordinates (scalars), vector of plot limits (length 2), number of markers, and axis name
    biplotaxis <- function(x, y, nmarkers, name, lims) {
      ## Get angle
      angle <- atan2(y, x)
      angle.orth <- pi/2 + angle
      
      ## Squared vector length
      len2 <- x^2 + y^2
      len <- sqrt(len2)
      
      ## Angles to plot corners
      lims.angles <- atan2(y = rep(lims, each = 2), x = lims[c(1, 2, 2, 1)])
      
      ## Determine which border is crossed by positive axis (0 = left, 1 = bottom, 2 = right, 3 = top)
      bord <- findInterval(angle, lims.angles)
      bord[bord == 4] <- 0
      
      ## Determine maximum inner product on axis, using the figure region range in either direction
      if (bord == 0 || bord == 2) {
        extremes <- c(lims[1] * (x + y^2 / x), lims[2] * (x + y^2 / x))
        maxinprod <- max(extremes)
        mininprod <- min(extremes)
      } else if (bord == 1 || bord == 3) {
        extremes <- c(lims[1] * (x^2 / y + y), lims[2] * (x^2 / y + y))
        maxinprod <- max(extremes)
        mininprod <- min(extremes)
      }
      
      ## Calculate markers
      markers <- pretty(c(mininprod, maxinprod), n = nmarkers)
      if (!label.0) markers <- markers[markers != 0]
      
      ## Get marker coordinates along the axis
      coordx <- x * markers / len2
      coordy <- y * markers / len2
      
      ## Get upper and lower coordinates for markers
      xupper <- coordx + tick.length * cos(angle.orth)
      xlower <- coordx - tick.length * cos(angle.orth)
      yupper <- coordy + tick.length * sin(angle.orth)
      ylower <- coordy - tick.length * sin(angle.orth)
      
      ## Get coordinates for marker labels
      ## Factor offset.tick.labels determines distance to marker
      if (bord == 0 || bord == 2) {
        xlabs <- coordx + offset.tick.labels * tick.length * cos(angle.orth)
        ylabs <- coordy + offset.tick.labels * tick.length * sin(angle.orth)
      } else if (bord == 1 || bord == 3) {
        xlabs <- coordx - offset.tick.labels * tick.length * cos(angle.orth)
        ylabs <- coordy - offset.tick.labels * tick.length * sin(angle.orth)
      }
      
      ## Calculate hjust and vjust
      hjust <- rep(0.5 + cos(angle.orth), length(coordx))
      vjust <- rep(0.5 + sin(angle.orth), length(coordy))
      
      ## Return
      out <- data.frame(x = coordx, y = coordy, axis = rep(name, length(coordx)),
                        xupper = xupper, xlower = xlower, yupper = yupper, ylower = ylower,
                        labels = markers, hjust = hjust, vjust = vjust, xlabs = xlabs, ylabs = ylabs)
      return(out)
    }
    
    ## Create all plots
    for (i in which) {
      
      ## Get subset of data
      dfC.cur <- dfC[dfC$Cluster == i, ]
      dfD.cur <- dfD[dfD$Cluster == i, ]
      
      ## Add points colours
      dfC.cur$Fill <- factor(seq_len(J))
      dfD.cur$Fill <- factor(seq_len(K))
      
      ## Do all with ggplot    
      plots[[i]] <- ggplot(data = dfC.cur, mapping = aes(x = x, y = y)) 
      
      #     ## Add lines for the origin
      #     plots[[i]] <- plots[[i]] + geom_hline(yintercept = 0, col = "lightgrey", linetype = 2) + 
      #         geom_vline(xintercept = 0, col = "lightgrey", linetype = 2)
      
      ## Add biplot axes, one at a time
      if (biplot.axes) {
        
        ## Biplot axis calculations for rows
        if (any(draw.axis$rows)) {
          dfMarkers.rows <- do.call(rbind,
                                    Map(biplotaxis, x = dfC.cur$x[draw.axis$rows], y = dfC.cur$y[draw.axis$rows], 
                                        nmarkers = nmarkers$rows[draw.axis$rows], name = rnms[draw.axis$rows],
                                        MoreArgs = list(lims = lims)))
        }
        
        ## Biplot axis calculations for columns
        if (any(draw.axis$columns)) {
          dfMarkers.columns <- do.call(rbind,
                                       Map(biplotaxis, x = dfD.cur$x[draw.axis$columns], y = dfD.cur$y[draw.axis$columns], 
                                           nmarkers = nmarkers$columns[draw.axis$columns], name = cnms[draw.axis$columns], 
                                           MoreArgs = list(lims = lims)))
        }
        if (any(draw.axis$rows) && any(draw.axis$columns)) dfMarkers <- rbind(dfMarkers.rows, dfMarkers.columns)
        if (any(draw.axis$rows) && !any(draw.axis$columns)) dfMarkers <- dfMarkers.rows
        if (!any(draw.axis$rows) && any(draw.axis$columns)) dfMarkers <- dfMarkers.columns
        
        ## Add biplot axis ablines
        if (any(draw.axis$rows)) {
          plots[[i]] <- plots[[i]] + 
            geom_abline(data = dfC.cur[draw.axis$rows, ], mapping = aes(intercept = 0, slope = y / x), 
                        colour = axis.col, size = axis.size)
        }
        if (any(draw.axis$columns)) {
          plots[[i]] <- plots[[i]] + 
            geom_abline(data = dfD.cur[draw.axis$columns, ], 
                        mapping = aes(intercept = 0, slope = y / x), colour = axis.col, size = axis.size)
        }
        
        if (any(unlist(draw.axis))) {
          ## Add markers to biplot axes
          plots[[i]] <- plots[[i]] + 
            geom_segment(data = dfMarkers, 
                         mapping = aes(x = xupper, xend = xlower, y = yupper, yend = ylower), 
                         colour = axis.col, size = axis.size)   
          
          ## Add marker labels to biplot axes
          plots[[i]] <- plots[[i]] + 
            geom_text(data = dfMarkers, 
                      mapping = aes(label = labels, x = xlabs, y = ylabs), # hjust = hjust, vjust = vjust
                      size = label.size, colour = axis.col)
        }
      }
      
      ## Add arrows for Cs and / or Ds, and labels
      if (segments[1]) {
        angles <- atan2(dfC.cur$y, dfC.cur$x)
        lens <- sqrt(rowSums(dfC.cur[, c("x", "y")]^2))
        
        ## Coordinates for axis labels
        xtitle <- (lens + offset.axis.title$rows) * cos(angles)
        ytitle <- (lens + offset.axis.title$rows) * sin(angles)
        if (labs.grey) dfTitle <- data.frame(x = xtitle, y = ytitle, Fit = dfC.cur$Fit, Label = rnms)
        else dfTitle <- data.frame(x = xtitle, y = ytitle, Fit = 1, Label = rnms)
        
        plots[[i]] <- plots[[i]] + 
          geom_segment(aes(xend = x, yend = y, x = 0, y = 0, alpha = Fit), 
                       arrow = axis.arrow) + 
          geom_text(data = dfTitle, mapping = aes(alpha = Fit, label = Label), size = axis.title.size)
      } 
      if (segments[2]) {
        angles <- atan2(dfD.cur$y, dfD.cur$x)
        lens <- sqrt(rowSums(dfD.cur[, c("x", "y")]^2))
        
        ## Coordinates for axis labels
        xtitle <- (lens + offset.axis.title$columns) * cos(angles)
        ytitle <- (lens + offset.axis.title$columns) * sin(angles)
        
        ## Add fit for alpha mapping if (labs.grey)
        if (labs.grey) dfTitle <- data.frame(x = xtitle, y = ytitle, Fit = dfD.cur$Fit, Label = cnms)
        else dfTitle <- data.frame(x = xtitle, y = ytitle, Fit = 1, Label = cnms)
        
        plots[[i]] <- plots[[i]] +
          geom_segment(data = dfD.cur, aes(xend = x, yend = y, x = 0, y = 0, alpha = Fit),  
                       arrow = axis.arrow) +
          geom_text(data = dfTitle, aes(alpha = Fit, label = Label), size = axis.title.size)
      }
      
      ## Add points for Cs and / or Ds, and labels
      if (!segments[1]) {
        plots[[i]] <- plots[[i]] + 
          geom_point(mapping = aes(alpha = Fit, fill = Fill), size = 4, shape = 21) +
          scale_fill_manual(values = points.col$rows, guide = FALSE) + 
          geom_text(mapping = aes(alpha = Fit, label = Label), vjust = 1.75, size = 4)
      }
      if (!segments[2]) {
        plots[[i]] <- plots[[i]] +
          geom_point(data = dfD.cur, mapping = aes(alpha = Fit, fill = Fill), size = 4, shape = 21) + 
          scale_fill_manual(values = points.col$columns, guide = FALSE) + 
          geom_text(data = dfD.cur, mapping = aes(alpha = Fit, label = Label), vjust = -0.75, size = 4)
      }
      
      ## Set theme elements, and add title
      plots[[i]] <- plots[[i]] + theme_bw(base_size = 16) + 
        theme(panel.grid = element_blank(), axis.title = element_blank()) +
        ggtitle(paste0("Interactions ", i, " (", round(100 * nvec[i] / N, 1),  "%)")) 
      
      ## Fix limits for all plots if fixed = "rows" or "columns" and hide axes
      if(fix.limits || x$fixed != "none"){
        plots[[i]] <- plots[[i]] + coord_fixed(xlim = lims, ylim = lims) + 
          theme(axis.ticks = element_blank(), axis.text = element_blank())
      } else {
        ## Otherwise apply limit.exp for expansion and keep x- and y-axes
        lims <- limit.exp * range(dfC.cur[, c("x", "y")], dfD.cur[, c("x", "y")])
        plots[[i]] <- plots[[i]]  + coord_fixed(xlim = lims, ylim = lims)
      }
      
      ## Set alpha to range from be mapped from 0 to 1 and pass additional arguments to theme
      plots[[i]] <- plots[[i]] + scale_alpha_continuous(limits = c(0, 1)) + theme(...)
    }
  }
  
  ## Do level plots of cluster means  
  if ("means" %in% plot.type) {
    
    ## Set up data frame for means
    dfMeans <- do.call(rbind, x$means)
    dfMeans <- as.data.frame(reshape2::melt(dfMeans, as.is = TRUE))
    colnames(dfMeans) <- c("Row", "Column", "Value")
    
    ## Change/add factor and ensure correct order
    dfMeans$Row <- factor(dfMeans$Row, levels = rev(rownames(x$means[[1]])))
    dfMeans$Column <- factor(dfMeans$Column, levels = colnames(x$means[[1]]))
    dfMeans$Cluster <- rep(rep(seq_len(nclust), each = J), K)
    
    ## Get limits for colour range (symmetric)
    lims <- c(-1, 1) * max(abs(dfMeans$Value))
    
    ## Do plots
    for (i in which) {
      plots[[i]] <- ggplot(data = dfMeans[dfMeans$Cluster == i, ], aes(y = Row, x = Column, fill = Value)) + 
        geom_tile(colour = "white") + 
        scale_fill_gradient2(limits = lims) +
        ggtitle(paste0("Mean: Interactions ", i, " (", round(100 * nvec[i] / N, 1),  "%)")) +
        theme_bw() + theme(...)
    }
  }
  
  ## Do level plots of cluster means  
  if ("estimates" %in% plot.type) {
    
    ## Set up data frame for means
    dfMeans <- do.call(rbind, Map(tcrossprod, x$C, x$D))
    dfMeans <- as.data.frame(reshape2::melt(dfMeans, as.is = TRUE))
    colnames(dfMeans) <- c("Row", "Column", "Value")
    
    ## Change/add factor and ensure correct order
    dfMeans$Row <- factor(dfMeans$Row, levels = rev(rownames(x$C[[1]])))
    dfMeans$Column <- factor(dfMeans$Column, levels = rownames(x$D[[1]]))
    dfMeans$Cluster <- rep(rep(seq_len(nclust), each = J), K)
    
    ## Get limits for colour range (symmetric)
    lims <- c(-1, 1) * max(abs(dfMeans$Value))
    
    ## Do plots
    for (i in which) {
      plots[[i]] <- ggplot(data = dfMeans[dfMeans$Cluster == i, ], aes(y = Row, x = Column, fill = Value)) + 
        geom_tile(colour = "white") + 
        scale_fill_gradient2(limits = lims) +
        ggtitle(paste0("Est. Mean: Interactions ", i, " (", round(100 * nvec[i] / N, 1),  "%)")) +
        theme_bw() + theme(...)
    }
  }
  
  ## Set plot names
  names(plots) <- paste0("plot", seq_len(nclust))
  
  ## Return plots
  if (arrange) do.call(gridExtra::grid.arrange, plots[which])
  else return(plots[which])
}

# if(biplot.axes) {
  #         if(x$fixed == "columns") {
  #           ## Determine markers for inner products
  #           inprods <- Map(tcrossprod, x$C, x$D)
  #           inprods <- do.call(rbind, inprods)
  #           markers <- apply(inprods, 2, pretty, n = nmarkers)
  #           
  #           ## Squared lengths of axes vectors
  #           lngs <- rowSums(dfD[, 1:2]^2)
  #           
  #           ## Get angles of all axes (plotting region is rectangular)
  #           angles <- atan2(dfD$y[seq_len(K)], dfD$x[seq_len(K)])
  #           lims.angles <- atan2(y = rep(lims, each = 2), x = lims[c(1, 2, 2, 1)])
  #           
  #           ## Determine which border is crossed by positive axis (1 = bottom, 2 = right, 3 = top, 4 = left)
  #           bords <- findInterval(angles, lims.angles)
  #           bords[bords == 0] <- 4
  #           
  #           ## Get coordinates on border for all axes
  #           dfLabs <- matrix(NA, nrow = K, ncol = 2)
  #           colnames(dfLabs) <- c("x", "y")
  #           rownames(dfLabs) <- cnms
  #           
  #           ## For border 1
  #           dfLabs[bords == 1, ] <- cbind(lims[1] / tan(angles[bords == 1]), lims[1])
  #           ## For border 2
  #           dfLabs[bords == 2, ] <- cbind(lims[2], lims[2] * tan(angles[bords == 2]))
  #           ## For border 3
  #           dfLabs[bords == 3, ] <- cbind(lims[2] / tan(angles[bords == 3]), lims[2])
  #           ## For border 4
  #           dfLabs[bords == 4, ] <- cbind(lims[1], lims[1] * tan(angles[bords == 4]))
  #           dfLabs <- as.data.frame(dfLabs)
  #           
  #           ## Axis labels: vjust and hjust according to location
  #           dfLabs$hjust <- 0.5
  #           dfLabs$vjust <- 0.5
  #           dfLabs$hjust[bords == 2] <- 1 #0
  #           dfLabs$hjust[bords == 4] <- 0 #1
  #           dfLabs$vjust[bords == 1] <- 0 #1
  #           dfLabs$vjust[bords == 3] <- 1 #0
  #           
  #           ## Coordinates for markers on each variable
  #           dfMarkers <- Map(f = function(x, y) outer(x, y, FUN = "/"), markers, as.list(lngs[dfD$Cluster == i]))
  #           dfMarkers <- Map(f = function(x, y) cbind(x, x) * matrix(unlist(y), nrow = nrow(x), ncol = 2, byrow = TRUE),
  #                            dfMarkers, split.data.frame(dfD[dfD$Cluster == i, c("x", "y")], f = seq_len(K)))
  #           dfMarkers <- as.data.frame(do.call(rbind, dfMarkers))
  #           colnames(dfMarkers) <- c("x", "y")
  #           dfMarkers$label <- unlist(markers)
  #         }
  #         
  #         ## Draw biplot axes
  #         plots[[i]] <- plots[[i]] + geom_abline(data = dfD[dfD$Cluster == i, ], aes(slope = y/x), 
  #                                                colour = axes.col) +
  #           geom_point(data = dfMarkers, colour = axes.col) + 
  #           geom_text(data = dfMarkers, aes(label = label), hjust = 1.5, size = 4, colour = axes.col)
  #           
  #         
  #         ## Add axis names
  #         for(j in seq_along(cnms)) {
  #           plots[[i]] <- plots[[i]] + 
  #             annotation_custom(grob = textGrob(label = cnms[j], hjust = dfLabs$hjust[j], vjust = dfLabs$vjust[j]), 
  #                               ymin = dfLabs$y[j], 
  #                               ymax = dfLabs$y[j], xmin = dfLabs$x[j], xmax = dfLabs$x[j])
  #         }
# }