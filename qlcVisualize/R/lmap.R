lmap <- function( points, data
            , main = NULL
						, draw = 5
						, levels = c(0.41, 0.46, 0.51)
            , labels = NULL
            , cex = 0.7
						, col = "rainbow"
						, add = FALSE
        # smoothing paramater for Krig
            , lambda = NA
        # parameters for legend
            , legend = TRUE
						, position = "bottomleft"
						, cex.legend = 0.7
						, font = ""
						, note = TRUE
        # write to file
						, file.out = NULL
				# parameters passed through to function "boundary"
						, ...
				) {

  # ============
	# preparations
  # ============

  if (is.data.frame(points)) {
    points <- as.matrix(points)
  }

  # treat data as matrix: points x levels
  if (is.null(dim(data))) {
    # make matrix from vector of forms
    data <- as.factor(data)
    data <- sapply(levels(data),function(x){as.numeric(data==x)})
	} else {
    # if a matrix, check for numeric
    data <- as.matrix(data)
    if (!is.numeric(data)) {
      stop("Matrix or dataframe can only contain numerical values")
    }
	}

  # check number of points and size of data
  if (nrow(points) != nrow(data)) {
    stop("Number of points does not correspond to the data supplied")
  }

  # points ignored because of no data
  ignore <- rowSums(data, na.rm = TRUE) == 0
  empty.points.present <- sum(ignore) > 0

  # normalize between 0 and 1: each point (row) adds up to 1
  sums <- rowSums(data, na.rm = TRUE)
  sums[sums == 0] <- 1
  data <- data/sums

  # check for multi-valued data
  single.valued.data <- sum(data[data != 1 & data != 0], na.rm = TRUE) == 0

  # =================================
  # which words to include in graphic
  # =================================

  if (is.numeric(draw) & length(draw) == 1) {
    # only most frequent levels are included
    freq <- colSums(data, na.rm = TRUE)
    ordered <- order(freq, decreasing = TRUE, na.last = NA)
    selection <- na.omit(ordered[1:draw])
  } else if (is.numeric(draw)) {
    # use numbers as column-numbers
    selection <- draw
  } else {
    # manually selected level-names in order selected
    selection <- sapply(draw,function(x){which(colnames(data)==x)})
  }

  # the rest is added as "other" in the last column
  others.present <- length(selection) < ncol(data)
  if (others.present) {
    other <- rowSums(data[ , -selection, drop = FALSE])
    data <- data[ , selection, drop = FALSE]
    data <- cbind(data, other = other)
  } else {
    data <- data[ , selection, drop = FALSE]
  }

  # ===========
  # set colours
  # ===========

  if (is.null(col)) {

    # grey when col = NULL

    if (single.valued.data) {
      col <- rep("black", times = length(selection))
      if (others.present) { col <- c(col, "grey") }
    } else {
      col <- grey( (0:length(selection)) / (length(selection)+1) )
      if (others.present) { col <- c(col, "white") }
    }

  } else {

    # either use in-built palette or manually specify colors

    palettes <- c("rainbow"
                  , "heat.colors"
                  , "terrain.colors"
                  , "topo.colors"
                  , "cm.colors")

    if (length(col) ==  1) {
      if (!is.na(pmatch(col, palettes))) {
        cols <- paste0(palettes[pmatch(col, palettes)]
                       , "(", length(selection), ")"
                       )
        col <- eval(parse(text = cols))
      }
    } else if (length(col) != length(selection)) {
      stop("Number of colors specified at 'col' should be the same as number of levels to plot specified at 'draw', or it should be a name of a built-in palette")
    }
    # add grey for other levels
    if (others.present) {
      col <- c(col, "grey")
    }
  }

  # determine interpolation area
  # add zero coordinates to get nicer plotted areas
  zeros <- boundary(points, plot = FALSE, ...)

	# ========
	# plotting
	# ========

	if (!add) {

	  # open plotting device for saving file
  	if (!is.null(file.out)) {
  		dev.new(width = 6
  			  	, height = 6
  				  , type = "pdf"
  				  , file = file.out
  		)
  	}

  	# prepare plotting frame
  	plot( rbind(points, zeros)
  	      , main = main
      		, type = "n"
      		, xlab = ""
      		, ylab = ""
      		, axes = FALSE
      		)
	}

	# =============
	# make contours
  # using package "fields"
	# =============

	for (i in 1:length(selection)) {

	  # determine height
	  h <- data[,i]
	  h0 <- rep.int(0, times = nrow(zeros))

	  # make countours
	  kriging <- fields::Krig(x = rbind(points, zeros)
	                  , Y = c(h, h0)
	                  , give.warnings = FALSE
	                  , lambda = lambda
	  )
	  surface <- fields::predictSurface(kriging)
	  contour(surface
	          , levels = levels
	          , lwd = rev(seq(2, 0.5, length.out = length(levels)))
	          , col = col[i]
	          , drawlabels = FALSE
	          , add = TRUE
	  )
	}

	# ==========
	# add legend
	# ==========

	if (legend) {
	  par(family = font)
	  if (single.valued.data & is.null(labels)) {

	    pch <- 1:length(selection)
	    if (others.present) { pch <- c(pch, 0) }

	    legend( position
	            , pch = pch
	            , legend = colnames(data)
	            , cex = cex.legend
	            , col = col
	    )
	  } else {
	    legend( position
	            , fill = col
	            , legend = colnames(data)
	            , cex = cex.legend
	    )
	  }
	  par(family = "")
	}

  # =====================
  # add note about levels
  # =====================

  if (note) {
    mtext(paste0("Levels drawn at "
                  , paste(100 * levels, collapse = "-")
                  , "%")
          , side = 1
          , cex = cex.legend
          , col = "grey"
          )
  }

	# =============================
  # plot labels or points or pies
	# =============================

  if (is.null(labels)) {
    # plot points/pies

    # add small grey points for points without data
    points( points[ignore, ,drop = FALSE]
            , pch = 20
            , col = "grey"
            , cex = cex
          )

    if (single.valued.data) {
      # plotting symbols if maximally one symbol per point
      for (i in 1:length(selection)) {
        points( points[data[,i] > 0, , drop = FALSE]
                , pch = i
                , cex = cex
                , col = col[i]
               )
      }
      if (ncol(data)>length(selection)) {
        points( points[data[,length(selection)+1] > 0, , drop = FALSE]
                , pch = 0
                , cex = cex
                , col = col[length(selection)+1])
      }
    } else {
      # otherwise plot pies. This implementation is slow!!!
      mapplots::draw.pie( x = points[!ignore,1]
                          , y = points[!ignore,2]
                          , z = data[!ignore,]
                          , radius = cex/20
                          , col = col
                          )
    }
  } else {
    # plot labels

    if (length(labels) != nrow(points)) {
      stop("Number of labels is not equal to the number of points")
    }
    bw <- rep("black", times = nrow(points))
    bw[ignore] <- "grey"

    par(family = font)
    text( points
          , labels = labels
          , col = bw
          , cex = cex
    )
    par(family = "")

  }

	# close plotting device when saving to file
	if (!is.null(file.out)) {
		dev.off()
	}
}
