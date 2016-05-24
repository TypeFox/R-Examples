#####################################
## Produce a 3D mosaic plot using rgl
#####################################

# TODO: provide formula interface
# TODO: handle zero margins (causes display to be erased in shapelist3d)
# DONE: handle zero cells 
# DONE: generalize the calculation of residuals
# DONE: allow display of type=c("observed", "expected")
# DONE: if ndim>3, provide for labels at max or min
# DONE: make object oriented and provide a loglm method

# mosaic3d: provide observed array of counts and either residuals, expected frequencies,
# or a loglin set of margins to fit

mosaic3d <- function(x, ...) {
	UseMethod("mosaic3d")
}

mosaic3d.loglm <-
function (x, type = c("observed", "expected"), 
    residuals_type = c("pearson", "deviance"), 
#    gp = shading_hcl, gp_args = list(), 
    ...) 
{
    residuals_type <- match.arg(tolower(residuals_type), c("pearson", "deviance"))
    if (is.null(x$fitted)) 
        x <- update(x, fitted = TRUE)
    expected <- fitted(x)
    residuals <- residuals(x, type = "pearson")
    observed <- residuals * sqrt(expected) + expected
    if (residuals_type == "deviance") 
        residuals <- residuals(x, type = "deviance")
#    gp <- if (inherits(gp, "grapcon_generator")) 
#        do.call("gp", c(list(observed, residuals, expected, x$df), 
#            as.list(gp_args)))
#    else gp
    mosaic3d.default(observed, residuals = residuals, expected = expected, 
        type = type, residuals_type = residuals_type, 
#        gp = gp, 
        ...)
}

mosaic3d.default <- function(x, expected=NULL, residuals=NULL, 
		type = c("observed", "expected"), residuals_type = NULL,
		shape=rgl::cube3d(alpha=alpha), alpha=0.5,
		spacing=0.1, split_dir=1:3, 
		shading=shading_basic, interpolate=c(2,4), zero_size=.05,
		label_edge,
		labeling_args=list(), newpage=TRUE, box=FALSE, ...) {
	
  if (!requireNamespace("rgl")) stop("rgl is required")

  type <- match.arg(type)
  if (is.null(residuals)) {
      residuals_type <- if (is.null(residuals_type)) 
          "pearson"
      else match.arg(tolower(residuals_type), c("pearson", 
          "deviance", "ft"))
  }

  ## convert structable object
  if (is.structable(x)) {
    x <- as.table(x)
  }

  ## table characteristics	
  levels <- dim(x)
  ndim <- length(levels)
  dn <- dimnames(x)
  if (is.null(dn))
    dn <- dimnames(x) <- lapply(levels, seq)
  vnames <- names(dimnames(x))
  if (is.null(vnames))
    vnames <- names(dn) <- names(dimnames(x)) <- LETTERS[1:ndim]

  ## replace NAs by 0
  if (any(nas <- is.na(x))) x[nas] <- 0

  ## model fitting:
  ## calculate expected if needed
  if ((is.null(expected) && is.null(residuals)) ||
      !is.numeric(expected)) {
      if (inherits(expected, "formula")) {
          fm <- loglm(expected, x, fitted = TRUE)
          expected <- fitted(fm)
          df <- fm$df
      } else {
          if (is.null(expected))
              expected <-  as.list(1:ndim)
          fm <- loglin(x, expected, fit = TRUE, print = FALSE)
          expected <- fm$fit
          df <- fm$df
      }
  }

  ## compute residuals
  if (is.null(residuals))
    residuals <- switch(residuals_type,
                        pearson = (x - expected) / sqrt(ifelse(expected > 0, expected, 1)),
                        deviance = {
                          tmp <- 2 * (x * log(ifelse(x == 0, 1, x / ifelse(expected > 0, expected, 1))) - (x - expected))
                          tmp <- sqrt(pmax(tmp, 0))
                          ifelse(x > expected, tmp, -tmp)
                        },
                        ft = sqrt(x) + sqrt(x + 1) - sqrt(4 * expected + 1)
                        )
  ## replace NAs by 0
  if (any(nas <- is.na(residuals))) residuals[nas] <- 0

	# switch observed and expected if required
  observed <- if (type == "observed") x else expected
  expected <- if (type == "observed") expected else x

	
	# replicate arguments to number of dimensions
	spacing <- rep(spacing, length=ndim)
	split_dir <- rep(split_dir, length=ndim)
	if(missing(label_edge)) 
		label_edge <- rep( c('-', '+'), each=3, length=ndim)
	
	zeros <- observed <= .Machine$double.eps
	
	shapelist <- shape
	# sanity check
	if (!inherits(shapelist, "shape3d")) stop("shape must be a shape3d object")
	
	if (newpage) rgl::open3d()
	for (k in 1:ndim) {
		marg <- margin.table(observed, k:1)
		if (k==1) {
			shapelist <- split3d(shapelist, marg, split_dir[k], space=spacing[k])
			label3d(shapelist, split_dir[k], dn[[k]], vnames[k], edge=label_edge[k], ...)
		}
		else {
			marg <- matrix(marg, nrow=levels[k])
			shapelist <- split3d(shapelist, marg, split_dir[k], space=spacing[k])
			names(shapelist) <- apply(as.matrix(expand.grid(dn[1:k])), 1, paste, collapse=":")
			L <- length(shapelist)
			label_cells <- if (label_edge[k]=='-') 1:levels[k] 
											else (L-levels[k]+1):L
			label3d(shapelist[label_cells], split_dir[k], dn[[k]], vnames[k], edge=label_edge[k], ...)
		}
	}

	# assign colors
	# TODO: allow alpha to control transparency of side walls
	col <- shading(residuals, interpolate=interpolate)

	# display, but exclude the zero cells
	rgl::shapelist3d(shapelist[!as.vector(zeros)], col=col[!as.vector(zeros)], ...)

	# plot markers for zero cells
	if (any(zeros)) {
		ctrs <- t(sapply(shapelist, center3d))
		rgl::spheres3d(ctrs[as.vector(zeros),], radius=zero_size)
	}

#	invisible(structable(observed))
	invisible(shapelist)
	
}


# basic shading_Friendly, adapting the simple code used in mosaicplot()

shading_basic <- function(residuals, interpolate=TRUE) {
	if (is.logical(interpolate)) 
		interpolate <- c(2, 4)
	else if (any(interpolate <= 0) || length(interpolate) > 5) 
		stop("invalid 'interpolate' specification")
	shade <- sort(interpolate)
	breaks <- c(-Inf, -rev(shade), 0, shade, Inf)
	colors <- c(hsv(0, s = seq.int(1, to = 0, length.out = length(shade) + 
									1)), hsv(4/6, s = seq.int(0, to = 1, length.out = length(shade) + 
									1)))
	colors[as.numeric(cut(residuals, breaks))]
}


# provide labels for 3D objects below/above their extent along a given dimension
# FIXME: kludge for interline gap between level labels and variable name
# TODO: how to pass & extract labeling_args, e.g., labeling_args=list(at='min', fontsize=10)

label3d <- function(objlist, dim, text, varname, offset=.05, adj, edge="-", gap=.1, labeling_args, ...) {
	if(missing(adj)) {
		if (dim < 3) adj <- ifelse(edge == '-', c(0.5, 1), c(0.5, 0))
		else         adj <- ifelse(edge == '-', c(1, 0.5), c(0, 0.5))
	}
	ranges <- lapply(objlist, range3d)
	loc <- t(sapply(ranges, colMeans))   # positions of labels on dimension dim
	min <- t(sapply(ranges, function(x) x[1,]))  # other dimensions at min values
	max <- t(sapply(ranges, function(x) x[2,]))  # other dimensions at max values
	xyz <- if (edge == '-') (min - offset) else (max + offset)
	xyz[,dim] <- loc[,dim]
	if(!missing(varname)) {
		loclab <- colMeans(loc)               # NB: doesn't take space into acct
		xyzlab <- if (edge == '-') 
			min[1,] - offset - gap
			else max[1,] + offset + gap
		xyzlab[dim] <- loclab[dim]
		xyz <- rbind(xyz, xyzlab)
		text <- c(text, varname)
	}
	result <- c(labels = rgl::texts3d(xyz, texts=text, adj=adj, ...))
	invisible(result)
}



