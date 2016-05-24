##' Combine \code{data.frame}s
##'
##' This function combines \code{data.frame}s by filling in missing
##' variables. This is useful for combining data from
##' sampled locations with prediction locations.
##'
##' If \code{fillwith} is a named object, its names must correspond to
##' the names of variables in the data frames. If a variable is
##' missing, then it is filled with the corresponding value in
##' \code{fillwith}. \code{fillwith} can contain only one unnamed
##' component which corresponds to the default filling.
##' @title Combine \code{data.frame}s
##' @param ... \code{data.frame}s or objects that can be coerced to
##' \code{data.frame}s
##' @param fillwith Which value to use for missing variables. This
##' could be a scalar, a named vector, or a named list with one value
##' in each component; see Details.
##' @param keepclass Whether to preserve the \code{\link[base]{class}}
##' of each variable. The elements in \code{fillwith} are coerced to
##' the corresponding variable's class.
##' @return A stacked \code{data.frame}.
##' @export 
##' @examples
##' \dontrun{
##' d1 <- data.frame(w = 1:3, z = 4:6 + 0.1)
##' d2 <- data.frame(w = 3:7, x = 1:5, y = 6:10)
##' (d12a <- stackdata(d1, d2))
##' lapply(d12a, class)
##' (d12b <- stackdata(d1, d2, fillwith = c(x = NA, y = 0, z = -99)))
##' lapply(d12b, class)
##' (d12c <- stackdata(d1, d2, fillwith = c(x = NA, y = 0, z = -99),
##'                    keepclass = TRUE))
##' lapply(d12c, class)
##' (d12d <- stackdata(d1, d2, fillwith = c(x = NA, 0)))
##'
##' data(rhizoctonia)
##' predgrid <- mkpredgrid2d(rhizoctonia[c("Xcoord", "Ycoord")],
##'                          par.x = 100, chull = TRUE, exf = 1.2)
##' rhizdata <- stackdata(rhizoctonia, predgrid$grid)
##' }
stackdata <- function (..., fillwith = NA, keepclass = FALSE) {
  fillNA <- function (d, allnames, fillwith) {
    miss <- allnames[!(allnames %in% names(d))]
    d[miss] <- fillwith[miss]
    d
  }
  input <- lapply(list(...), data.frame)
  nmall <- unique(unlist(lapply(input, names)))
  if (all(is.na(fillwith))) {
    fillwith <- rep(NA, length(nmall))
    names(fillwith) <- nmall
    keepclass <- FALSE
  } else {
    if (is.numeric(fillwith) | is.logical(fillwith)) {
      if (length(fillwith) == 1) {
        fillwith <- rep(fillwith, length(nmall))
        fillwith <- as.list(fillwith)
        names(fillwith) <- nmall
      } else {
        fillwith <- lapply(nmall, function(nm)
          do.call("switch", c(list(nm), as.list(fillwith))))
        names(fillwith) <- nmall
      }
    } else if (!is.list(fillwith)) {
      stop ("Argument fillwith must be either numeric or list")
    } else {
      fillwith <- lapply(nmall, function(nm)
        do.call("switch", c(list(nm), fillwith)))
      names(fillwith) <- nmall
    }      
  }
  fmiss <- nmall[!(nmall %in% names(fillwith))]
  fillwith[fmiss] <- NA
  fillwith <- fillwith[nmall]
  if(keepclass) {
    types <- lapply(stackdata(..., fillwith = NA), class)
    typef <- lapply(types, function (c) match.fun(paste("as.", c, sep = "")))
    fillwith <- lapply(nmall, function(nm) typef[[nm]](fillwith[[nm]]))
    names(fillwith) <- nmall
    # fillwith <- mapply("class<-", fillwith, types, SIMPLIFY = FALSE)
  }
  newdata <- lapply(input, fillNA, allnames = nmall, fillwith = fillwith)
  out <- do.call("rbind", newdata)
  out
}

##' This function creates a grid for prediction.
##'
##' If \code{chull} this function first calculates the convex hull of
##' the points. If \code{exf} is not 1 the borders are expanded. Then
##' the function calls \code{\link[sp]{point.in.polygon}} to select
##' points that fall inside the borders.
##' @title Make prediction grid
##' @param pnts.x x coordinate of the domain. Could also be a
##' two-column matrix containing the x and y coordinates
##' @param pnts.y y coordinate of the domain. Should be
##' omitted or set to \code{NULL} if the argument \code{pnts.x} is a
##' two-column matrix.
##' @param par.x A scalar parameter for the x component of the new
##' grid. This parameter corresponds to either the \code{by} or the
##' \code{length.out} arguments of the function
##' \code{\link[base]{seq}}. Could also be a vector of two elements
##' containing the parameter for x and y.
##' @param par.y As in \code{par.x} for the y component of the new
##' grid. Should be omitted or set to \code{NULL} if the argument
##' \code{par.x} is a two-dimensional vector.
##' @param isby If \code{TRUE}, the arguments \code{par.x} and
##' \code{par.y} correspond to the \code{by} argument of the function
##' \code{\link[base]{seq}}, otherwise they correspond to
##' \code{length.out}.
##' @param chull Whether to calculate the convex hull of the points.
##' Set this to \code{TRUE} if \code{pnts.x} and \code{pnts.y} denote
##' the sampled locations. If they correspond to the borders of the
##' domain, it is recommended to set this to \code{FALSE}.
##' @param exf An expansion factor of the convex hull of
##' \code{cbind(pnts.x, pnts.y)}. Must be positive. If larger or
##' smaller than 1, the convex hull is respectively expanded or
##' contracted.
##' @return A list with components
##' \itemize{
##' \item \code{grid} A two-column matrix with the prediction grid
##' \item \code{xycoord} A list with components "x" and "y"
##' containing the sequence of points used to create the grid
##' \item \code{xygrid} A matrix with the full square grid derived
##' from \code{xycoord}
##' \item \code{borders} The (expanded) borders of the domain
##' }
##' @seealso \code{\link[geoR]{pred_grid}}
##' @importFrom sp point.in.polygon
##' @export
##' @examples
##' \dontrun{
##' data(rhizoctonia)
##' predgrid <- mkpredgrid2d(rhizoctonia[c("Xcoord", "Ycoord")],
##'                          par.x = 100, chull = TRUE, exf = 1.2)
##' plot(predgrid$borders, type = "l")         # Domain for prediction
##' points(predgrid$grid, pch = 20, cex = .3)  # Prediction locations
##' points(rhizoctonia[c("Xcoord", "Ycoord")]) # Observed locations
##' }
mkpredgrid2d <- function (pnts.x, pnts.y, par.x, par.y, isby = FALSE,
                          chull = FALSE, exf = 1) {
  if (exf <= 0) stop ("Argument exf must be positive")
  if (missing(pnts.y)) pnts.y <- NULL
  if (missing(par.y)) par.y <- NULL
  if (!is.null(pnts.x)) pnts.x <- as.matrix(pnts.x)
  if (!is.null(pnts.y)) pnts.y <- as.matrix(pnts.y)
  ph <- cbind(pnts.x, pnts.y)
  nm <- colnames(ph)
  par <- c(par.x[1], par.y[1])
  d <- NCOL(ph)
  if (d != 2) stop ("Can only generate 2-dimensional grids")
  if (isTRUE(chull)) {
    ph <- ph[chull(ph), , drop = FALSE]
  }
  if (exf != 1) { ## Extend the covex hull by a factor exf
    centr <- colMeans(ph)     # Central coordinate of the polygon
    phu <- cbind(ph[, 1, drop = FALSE] - centr[1],
                 ph[, 2, drop = FALSE] - centr[2]) # Uncenter
    phu_r <- exf*sqrt(phu[, 2, drop = FALSE]^2 +
                        phu[, 1, drop = FALSE]^2) # Covert to polar
    phu_u <- atan2(phu[, 2, drop = FALSE], phu[, 1, drop = FALSE])
    phu_x <- phu_r*cos(phu_u) # Convert back to cartesian
    phu_y <- phu_r*sin(phu_u)
    ph[, 1] <- phu_x + centr[1]
    ph[, 2] <- phu_y + centr[2]
  }
  ft <- apply(ph, 2, range)
  if (isTRUE(!isby)) {
    par <- ceiling(par)
    par <- ((ft[2, , drop = FALSE] - ft[1, , drop = FALSE])/(par - 1))
  } else {
    par <- rep(par, length.out = d)
  }
  xycoord <- lapply(1:d, function (i) seq(ft[1, i], ft[2, i], par[i]))
  names(xycoord) <- nm
  eg <- as.matrix(expand.grid(xycoord, KEEP.OUT.ATTRS = FALSE))
  dimnames(eg) <- list(NULL, nm)
  iin <- sp::point.in.polygon(eg[, 1, drop = FALSE],
                              eg[, 2, drop = FALSE],
                              ph[, 1, drop = FALSE],
                              ph[, 2, drop = FALSE]) > 0
  grid <- eg[iin, , drop = FALSE]
  out <- list(grid = grid, xycoord = xycoord, xygrid  = eg, borders = ph)
  out
}
