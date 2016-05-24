plot.ramps <- function(x, type = c("i", "c", "w"), col = tim.colors(64),
   func = mean, sites = FALSE, database = NULL, regions = ".",
   resolution = c(64, 64), bw = 1, ...)
{
   if (ncol(x$z) == 0) stop("No monitored latent spatial parameters to plot")
   if (length(resolution) < 2) stop("'resolution' must be a two-element",
                                    " numerical vector of width x height")

   coords <- x$coords[x$control$z$monitor,,drop = FALSE]
   if (sites) {
      val <- as.matrix(x$kmat) == 1
      idx <- unique(col(val)[val])
      sites <- x$coords[idx, c(1,2), drop = FALSE]
   } else {
      sites <- matrix(0, 0, 2)
   }
   plot.default(x$z, coords, type=match.arg(type), col=col, func=func,
                sites=sites, database=database, regions=regions,
                resolution=resolution, bw=bw, ...)
}


plot.predict.ramps <- function(x, type = c("i", "c", "w"), col = tim.colors(64),
   func = mean, database = NULL, regions = ".", resolution = c(64, 64), bw = 1,
   ...)
{
   if (length(resolution) < 2) stop("'resolution' must be a two-element",
                                    " numerical vector of width x height")

   plot.default(x, attr(x, "coords"), type=match.arg(type), col=col, func=func,
                sites=matrix(0, 0, 2), database=database, regions=regions,
                resolution=resolution, bw=bw, ...)
}


plot.default <- function(z, coords, type, col, func, sites, database, regions,
   resolution, bw, ...)
{
   mu <- apply(z, 2, func)

   if (!is.vector(mu))
      stop("'func' must return a summary statistic of the MCMC sampler output")

   val <- as.image(mu, x = coords, nrow = resolution[2], ncol = resolution[1])
   fit <- image.smooth(val$z, dx = val$x[2] - val$x[1], dy = val$y[2] - val$y[1],
             theta = bw * sqrt(diff(range(val$x))^2 + diff(range(val$y))^2) / 100)

   args <- list(x = val$x, y = val$y, z = fit$z, col = col, ...)
   val <- colnames(coords)
   if (is.null(args$xlab)) args$xlab <- val[1]
   if (is.null(args$ylab)) args$ylab <- val[2]

   switch(type,
      c = foo <- function(...) {
             args <- list(...)

             args1 <- args
             args1$nlevels <- NULL
             args1$levels <- NULL
             args1$labels <- NULL
             args1$labcex <- NULL
             args1$drawlabels <- NULL
             args1$method <- NULL
             args1$vfont <- NULL
             do.call("image", args1)

             args2 <- args["z"]
             args2$x <- args$x
             args2$y <- args$y
             args1$nlevels <- args$nlevels
             args1$levels <- args$levels
             args1$labels <- args$labels
             args2$zlim <- args$zlim
             args2$xlim <- args$xlim
             args2$ylim <- args$ylim
             args2$labcex <- args$labcex
             args2$drawlabels <- args$drawlabels
             args2$method <- args$method
             args2$vfont <- args$vfont
             args2$axes <- FALSE
             args2$frame.plot <- FALSE
             args2$lty <- args$lty
             args2$lwd <- args$lwd
             args2$add <- TRUE
             do.call("contour", args2)
          },
      i = foo <- image.plot,
      w = {
         foo <- drape.plot
         if (is.null(args$zlab)) args$zlab <- "z"
      }
   )

   if (is.null(database)) {
      do.call("foo", args)
   } else {
      outline <- map(database, regions, plot=F)

      val <- as.matrix(expand.grid(args$x, args$y))
      n <- nrow(val)
      val <- rbind(val, sites)
      val <- strsplit(map.where(database, val[,1], val[,2]), ",")
      idx <- !(unlist(lapply(val, function(x) x[[1]])) %in% regions)

      if (all(idx)) {
         stop("No spatial coordinates to map within the specified region")
      } else {
         args$z[idx[1:n]] <- NA
         sites[idx[-(1:n)],] <- NA
      }

      val <- range(outline$x, na.rm = T)
      if (is.null(args$xlim))
         args$xlim <- val + c(-0.025, 0.025) * abs(diff(val))
      val <- range(outline$y, na.rm = T)
      if (is.null(args$ylim))
         args$ylim <- val + c(-0.025, 0.025) * abs(diff(val)) 

      val <- do.call("foo", args)
      if (type == "w") {
         outline <- trans3d(outline$x, outline$y, min(args$z, na.rm=TRUE), val)
         if (nrow(sites) > 0)
            sites <- trans3d(sites[,1], sites[,2], min(args$z, na.rm=TRUE), val)
      }
      lines(outline)
   }
   points(sites)
}
