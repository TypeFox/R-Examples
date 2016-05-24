# Part of the mi package for multiple imputation of missing data
# Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015 Trustees of Columbia University
# Copyright (C) 2011 Douglas Bates and Martin Maechler
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

setMethod("image", "dgTMatrix", # slight hack of a method in the library(Matrix)
          function(x,
                   xlim = .5 + c(0, di[2]),
                   ylim = .5 + c(di[1], 0),
                   aspect = "iso", ## was default "fill"
                   sub = sprintf("Dimensions: %d x %d", di[1], di[2]),
                   xlab = "Column", ylab = "Row", cuts = 15,
                   useAbs = NULL, colorkey = !useAbs, col.regions = NULL,
                   lwd = NULL, ...)
          {
            di <- x@Dim
            xx <- x@x
            if(missing(useAbs)) ## use abs() when all values are non-neg
              useAbs <- min(xx, na.rm=TRUE) >= 0
            else if(useAbs)
              xx <- abs(xx)
            rx <- range(xx, finite=TRUE)
            if(is.null(col.regions))
              col.regions <-
                if(useAbs) {
                  grey(seq(from = 0.7, to = 0, length = 100))
                } else { ## no abs(.), rx[1] < 0
                  nn <- 100
                  n0 <- min(nn, max(0, round((0 - rx[1])/(rx[2]-rx[1]) * nn)))
                  col.regions <-
                    c(colorRampPalette(c("blue3", "gray80"))(n0),
                      colorRampPalette(c("gray75","red3"))(nn - n0))
                }
            if(!is.null(lwd) && !(is.numeric(lwd) && all(lwd >= 0))) # allow lwd=0
              stop("'lwd' must be NULL or non-negative numeric")
            lattice::levelplot(x@x ~ (x@j + 1L) * (x@i + 1L),
                      sub = sub, xlab = xlab, ylab = ylab,
                      xlim = xlim, ylim = ylim, aspect = aspect,
                      colorkey = colorkey, col.regions = col.regions, cuts = cuts,
                      # 		    par.settings = list(background = list(col = "transparent")),
                      panel = function(x, y, z, subscripts, at, ..., col.regions)
                      {
                        x <- as.numeric(x[subscripts])
                        y <- as.numeric(y[subscripts])
                        numcol <- length(at) - 1
                        num.r <- length(col.regions)
                        col.regions <-
                          if (num.r <= numcol)
                            rep(col.regions, length = numcol)
                        else col.regions[1+ ((1:numcol-1)*(num.r-1)) %/% (numcol-1)]
                        zcol <- rep.int(NA_integer_, length(z))
                        for (i in seq_along(col.regions))
                          zcol[!is.na(x) & !is.na(y) & !is.na(z) &
                            at[i] <= z & z < at[i+1]] <- i
                        zcol <- zcol[subscripts]
                        if (any(subscripts)) {
                          if(is.null(lwd)) {
                            wh <- grid::current.viewport()[c("width", "height")]
                            ## wh : current viewport dimension in pixel
                            wh <- c(grid::convertWidth(wh$width, "inches",
                                                       valueOnly=TRUE),
                                    grid::convertHeight(wh$height, "inches",
                                                        valueOnly=TRUE)) *
                                                          par("cra") / par("cin")
                            pSize <- wh/di ## size of one matrix-entry in pixels
                            pA <- prod(pSize) # the "area"
                            p1 <- min(pSize)
                            lwd <- ## crude for now
                              if(p1 < 2 || pA < 6) 0.01 # effectively 0
                            else if(p1 >= 4) 1
                            else if(p1 > 3) 0.5 else 0.2
                          } else stopifnot(is.numeric(lwd), all(lwd >= 0)) # allow 0
                          
                          grid::grid.rect(x = x, y = y, width = 1, height = 1,
                                           default.units = "native",
                                           gp = grid::gpar(fill = ifelse(is.na(zcol), "black", col.regions[zcol]),
                                                            lwd = lwd, col = if(lwd < .01) NA else NA))
                        }
                      }, ...)
          })

setMethod("image", signature(x = "missing_data.frame"), def =
  function (x, y.order = FALSE, x.order = FALSE,  clustered = TRUE, grayscale = FALSE, ...) {
    data <- lapply(x@variables, FUN = function(z) if(is(z, "irrelevant")) NULL else is.na(z) * 1)
    data <- as.matrix(as.data.frame(data[!sapply(data, is.null)]))
    index <- seq(nrow(data))
    x.at <- 1:nrow( data )
    x.lab <- index
    if( x.order ) {
      orderIndex <- order(colSums(data), decreasing = TRUE)
      sub <- "Ordered by number of missing items per variable" 
    }
    if( y.order ) {
      orderIndex <- order(rowSums(data), decreasing = FALSE)
      index <- row.names( data )
      sub <- "Ordered by number of missing items per observation" 
      x.at <- NULL
      x.lab <- FALSE
    }
    if(clustered){
      orderIndex <-  order.dendrogram(as.dendrogram(hclust(dist(data, method = "binary"), method="mcquitty")))
      sub <- "Clustered by missingness"
    }
    if(!grayscale) {
      data <- lapply(x@variables, FUN = function(z) {
        y <- z@data
        if(is(z, "irrelevant")) return(NULL)
        else if(is(z, "continuous"))  return( (y - mean(y, na.rm = TRUE)) / (2 * sd(y, na.rm = TRUE)) )
        else if(is(z, "count"))  return( (y - mean(y, na.rm = TRUE)) / (2 * sd(y, na.rm = TRUE)) )
        else if(is(z, "categorical")) {
          y <- z@data
          if(is(z, "binary")) {
            y <- y == max(y, na.rm = TRUE)
            return( (y - 0.5) * 2 )
          }
          else {
            the_range <- seq(from = -.99, to = 1, length.out = length(unique(na.omit(y))))
            return(the_range[as.integer(as.factor(y))])
          }
        }
        else return( (y - mean(y, na.rm = TRUE)) / (2 * sd(y, na.rm = TRUE)) ) 
      })
      data <- as.matrix(as.data.frame(data[!sapply(data, is.null)]))      
    }
    
    if(y.order) X <- Matrix(data[,orderIndex])
    else        X <- Matrix(data[orderIndex,])
    
    if(grayscale) {
      plot(image(X, aspect = "fill", xlab = "Standardized Variable", ylab = "Observation Number", sub = sub,
                 scales = list(x = list(at = 1:ncol(data), labels = colnames(data), rot = 90, abbreviate = TRUE, minlength = 8)), 
                 main = "Dark represents missing data", colorkey = FALSE, alpha.regions = 1, ...))
      return(invisible(NULL))
    }
    nn <- 100
    rx <- range(X, finite = TRUE)
    n0 <- min(nn, max(0, round((0 - rx[1])/(rx[2]-rx[1]) * nn)))
    col.regions <- heat.colors(17)
    breaks <- seq(from = rx[1] - 1e-8, to = rx[2] + 1e-8, length.out = 16)
    plot(image(X, aspect = "fill", xlab = "Standardized Variable", ylab = "Observation Number", sub = sub, at = breaks,
               scales = list(x = list(at = 1:ncol(data), labels = colnames(data), rot = 90, abbreviate = TRUE, minlength = 8)), 
               main = "Dark represents missing data", colorkey = TRUE, col.regions = col.regions, alpha.regions = 1, ...))
    
    return(invisible(NULL))
  })

setMethod("image", signature(x = "mdf_list"), def =
  function (x, y.order = FALSE, x.order = FALSE,  clustered = TRUE, grayscale = FALSE, ask = TRUE, ...) {
    if (.Device != "null device") {
      oldask <- grDevices::devAskNewPage(ask = ask)
      if (!oldask) on.exit(grDevices::devAskNewPage(oldask), add = TRUE)
      op <- options(device.ask.default = ask)
      on.exit(options(op), add = TRUE)
    }
    sapply(x, FUN = image, y.order = y.order, x.order = x.order, clustered = clustered, grayscale = grayscale, ...)
    return(invisible(NULL))
  })

setMethod("image", signature(x = "mi"), def =
  function (x, y.order = FALSE, x.order = FALSE,  clustered = TRUE, ...) {
    data <- lapply(x@data[[1]]@variables, FUN = function(z) if(is(z, "irrelevant")) NULL else is.na(z) * 1)
    data <- as.matrix(as.data.frame(data[!sapply(data, is.null)]))
    if( x.order ) {
      orderIndex <- order(colSums(data), decreasing = TRUE)
      sub <- "Ordered by number of missing items per variable" 
    }
    if( y.order ) {
      orderIndex <- order(rowSums(data), decreasing = FALSE)
      index <- row.names( data )
      sub <- "Ordered by number of missing items per observation" 
    }
    if(clustered){
      orderIndex <-  order.dendrogram(as.dendrogram(hclust(dist(data, method = "binary"), method="mcquitty")))
      sub <- "Clustered by missingness"
    }
    foo <- function(z, raw = FALSE) {
      y <- if(raw) z@raw_data else z@data
      # 	    y <- z@data
      if(is(z, "irrelevant")) return(NULL)
      else if(is(z, "continuous"))  return( (y - mean(y, na.rm = TRUE)) / (2 * sd(y, na.rm = TRUE)) )
      else if(is(z, "count"))  return( (y - mean(y, na.rm = TRUE)) / (2 * sd(y, na.rm = TRUE)) )
      else if(is(z, "categorical")) {
        y <- if(raw) as.numeric(z@raw_data) else z@data
        if(is(z, "binary")) {
          y <- y == max(y, na.rm = TRUE)
          return( (y - 0.5) * 2 )
        }
        else {
          the_range <- seq(from = -.99, to = 1, length.out = length(unique(na.omit(y))))
          return(the_range[as.integer(as.factor(y))])
        }
      }
      else return( (y - mean(y, na.rm = TRUE)) / (2 * sd(y, na.rm = TRUE)) ) 
    }
    temp <- lapply(x@data[[1]]@variables, FUN = foo)
    temp <- as.matrix(as.data.frame(temp[!sapply(temp, is.null)]))
    temp[data == 1] <- NA_real_
    data <- temp
    if(y.order) data <- data[,orderIndex]
    else        data <- data[orderIndex,]
    X0 <- Matrix(data)
    data <- 0
    chains <- min(3, length(x@data))
    for(i in seq_along(x@data)) {
      temp <- lapply(x@data[[i]]@variables, FUN = foo, raw = FALSE)
      temp <- as.matrix(as.data.frame(temp[!sapply(temp, is.null)]))
      data <- data +  temp / chains 
    }
    if(y.order) data <- data[,orderIndex]
    else        data <- data[orderIndex,]
    X1 <- Matrix(data)
    X <- rbind2(X0, X1)
    breaks <- seq(from = min(X, na.rm = TRUE), to = max(X, na.rm = TRUE), length.out = 15)
    plot(image(X0, aspect = "fill", xlab = "", ylab = "Observation Number", sub = "", at = breaks,
               scales = list(x = list(at = 1:ncol(data), labels = colnames(data), rot = 90, abbreviate = TRUE, minlength = 5)), 
               main = "Original data", colorkey = TRUE, col.regions = heat.colors(17), ...), 
         split = c(1,1,1,2))
    
    plot(image(X1, aspect = "fill", xlab = "", ylab = "Observation Number", sub = "", at = breaks,
               scales = list(x = list(at = 1:ncol(data), labels = colnames(data), rot = 90, abbreviate = TRUE, minlength = 5)), 
               main = "Average completed data", colorkey = TRUE, col.regions = heat.colors(17), ...), 
         newpage = FALSE, split = c(1,2,1,2))
    return(invisible(NULL))
  })

setMethod("image", signature(x = "mi_list"), def =
  function (x, y.order = FALSE, x.order = FALSE,  clustered = TRUE, grayscale = FALSE, ask = TRUE, ...) {
    if (.Device != "null device") {
      oldask <- grDevices::devAskNewPage(ask = ask)
      if (!oldask) on.exit(grDevices::devAskNewPage(oldask), add = TRUE)
      op <- options(device.ask.default = ask)
      on.exit(options(op), add = TRUE)
    }
    sapply(x, FUN = image, y.order = y.order, x.order = x.order, clustered = clustered, grayscale = grayscale, ...)
    return(invisible(NULL))
  })

.binnedplot <-
  function (x, y, nclass = NULL, xlab = "Expected Values", ylab = "Average residual", 
            main = "", cex.pts = 0.8, col.pts = "blue", col.int = "gray") 
  {
    n <- length(x)
    if (is.null(nclass)) {
      if (n >= 100) {
        nclass = floor(sqrt(length(x)))
      }
      if (n > 10 & n < 100) {
        nclass = 10
      }
      if (n <= 10) {
        nclass = floor(n/2)
      }
    }
    aa <- data.frame(arm::binned.resids(x, y, nclass)$binned)
    #     aa <- aa[!is.na(aa$X2se),] ## FIXME: remove once Yu-Sung fixes arm::binned.resids
    plot(range(aa$xbar), range(aa$ybar, aa$X2se, -aa$X2se), xlab = xlab, 
         ylab = ylab, type = "n", main = main, mgp = c(2, 1, 0), tcl = .05)
    abline(0, 0, lty = 2)
    lines(aa$xbar, aa$X2se, col = col.int)
    lines(aa$xbar, -aa$X2se, col = col.int)
    points(aa$xbar, aa$ybar, pch = 19, cex = cex.pts, col = col.pts)
  }

.binnedpoints <-
  function (x, y, nclass = NULL, cex.pts = 0.8, col.pts = "red") 
  {
    n <- length(x)
    if (is.null(nclass)) {
      if (n >= 100) {
        nclass = floor(sqrt(length(x)))
      }
      if (n > 10 & n < 100) {
        nclass = 10
      }
      if (n <= 10) {
        nclass = floor(n/2)
      }
    }
    if(n > 5) {
      aa <- data.frame(arm::binned.resids(x, y, nclass)$binned)
      points(aa$xbar, aa$ybar, pch = 19, cex = cex.pts, col = col.pts)
    }
    return(invisible(NULL))
  }

setMethod("plot", signature(x = "missing_data.frame", y = "missing_variable"), def =
  function(x, y, ...) {
    NAs <- is.na(y@raw_data)
    z <- y@data 
    hist(y)
    yhat <- y@fitted
    the_range <- range(c(yhat, z))
    plot(the_range, the_range, type = "n", xlab = "Expected Values", ylab = "Completed",
         mgp = c(2, 1, 0), tcl = .05)
    abline(0, 1, lty = 2, col = "lightgray")
    points(yhat, z, col = ifelse(NAs, "red", "blue"), pch = ".", cex = 2)
    lines(lowess(x = yhat[!NAs], y = z[!NAs]), col = "blue")
    .binnedplot(yhat[!NAs], (y@data - yhat)[!NAs])
    .binnedpoints(yhat[NAs], (y@data - yhat)[NAs])
    return(invisible(NULL))
  })

setMethod("plot", signature(x = "missing_data.frame", y = "categorical"), def =
  function(x, y, ...) {
    NAs <- is.na(y@raw_data)
    z <- y@data
    hist(y)
    s <- nrow(y@parameters) + 1
    yhat <- y@fitted
    if(length(yhat) == 0) { # embedded
      varname <- y@variable_name
      varname <- strsplit(varname, ":")[[1]][1]
      to_drop <- x@index[[varname]]
      X <- x@X[,-to_drop]
      s <- nrow(y@parameters) + 1
      model <- fit_model(y, data = x, s = s, warn = TRUE, X = X)
      yhat <- fitted(model)
    }
    if(is.matrix(yhat)) yhat <- yhat %*% (1:ncol(yhat))
    the_range <- range(c(yhat, z)) #+ c(-.1, .1)
    plot(the_range, the_range, type = "n", xlab = "Expected Values", ylab = "Completed (jittered)",
         mgp = c(2, 1, 0), tcl = .05)
    abline(0, 1, lty = 2, col = "lightgray")
    points(yhat, jitter(z), col = ifelse(NAs, "red", "blue"), pch = ".", cex = 2)
    .binnedplot(yhat[!NAs], (y@data - yhat)[!NAs])
    .binnedpoints(yhat[NAs], (y@data - yhat)[NAs])
    return(invisible(NULL))
  })

setMethod("plot", signature(x = "missing_data.frame", y = "binary"), def =
  function(x, y, ...) {
    
    NAs <- is.na(y@raw_data)
    z <- y@data - 1L
    hist(y)
    s <- nrow(y@parameters) + 1
    yhat <- y@fitted
    if(length(yhat) == 0) { # embedded
      varname <- y@variable_name
      varname <- strsplit(varname, ":")[[1]][1]
      to_drop <- x@index[[varname]]
      X <- x@X[,-to_drop]
      model <- fit_model(y = y, data = x, s = s, warn = TRUE, X = X)
      yhat <- fitted(model)
    }
    the_range <- range(c(yhat, z)) #+ c(-.1, .1)
    plot(the_range, the_range, type = "n", xlab = "Expected Values", ylab = "Completed (jittered)",
         mgp = c(2, 1, 0), tcl = .05)
    abline(0, 1, lty = 2, col = "lightgray")
    points(yhat, jitter(z), col = ifelse(NAs, "red", "blue"), pch = ".", cex = 2)
    .binnedplot(yhat[!NAs], (y@data - 1 - yhat)[!NAs])
    .binnedpoints(yhat[NAs], (y@data - 1 - yhat)[NAs])
    return(invisible(NULL))
  })

setMethod("plot", signature(x = "allcategorical_missing_data.frame", 
                            y = "categorical"), def =
            function(x, y, ...) {
  NAs <- is.na(y@raw_data)
  z <- y@raw_data
  hist(y)
  latents <- x@latents@data
  yhat <- t(sapply(latents[!NAs], FUN = function(l) y@fitted[l,]))
  tab_obs <- table(z[!NAs])
  tab_model <- table(apply(yhat, 1, FUN = function(p) which(rmultinom(1, 1, p) == 1)))  
  the_range <- c(0, max(c(tab_obs, tab_model)))
  barplot(tab_obs, beside = TRUE, xlab = "Observed Values", ylim = the_range)
  names(tab_model) <- levels(z)
  barplot(tab_model, beside = TRUE, xlab = "Expected Values", ylim = the_range)
  return(invisible(NULL))
})

setMethod("plot", signature(x = "allcategorical_missing_data.frame", 
                            y = "binary"), def =
            function(x, y, ...) {
  NAs <- is.na(y@raw_data)
  z <- y@raw_data
  hist(y)
  latents <- x@latents@data
  yhat <- t(sapply(latents[!NAs], FUN = function(l) y@fitted[l,]))
  tab_obs <- table(z[!NAs])
  tab_model <- table(apply(yhat, 1, FUN = function(p) which(rmultinom(1, 1, p) == 1)))  
  the_range <- c(0, max(c(tab_obs, tab_model)))
  barplot(tab_obs, beside = TRUE, xlab = "Observed Values", ylim = the_range)
  names(tab_model) <- levels(z)
  barplot(tab_model, beside = TRUE, xlab = "Expected Values", ylim = the_range)
  return(invisible(NULL))
})

setMethod("plot", signature(x = "missing_data.frame", y = "semi-continuous"), def =
  function(x, y, ...) {
    NAs <- is.na(y@raw_data)
    z <- y@data 
    hist(y)
    s <- nrow(y@parameters) + 1
    yhat <- z
    yhat[complete(y@indicator, m = 0L, to_factor = TRUE) == 0] <- y@fitted #fitted(model)
    the_range <- range(c(yhat, z))
    plot(the_range, the_range, type = "n", xlab = "Expected Values", ylab = "Completed",
         mgp = c(2, 1, 0), tcl = .05)
    abline(0, 1, lty = 2, col = "lightgray")
    points(yhat, z, col = ifelse(NAs, "red", "blue"), pch = ".", cex = 2)
    lines(lowess(x = yhat[!NAs], y = z[!NAs]), col = "blue")
    .binnedplot(yhat[!NAs], (y@data - yhat)[!NAs])
    .binnedpoints(yhat[NAs], (y@data - yhat)[NAs])
    return(invisible(NULL))
  })

setMethod("plot", signature(x = "mi", y = "ANY"), def = 
  function(x, y, ask = TRUE, header = character(0), ...) {
    if(missing(y))           select <- 1:ncol(x@data[[1]])
    else if(is.logical(y))   select <- which(y)
    else if(is.character(y)) select <- which(colnames(x@data[[1]]) %in% y)
    else if(is.numeric(y))   select <- which(1:nrow(x@data[[1]]) %in% y)
    
    for(i in seq_along(x@data[[1]]@variables)) {
      if(x@data[[1]]@no_missing[i]) next
      else if(is(x@data[[1]]@variables[[i]], "irrelevant")) next
      else if(x@data[[1]]@variables[[i]]@imputation_method == "mcar") {
        warning(x@data[[1]]@variables[[i]]@variable_name, " not plotted because it assumes MCAR")
        next
      }
      if(!(i %in% select)) next
      l <- min(3, length(x@data))
      
      if (.Device != "null device") {
        oldask <- grDevices::devAskNewPage(ask = ask)
        if (!oldask) on.exit(grDevices::devAskNewPage(oldask), add = TRUE)
        op <- options(device.ask.default = ask)
        on.exit(options(op), add = TRUE)
      }
      
      par(mfrow = c(l,3), mar = c(5,4,1,1) + .1)
      if(is(x@data[[1]]@variables[[i]], "semi-continuous")) {
        for(j in 1:l) plot(x@data[[j]], x@data[[j]]@variables[[i]]@indicator, ...)
        title(main = paste("\n", header, x@data[[1]]@variables[[i]]@indicator@variable_name, sep = ""), outer = TRUE)
      }
      for(j in 1:l) plot(x@data[[j]], x@data[[j]]@variables[[i]], ...)      
      new_header <- paste(header, x@data[[1]]@variables[[i]]@variable_name)
      if(is(x@data[[1]]@variables[[i]], "continuous")) {
        trans <- .show_helper(x@data[[1]]@variables[[i]])$transformation[1]
        new_header <- paste("\n", new_header, " (", trans, ")", sep = "")
      }
      else new_header <- paste("\n", new_header, sep = "")
      title(main = new_header, outer = TRUE)
    }
    return(invisible(NULL))
  })

setMethod("plot", signature(x = "mi_list", y = "ANY"), def = 
  function(x, y, ask = TRUE, ...) {
    if(missing(y)) for(i in seq_along(x)) plot(x[[i]], ask = ask, header = paste(names(x)[i], ": ", sep = ""), ...)
    else for(i in seq_along(x)) plot(x[[i]], y = y, ask = ask, header = paste(names(x)[i], ": ", sep = ""), ...)
    return(invisible(NULL))
  })
