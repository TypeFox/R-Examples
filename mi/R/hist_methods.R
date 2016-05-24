# Part of the mi package for multiple imputation of missing data
# Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015 Trustees of Columbia University
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


setMethod("hist", signature(x = "missing_variable"), def = 
  function(x, ...) {
    y <- x@data
    NAs <- is.na(x)
    
    h_all  <- hist(y, plot = FALSE)
    plot(h_all, border = "lightgray", main = "", xlab = if(x@done) "Completed" else "Observed", axes = FALSE,
         mgp = c(2, 1, 0), tcl = .05, col = if(x@done) "lightgray" else "blue", freq = TRUE, ...)
    axis(1, lwd = 0)
    axis(2)
    if(x@done) {
      h_obs  <- hist(y[!NAs], breaks = h_all$breaks, plot = FALSE)
      h_miss <- hist(y[NAs], breaks = h_all$breaks, plot = FALSE)
      
      segments(h_obs$breaks[1], 0,  y1 = h_obs$counts[1], col = "blue")
      segments(h_miss$breaks[1], 0,  y1 = h_miss$counts[1], col = "red")
      
      segments(h_obs$breaks[1], y0 = h_obs$counts[1], x1  = h_obs$breaks[2], col = "blue")
      segments(h_miss$breaks[1], y0 = h_miss$counts[1], x1 = h_miss$breaks[2], col = "red")
      
      for(i in 2:(length(h_obs$breaks)-1)) {
        segments(x0 = h_obs$breaks[i], y0 = h_obs$counts[i-1], y1 = h_obs$counts[i], col = "blue")
        segments(x0 = h_miss$breaks[i], y0 = h_miss$counts[i-1], y1 = h_miss$counts[i], col = "red")
        
        segments(x0 = h_obs$breaks[i], y0 = h_obs$counts[i], x1 = h_obs$breaks[i+1], col = "blue")
        segments(x0 = h_miss$breaks[i], y0 = h_miss$counts[i], x1 = h_miss$breaks[i+1], col = "red")
      }
      
      segments(x0 = h_obs$breaks[i+1], y0 = h_obs$counts[i], y1 = 0, col = "blue")
      segments(x0 = h_miss$breaks[i+1], y0 = h_miss$counts[i], y1 = 0, col = "blue")
      
      if(.MI_DEBUG) stopifnot(all(h_all$counts == (h_obs$counts + h_miss$counts)))
    }
    return(invisible(NULL))
  })

setMethod("hist", signature(x = "semi-continuous"), def = 
  function(x, ...) {
    con <- complete(x@indicator, 0L) == 0
    y <- x@data[con]
    NAs <- is.na(x)[con]
    
    h_all  <- hist(y, plot = FALSE)
    plot(h_all, freq = TRUE, border = "lightgray", main = "", xlab = if(x@done) "Completed" else "Observed", axes = FALSE,
         mgp = c(2, 1, 0), tcl = .05, col = if(x@done) "lightgray" else "blue", xlim = range(x@data, na.rm = TRUE), ...)
    axis(1, lwd = 0)
    axis(2)
    if(x@done) {
      h_obs  <- hist(y[!NAs], breaks = h_all$breaks, plot = FALSE)
      h_miss <- hist(y[NAs], breaks = h_all$breaks, plot = FALSE)
      
      segments(h_obs$breaks[1], 0,  y1 = h_obs$counts[1], col = "blue")
      segments(h_miss$breaks[1], 0,  y1 = h_miss$counts[1], col = "red")
      
      segments(h_obs$breaks[1], y0 = h_obs$counts[1], x1  = h_obs$breaks[2], col = "blue")
      segments(h_miss$breaks[1], y0 = h_miss$counts[1], x1 = h_miss$breaks[2], col = "red")
      
      for(i in 2:(length(h_obs$breaks)-1)) {
        segments(x0 = h_obs$breaks[i], y0 = h_obs$counts[i-1], y1 = h_obs$counts[i], col = "blue")
        segments(x0 = h_miss$breaks[i], y0 = h_miss$counts[i-1], y1 = h_miss$counts[i], col = "red")
        
        segments(x0 = h_obs$breaks[i], y0 = h_obs$counts[i], x1 = h_obs$breaks[i+1], col = "blue")
        segments(x0 = h_miss$breaks[i], y0 = h_miss$counts[i], x1 = h_miss$breaks[i+1], col = "red")
      }
      
      segments(x0 = h_obs$breaks[i+1], y0 = h_obs$counts[i], y1 = 0, col = "blue")
      segments(x0 = h_miss$breaks[i+1], y0 = h_miss$counts[i], y1 = 0, col = "blue")
      
      NAs <- is.na(x)[!con]
      tab <- table(x@data[!con], NAs)
      for(i in 1:NROW(tab)) {
        segments(x0 = as.numeric(rownames(tab)[i]), y0 = 0, y1 = sum(tab[i,]), col = "lightgray", lty = "dashed")
        segments(x0 = as.numeric(rownames(tab)[i]), y0 = 0, y1 = tab[i,1], col = "blue", lty = "dashed")
        if(ncol(tab) == 2) segments(x0 = as.numeric(rownames(tab)[i]), y0 = 0, y1 = tab[i,2], col = "red", lty = "dashed")
      }
      
      if(.MI_DEBUG) stopifnot(all(h_all$counts == (h_obs$counts + h_miss$counts)))
    }
    else {
      tab <- table(x@data[!con])
      for(i in 1:NCOL(tab)) segments(x0 = as.numeric(names(tab)[i]), y0 = 0, y1 = tab[i], col = "blue", lty = "dashed")
    }
    return(invisible(NULL))
  })

setMethod("hist", signature(x = "categorical"), def = 
  function(x, ...) {
    y <- x@data
    values <- sort(unique(y))
    breaks <- c(min(values) - 0.5, values + 0.5)
    
    values <- unique(y)
    values <- sort(values[!is.na(values)])
    breaks <- c(sapply(values, FUN = function(x) c(x - .25, x + .25)))
    NAs <- is.na(x)
    
    h_all  <- hist(y, breaks, plot = FALSE)
    #   h_all$counts[h_all$counts == 0] <- NA_integer_
    plot(h_all, border = "lightgray", axes = FALSE, main = "", xlab = if(x@done) "Completed" else "Observed", 
         mgp = c(2, 1, 0), tcl = .05, col = if(x@done) "lightgray" else "blue", freq = TRUE, ylim = range(h_all$counts, na.rm = TRUE), ...)
    axis(1, at = values, labels = levels(x@raw_data), lwd = 0)
    axis(2)
    
    if(x@done) {
      h_obs  <- hist(y[!NAs], breaks, plot = FALSE)
      h_miss <- hist(y[NAs], breaks, plot = FALSE)
      
      counts_obs <- h_obs$counts
      counts_obs <- counts_obs
      
      counts_miss <- h_miss$counts
      counts_miss <- counts_miss
      
      segments(breaks[1], 0,  y1 = counts_obs[1], col = "blue")
      segments(breaks[1], 0,  y1 = counts_miss[1], col = "red")
      
      if(counts_obs[1])  segments(breaks[1], y0 = counts_obs[1], x1  = breaks[2], col = "blue")
      if(counts_miss[1]) segments(breaks[1], y0 = counts_miss[1], x1 = breaks[2], col = "red")
      
      for(i in 2:(length(breaks)-1)) {
        segments(x0 = breaks[i], y0 = counts_obs[i-1], y1 = counts_obs[i], col = "blue")
        segments(x0 = breaks[i], y0 = counts_miss[i-1], y1 = counts_miss[i], col = "red")
        
        if(counts_obs[i]) segments(x0 = breaks[i], y0 = counts_obs[i], x1 = breaks[i+1], col = "blue")
        if(counts_miss[i]) segments(x0 = breaks[i], y0 = counts_miss[i], x1 = breaks[i+1], col = "red")
      }
      
      segments(x0 = breaks[i+1], y0 = counts_obs[i], y1 = 0, col = "blue")
      segments(x0 = breaks[i+1], y0 = counts_miss[i], y1 = 0, col = "red")
      
      if(.MI_DEBUG) stopifnot(all(h_all$counts == (h_obs$counts + h_miss$counts)))
    }
    
    return(invisible(NULL))
  })

setMethod("hist", signature(x = "binary"), def = 
  function(x, ...) {
    y <- x@data
    if(max(y, na.rm = TRUE) > 1) y <- y - 1L
    values <- 0:1
    breaks <- c(-.5, .5, 1.5)
    breaks <- c(-.25, .25, .75, 1.25)
    NAs <- is.na(x)
    
    h_all  <- hist(y, breaks, plot = FALSE)
    #   h_all$counts[h_all$counts == 0] <- NA_integer_
    plot(h_all, border = "lightgray", axes = FALSE, main = "", xlab = if(x@done) "Completed" else "Observed", 
         mgp = c(2, 1, 0), tcl = .05, col = if(x@done) "lightgray" else "blue", freq = TRUE, ylim = range(h_all$counts, na.rm = TRUE), ...)
    axis(1, at = values, lwd = 0)
    axis(2)
    
    if(x@done) {
      h_obs  <- hist(y[!NAs], breaks, plot = FALSE)
      h_miss <- hist(y[NAs], breaks, plot = FALSE)
      
      counts_obs <- h_obs$counts
      counts_obs <- counts_obs
      
      counts_miss <- h_miss$counts
      counts_miss <- counts_miss
      
      segments(breaks[1], 0,  y1 = counts_obs[1], col = "blue")
      segments(breaks[1], 0,  y1 = counts_miss[1], col = "red")
      
      if(counts_obs[1])  segments(breaks[1], y0 = counts_obs[1], x1  = breaks[2], col = "blue")
      if(counts_miss[1]) segments(breaks[1], y0 = counts_miss[1], x1 = breaks[2], col = "red")
      
      for(i in 2:(length(breaks)-1)) {
        segments(x0 = breaks[i], y0 = counts_obs[i-1], y1 = counts_obs[i], col = "blue")
        segments(x0 = breaks[i], y0 = counts_miss[i-1], y1 = counts_miss[i], col = "red")
        
        if(counts_obs[i])  segments(x0 = breaks[i], y0 = counts_obs[i], x1 = breaks[i+1], col = "blue")
        if(counts_miss[i]) segments(x0 = breaks[i], y0 = counts_miss[i], x1 = breaks[i+1], col = "red")
      }
      
      segments(x0 = breaks[i+1], y0 = counts_obs[i], y1 = 0, col = "blue")
      segments(x0 = breaks[i+1], y0 = counts_miss[i], y1 = 0, col = "red")
      
      if(.MI_DEBUG) stopifnot(all(h_all$counts == (h_obs$counts + h_miss$counts)))
    }
    
    return(invisible(NULL))
  })

setMethod("hist", signature(x = "missing_data.frame"), def = 
  function(x, ask = TRUE, ...) {
    k <- sum(!x@no_missing)
    if (.Device != "null device" && x@done) {
      oldask <- grDevices::devAskNewPage(ask = ask)
      if (!oldask) on.exit(grDevices::devAskNewPage(oldask), add = TRUE)
      op <- options(device.ask.default = TRUE)
      on.exit(options(op), add = TRUE)
    }
    par(mfrow = n2mfrow(k))
    for(i in 1:x@DIM[2]) {
      if(x@no_missing[i]) next
      hist(x@variables[[i]])
      header <- x@variables[[i]]@variable_name
      if(is(x@variables[[i]], "continuous")) {
        trans <- .show_helper(x@variables[[i]])$transformation[1]
        header <- paste("\n", header, " (", trans, ")", sep = "")
      }
      title(main = header)
    }
    return(invisible(NULL))
  })

setMethod("hist", signature(x = "mdf_list"), def = 
  function(x, ask = TRUE, ...) {
    if (.Device != "null device") {
      oldask <- grDevices::devAskNewPage(ask = ask)
      if (!oldask) on.exit(grDevices::devAskNewPage(oldask), add = TRUE)
      op <- options(device.ask.default = ask)
      on.exit(options(op), add = TRUE)
    }
    sapply(x, FUN = hist, ...)
    return(invisible(NULL))
  })

setMethod("hist", signature(x = "mi"), def = 
  function(x, m = 1:length(x), ask = TRUE, ...) {
    for(i in m) hist(x@data[[i]], ask = ask, ...)
    return(invisible(NULL))
  })

setMethod("hist", signature(x = "mi_list"), def = 
  function(x, m = 1:length(x), ask = TRUE, ...) {
    if (.Device != "null device") {
      oldask <- grDevices::devAskNewPage(ask = ask)
      if (!oldask) on.exit(grDevices::devAskNewPage(oldask), add = TRUE)
      op <- options(device.ask.default = ask)
      on.exit(options(op), add = TRUE)
    }
    sapply(x, FUN = hist, m = m, ask = ask, ...)
    return(invisible(NULL))
  })
