PlotHull <- function(dat, xlim = c(-180, 180), ylim = c(-90,90), 
                     col = rgb(255, 0, 0, 50, maxColorValue = 255), 
                     border = rgb(255, 0, 0, 50, maxColorValue = 255), 
                     type = "overlay", select = "all", ...) {
  if (sum(is(dat) == "list") > 0 & sum(is(dat[[1]]) == "SpatialPolygons") > 0) {
    if (select[1] == "all") {
      if (type == "overlay") {
        map("world", xlim = xlim, ylim = ylim)
        axis(1)
        axis(2)
        box("plot")
        for (i in 1:length(dat)) {
          plot(dat[[i]], add = T, col = col, border = border, ...)
        }
      }
      if (type == "separate") {
        for (i in 1:length(dat)) {
          map("world", xlim = xlim, ylim = ylim)
          axis(1)
          axis(2)
          box("plot")
          title(names(dat)[[i]])
          plot(dat[[i]], add = T, col = col, border = border, ...)
        }
      }
    }
    if (length(select) > 1 | (length(select) == 1 & select[1] != "all")) {
      if (type == "overlay") {
        map("world", xlim = xlim, ylim = ylim)
        axis(1)
        axis(2)
        box("plot")
        title(names(dat)[names(dat) %in% select])
        lapply(dat[names(dat) %in% select], function(x) plot(x, add = T, col = col, border = border, ...))
      }
      if (type == "separate") {
        for (i in 1:length(select)) {
          map("world", xlim = xlim, ylim = ylim)
          axis(1)
          axis(2)
          box("plot")
          title(select)
          lapply(dat[names(dat) %in% select], function(x) plot(x, add = T, col = col, border = border, ...))
        }
      }
    }
  } else {
    warning("'dat' must be a list of \"SpatialPolygons\"")
  }
} 