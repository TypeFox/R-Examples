# A function used to set or query data and parameters.

Data <- local({

  # Store data locally
  dat <- list()

  # Set default values
  default <- list("nlevels"       = 20,
                  "width"         = 7,
                  "cex.pts"       = 1,
                  "default.dir"   = getwd(),
                  "sep"           = "\t",
                  "rkey"          = 0,
                  "show.poly"     = 0,
                  "img.contour"   = 0,
                  "show.lines"    = 0,
                  "show.points"   = 0,
                  "vuni"          = 0,
                  "show.2.axes"   = 0,
                  "minor.ticks"   = 0,
                  "ticks.inside"  = 0,
                  "rm.pnt.line"   = 0,
                  "grid.res"      = list(x=NA, y=NA),
                  "grid.mba"      = list(n=NA, m=NA, h=11),
                  "color.palette" = grDevices::terrain.colors
              )

  ## Main program

  function(option, value, which.attr=NULL, clear.proj=FALSE, clear.data=FALSE,
           replace.all=NULL) {

    # Replace all values
    if (is.list(replace.all)) {
      dat <<- replace.all
      return(invisible())
    }

    # Save parameters
    if (clear.proj | clear.data) {
      save.params <- c("default.dir", "win.loc", "csi", "width", "cex.pts")
      if (clear.data)
        save.params <- c(save.params, "nlevels", "asp.yx", "asp.zx",
                         "vmax", "vxby", "vyby", "rkey", "show.poly",
                         "img.contour", "show.lines", "show.points",
                         "vuni", "date.fmt", "polys", "proj.file",
                         "show.2.axes", "minor.ticks", "ticks.inside",
                         "color.palette", "rm.pnt.line")
      save.params <- save.params[save.params %in% names(dat)]
      dat <<- sapply(save.params, function(i) list(dat[[i]]))
      return(invisible())
    }

    # Return all data
    if (missing(option))
      return(dat)

    # Check indices for numeric option elements
    if (is.numeric(option)) {
      option <- sapply(option, as.integer)
      option.new <- option[1L]
      if (option.new > length(dat))
        option.new <- NULL
      if (!is.null(option.new) && length(option) > 1) {
        for (i in 2:length(option)) {
          if (option[i] > length(dat[[option.new[-i]]]))
            break
          else
            option.new <- c(option.new, option[i])
        }
      }

    # Determine numeric indices from character option element
    } else {
      idx <- match(option[1], names(dat))
      option.new <- idx
      if (is.na(option.new))
        option.new <- NULL
      if (!is.null(option.new) && length(option) > 1) {
        for (i in 2:length(option)) {
          idx <- match(option[i], names(dat[[option.new[-i]]]))
          if (is.na(idx))
            break
          option.new <- c(option.new, idx)
        }
      }
    }

    # Determine number of options
    noption <- length(option)
    noption.new <- length(option.new)

    # Return value
    if (missing(value)) {
      if (noption.new < noption) {
        if (noption == 1 && option %in% names(default)) {
          return(default[[option]])
        }
        return(NULL)
      }
      if (is.null(which.attr))
        return(dat[[option.new]])
      else
        return(attr(dat[[option.new]], which.attr, exact=TRUE))

    # Set value
    } else {
      if (noption.new == noption || (noption.new == (noption - 1)
          && is.list(if (is.null(option.new)) dat else dat[[option.new]]))) {
        if (is.null(which.attr))
          dat[[option]] <<- value
        else
          attr(dat[[option]], which.attr) <<- value
      }
    }
  }
})
