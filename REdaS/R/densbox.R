densbox <- function(formula, data, rug = FALSE, from, to, gsep = .5, kernel, bw, main, ylab, var_names, box_out = TRUE, horizontal = FALSE, ...){

  if(horizontal) stop('"horizontal = TRUE" not implemented yet.', call. = FALSE)

  if(missing(...)) ddd <- list() else ddd <- list(...)
  for(i in c("gp", "gp.main", "gp.dens", "gp.varn", "gp.xlab", "gp.box")){ if(i %in% names(ddd)){ next } else { ddd[[i]] <- list() } }

  if(missing(kernel)) kernel <- NULL
  kde_settings <- list(kernel = match.arg(kernel, c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine")))
  if(!missing(from)) kde_settings$from <- from
  if(!missing(to)) kde_settings$to <- to
  if(!missing(bw)) kde_settings$bw <- bw

  if(missing(main)) main <- "Density\u00adBox\u00adPlot" else if(length(main) != 1L) stop('"main" must have exactly one element.')
  if(missing(ylab)) ylab <- paste("Values of", deparse(formula[[2]])) else if(length(ylab) != 1L) stop('"ylab" must have exactly one element.')

  ### convert all RHS variables in data to factors
  if(!missing(data)){
    for(i in attributes(terms(formula))$term.labels){
      if(i %in% names(data)) data[[i]] <- as.factor(data[[i]])
    }
  }
  mf           <- if(missing(data)) model.frame(formula, drop.unused.levels = TRUE) else model.frame(formula, data, drop.unused.levels = TRUE)
  y            <- model.response(mf)
  formula[[2]] <- NULL
  if(length(attributes(terms(formula))$order) == 0L){
    x <- data.frame(factor(rep(1, length(y)), labels=""))
    intercept_only <- TRUE
  } else {
    x <- model.frame(formula, data=mf, drop.unused.levels = TRUE)
    intercept_only <- FALSE
  }

  x_vn         <- ncol(x)
  along_x_var  <- seq_len(x_vn)
  y_da         <- tapply(y, x[, rev(along_x_var)], "c", simplify=FALSE)
  data_attr    <- attributes(y_da)
  data_a_dim   <- data_attr$dim
  # variable names        
  if(!missing(var_names)){
    if(length(var_names) == 1L && is.logical(var_names) && var_names){
      var_names <- names(data_attr$dimnames)
    } else if(length(var_names) == 1L && is.logical(var_names) && !var_names){
      var_names <- rep("", length(data_attr$dimnames))
    } else {
      if(length(var_names) == length(data_attr$dimnames)){
        var_names <- rev(var_names)
      } else if(length(var_names) != length(data_attr$dimnames)){
        stop('if "var_names" are specified as a character object, the number of variable names must match the number of variables.')
      }
    }
  } else {
    var_names <- rep("", length(data_attr$dimnames))
  }
  x_all_labels <- lapply(seq_len(x_vn), function(i){ rep(data_attr$dimnames[[i]], length.out=prod(data_a_dim[i:x_vn])) })

  x_at <- seq_along(y_da)
  x_offset <- lapply(seq_along(data_a_dim), function(i){ rep(seq_len(length(y_da)/prod(data_a_dim[1:i])), each=prod(data_a_dim[1:i])) })
  for(i in seq_along(x_offset)){
    x_at <- x_at + gsep * (x_offset[[i]] - 1.0)
  }
  x_offset <- c(list(seq_along(y_da)), x_offset)

  # compute statistics for the plot
  db_stats <- lapply(y_da, function(x){
    if(length(x) < 5){
      list(NULL)
    } else {
      list(
        boxp  = boxplot.stats(x, do.conf = FALSE, do.out = box_out, coef = ifelse(box_out, as.list(args(boxplot.stats))[["coef"]], 0)),
        dens  = do.call("density", modifyList(list(x = x), kde_settings)),
        means = mean(x)
      )
    }
  })

  range_dens <- range(unlist(lapply(db_stats, function(a_list){
    if(!is.null(a_list[[1]])){ range(a_list$dens$x) } else numeric()
  })))
  ymax_dens <- 2*unlist(lapply(db_stats, function(a_list){
    if(!is.null(a_list[[1]])){ max(a_list$dens$y) } else numeric(1)
  }))

  y_ra <- range(c(unlist(y_da), range_dens))
  y_at <- axisTicks(range(unlist(y_ra)), log = FALSE)
  y_ra <- c(min(c(y_ra[1L], y_at[1L])), max(c(y_ra[2L], y_at[length(y_at)])))

  ###  
  ###  
  ###  
  ###  
  ### GRID PLOT VERTICAL
  grid.newpage()

  ### PROCESS OPTIONS (to be implemented & documented ...)
  gp      <- modifyList(gpar(fontsize = 10, lineheight = 1.1), ddd[["gp"]])
  gp.main <- modifyList(gpar(fontface="bold", fontsize = gp$fontsize*1.2), ddd[["gp.main"]])
  gp.dens <- modifyList(gpar(fill="#000000", lty = "blank"), ddd[["gp.dens"]])
  gp.varn <- modifyList(gpar(), ddd[["gp.varn"]])
  gp.ylab <- modifyList(gpar(), ddd[["gp.xlab"]])
  gp.xlab <- modifyList(gpar(), ddd[["gp.xlab"]])
  gp.box  <- modifyList(gpar(fill=gray(.75, .5), col=gray(0)), ddd[["gp.box"]])
  
  
  
  ### main viewport
  pushViewport(viewport(gp = gp, clip = "on"))
  
  

  ### set up heights and widths
    # main title height
    mainH <- convertUnit(grobHeight(textGrob(main, gp = gp.main)) + unit(20, "bigpts"), "bigpts")
    # width of the y-label
    ylabW <- convertUnit(grobWidth(textGrob(ylab, rot = 90)) + unit(10, "bigpts"), "bigpts")
    # y-axis number width
    y_axW <- unit(max(unlist(lapply(y_at, function(temp_y_at){
      convertUnit(grobWidth( textGrob(label = as.expression(temp_y_at), just = c(1, .5), gp = gp.ylab) ), "bigpts")
    }))), "bigpts") + unit(15, "bigpts")

    if(intercept_only){
      xlabH <- unit(2, "bigpts")
    } else {
      xlabH <- unlist(lapply(x_all_labels, function(lab_level){
        max(unlist(lapply(lab_level, function(this_label){
          convertUnit(grobHeight( textGrob(label = this_label, gp = gp.xlab) ), "bigpts")
        })))
      }))
      # control for variable names' heights
      varnH <- unlist(lapply(var_names, function(this_label){
        convertUnit(grobHeight( textGrob(label = this_label, gp = gp.varn) ), "bigpts")
      }))
      #
      xlabH <- unit(apply(cbind(unclass(xlabH), unclass(varnH)), 1L, max), "bigpts")
      # spacers for the tree structure
      if(length(xlabH) == 1L){
        xlabH_single <- c(7.5, xlabH, 0, 2)
      } else if(length(xlabH) > 1L){
        xlabH_single <- c(as.numeric(rbind(15, xlabH)), 0, 2)
        xlabH_single[1L] <- 7.5
      } else {
        xlabH_single <- c(7.5, xlabH, 0, 2)
      }
      xlabH_single <- unit(c(0, xlabH_single), "bigpts")
      xlabH <- unit(sum(unclass(xlabH_single)), "bigpts")
    }

      x_margin      <- ylabW + y_axW + unit(5, "bigpts")
      canvas_width  <- unit(1, "npc") - unit(2, "bigpts") - (ylabW + y_axW + unit(5, "bigpts"))
      y_margin      <- xlabH
      canvas_height <- unit(1, "npc") - (xlabH + mainH)
    
    

  ### set up main canvas
    pushViewport(viewport(
      x      = x_margin,
      y      = y_margin,
      width  = canvas_width,
      height = canvas_height,
      xscale = range(x_at)+c(-.75,.5),
      yscale = y_ra+abs(diff(y_ra))*c(-.05,.05),
      just   = c(0, 0),
      name   = "mainCanvas",
      clip   = "off"
    ))
    
    
  ### x-annotation (really nice tree-structure)
  if(!intercept_only){
    # tree + variable descriptions
    for(i in seq_along(data_a_dim)){
      idx <- seq.int(3L, length(xlabH_single), 2L)
      these_axes <- unlist(tapply(x_at, x_offset[[i]], mean, simplify=FALSE))
      
      if(i == 1L) last_y <- unit(0, "bigpts")
      shift_up <- unit(0, "bigpts")
      if(i == length(data_a_dim)) shift_up <- 0.5*xlabH_single[length(xlabH_single) - 2L]
      
      grid.segments(
        x0 = unit(these_axes, "native"),
        x1 = unit(these_axes, "native"),
        y0 = last_y,
        y1 = last_y <- unit(-sum(unclass(xlabH_single[seq_len(idx[i])])), "bigpts") - 0.5*xlabH_single[idx[i]+1] + shift_up
      )

      join_lines <- apply(matrix(these_axes, nrow=data_a_dim[i]), 2, range)
      for(j in seq_len(ncol(join_lines))){
        grid.segments(
          x0 = unit(join_lines[1L,j], "native"),
          x1 = unit(join_lines[2L,j], "native"),
          y0 = unit(-sum(unclass(xlabH_single[seq_len(idx[i])])), "bigpts") - 0.5*xlabH_single[idx[i]+1] + shift_up,
          y1 = unit(-sum(unclass(xlabH_single[seq_len(idx[i])])), "bigpts") - 0.5*xlabH_single[idx[i]+1] + shift_up
        )
      }
    }

    # label boxes
    for(i in seq_along(data_a_dim)){
      idx <- seq.int(3L, length(xlabH_single), 2L)
      these_axes <- unlist(tapply(x_at, x_offset[[i]], mean, simplify=FALSE))
      for(j in seq_along(these_axes)){
        grid.rect(
          x = unit(these_axes[j], "native"),
          y = unit(-sum(unclass(xlabH_single[seq_len(idx[i]-1)])), "bigpts") - 0.5*xlabH_single[idx[i]],
          width  = grobWidth(textGrob(label=x_all_labels[[i]][j])) + unit(5, "bigpts"),
          height = grobHeight(textGrob(label=x_all_labels[[i]][j])) + unit(5, "bigpts"),
          gp = gpar(col = "white", fill = "white")
        )
      }
    }
    
    # labels
    for(i in seq_along(data_a_dim)){
      idx <- seq.int(3L, length(xlabH_single), 2L)
      these_axes <- unlist(tapply(x_at, x_offset[[i]], mean, simplify=FALSE))

      ### variable descriptions
      grid.text(
        label = var_names[i],
        x     = unit(0, "npc") - x_margin + unit(2, "bigpts"),
        y     = unit(-sum(unclass(xlabH_single[seq_len(idx[i]-1)])), "bigpts") - 0.5*xlabH_single[idx[i]],#unit(-sum(unclass(xlabH_single[seq_len(idx[i]-1L)])), "bigpts") - 0.5*xlabH_single[idx[i]] + shift_up,
        just  = c(0, .5),
        gp    = gp.varn
      )

      grid.text(
        label = x_all_labels[[i]],
        x     = unit(these_axes, "native"),
        y     = unit(-sum(unclass(xlabH_single[seq_len(idx[i]-1)])), "bigpts") - 0.5*xlabH_single[idx[i]],
        just  = c(.5, .5),
        gp    = gp.xlab,
        check.overlap = TRUE,
      )
    }
  }
  


  ### annotations
    # main
    grid.text(
      label = main,
      x = unit(0.5, "npc"),
      y = unit(1, "npc") + 0.5*mainH,
      gp = gp.main
    )
    
    # y-label
    grid.text(
      label = ylab,
      x = unit(0, "npc") - x_margin + 0.5*ylabW,
      y = unit(0.5, "npc"),
      rot = 90
    )

    # y-axes and y-ticks
    grid.text(
      as.expression(y_at),
      x = unit(-7.5, "bigpts"),
      y = unit(y_at, "native"),
      just = c(1, .5),
      check.overlap = TRUE,
      gp = gp.ylab
    )
    grid.segments(
      x0 = unit(0, "npc") - unit(5, "bigpts"),
      x1 = unit(0, "npc"),
      y0 = unit(y_at, "native"),
      y1 = unit(y_at, "native")
    )
  
  
  
  ### data
  # horizontal lines
  for(i in seq_along(y_at)){
    grid.segments(
      x0 = unit(0, "npc"),
      x1 = unit(1, "npc"),
      y0 = unit(y_at[i], "native"),
      y1 = unit(y_at[i], "native"),
      gp = gpar(lty=2, col=gray(.875)))
  }
  # cycle through all groups
  for(i in seq_along(x_at)){
    # vertical lines
    grid.segments(
      x0 = unit(x_at[i], "native"),
      x1 = unit(x_at[i], "native"),
      y0 = unit(0, "npc"),
      y1 = unit(1, "npc"),
      gp = gpar(lty = 3, col = gray(.750))
    )

    # if the current group has >= 5 obs.
    if(!is.null(db_stats[[i]][[1L]])){
      # boxes
      ystats <- c(db_stats[[i]]$boxp$stats, db_stats[[i]]$means)
      
      # box-background
      grid.rect(
        x      = unit(x_at[i], "native"),
        y      = unit(ystats[2L], "native"),
        width  = unit(.25, "native"),
        height = unit(abs(diff(ystats[c(2L, 4L)])), "native"),
        just   = c(0, 0),
        gp     = modifyList(gp.box, list(lty = "blank"))
      )

      # rug?
      if(rug){
        grid.segments(
          x0 = unit(x_at[i], "native"),
          x1 = unit(x_at[i] + .1, "native"),
          y0 = unit(y_da[[i]], "native"),
          y1 = unit(y_da[[i]], "native"),
          gp = gpar(col=gray(0, .5), lwd=.5)
        )
      }

      # min, median, max
      grid.segments(
        x0 = unit(x_at[i], "native"),
        x1 = unit(x_at[i] + c(0,.2,.25,.2), "native"),
        y0 = unit(ystats[c(1L, 1L, 3L, 5L)], "native"),
        y1 = unit(ystats[c(5L, 1L, 3L, 5L)], "native"),
        gp = gpar(col=gray(0), lwd=1)
      )
      # box-border
      grid.rect(
        x      = unit(x_at[i], "native"),
        y      = unit(ystats[2L], "native"),
        width  = unit(.25, "native"),
        height = unit(abs(diff(ystats[c(2L, 4L)])), "native"),
        just   = c(0, 0),
        gp     = modifyList(gp.box, list(fill = "transparent"))
      )
      # mean
      grid.segments(
        x0 = unit(x_at[i], "native"),
        x1 = unit(x_at[i] + .3, "native"),
        y0 = unit(ystats[6L], "native"),
        y1 = unit(ystats[6L], "native"),
        gp = gpar(lty=2, col=gray(0), lwd=1)
      )
      # outliers
      if((y_box_out_length <- length(db_stats[[i]]$boxp$out)) > 0){
        grid.points(
          x    = unit(rep(x_at[i], y_box_out_length), "native") + unit(4, "bigpts"),
          y    = unit(db_stats[[i]]$boxp$out, "native"),
          pch  = 8L,
          size = unit(4, "bigpts"),
          gp   = gpar(col = gray(0))
        )
      }

      # densities
      dens_len <- length(db_stats[[i]]$dens$x)
      grid.polygon(
        x = -c(0, db_stats[[i]]$dens$y, 0)/ymax_dens[i] + x_at[i],
        y = db_stats[[i]]$dens$x[c(1, seq_len(dens_len), dens_len)],
        gp = gp.dens,
        default.units="native"
      )
    # if the current group has < 5 obs.
    } else {
      grid.text(
        label = expression(""<5),
        x     = unit(x_at[i], "native"),
        y     = unit(.5, "npc")
      )
      warning("at least one group with < 5 observations was present.")
    }
  }

  ### finalize
  grid.rect()
  popViewport()

}
