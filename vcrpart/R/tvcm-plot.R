##' -------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          rbuergin@gmx.ch
##' Date:            2016-02-16
##'
##' Description:
##' Plot functions for 'tvcm' objects.
##'
##' References:
##' party:           http://CRAN.R-project.org/package=party
##' partykit:        http://CRAN.R-project.org/package=partykit
##'
##' Overview:
##' plot.tvcm:       generic plot for tree.tvcm objects
##' panel_partdep:   partial coefficient plots
##' panel_get_main:  extracts the title for each terminal node
##' panel_coef:      grapcon generator for coefficient plots
##' panel_empty:     grapcon generator for empty terminal node plots
##'
##' Last modifications:
##' 2016-02-16: modified titles for 'panel_coef'.
##' 2016-02-08: add warning in cases 'conf.int = TRUE' in 'panel_coef'.
##' 2014-09-08: replace 'do.call' by 'eval'
##' 2014-09-06: solve bugs in 'panel_partdep'
##' 2014-09-06: - add 'type = "cv"' option for the cases where
##'               cross validation was incorporated in the
##'               partitioning stage
##' 2014-08-30: - automatic titles
##'             - write '(no split)' in panels where no split
##'               is applied
##' 2014-08-27: added the 'yadj' argument for 'panel_coef'
##' 2014-07-02: make a loop over the partitions
##' 2014-07-01: add checks for each function
##' -------------------------------------------------------- #

plot.tvcm <- function(x, type = c("default", "coef", 
                           "simple", "partdep", "cv"),
                      main, part = NULL,
                      drop_terminal = TRUE,
                      tnex, newpage = TRUE, ask = NULL, 
                      pop = TRUE, gp = gpar(), ...) {

    ## checks
    type <- match.arg(type)
    stopifnot(is.logical(drop_terminal) && length(drop_terminal) == 1L)
    if (!missing(tnex)) stopifnot(is.numeric(tnex) && length(tnex) == 1L)
    stopifnot(is.logical(newpage) && length(newpage) == 1L)
    stopifnot(is.logical(pop) && length(pop) == 1L)
    stopifnot(class(gp) %in% c("gpar", "list"))
    
    if (type == "partdep") {
        
        if (missing(main)) main <- NULL
        call <- list(as.name("panel_partdep"),
                     object = quote(x), ask = ask, main = main)
        call <- append(call, list(...))
        mode(call) <- "call"
        eval(call)
        
    } else if (type == "cv") {
        
        if (missing(main)) main <- NULL
    if (is.null(x$info$cv)) {
        warning("no information on cross validation.")
    } else {
        plot(x$info$cv, main = main, ...)
    }
        
    } else {
        
        ## tree plots
        if (is.null(part)) part <- seq_along(x$info$node)
        if (is.character(part)) {
            part <- which(LETTERS[seq_along(x$info$node)] %in% part)
        } else if (is.numeric(part)) {
            part <- as.integer(part)
        } else {
            stop("'part' must be a 'character' or a 'integer'.")
        }
        if (length(part) < 1L) stop("no valid 'part' specified.")
        
        ## whether an input is expected before plotting the next tree
        if (is.null(ask))
            ask <- ifelse(length(part) == 1L, FALSE, TRUE)
        
        ## repeat the title
        if (missing(main)) {
            main <- tvcm_print_vclabs(x$info$formula, TRUE)[part]
        } else if (is.character(main)){
            main <- rep(main, length.out = length(part))
        } else {
            main <- NULL
        }
        
        ## terminal panel
        tp_fun <-
            switch(type,
                   "default" =
                       if (max(sapply(x$info$node[part], width))> 4L) {
                           panel_empty
                       } else {
                           panel_coef
                       },
                   "coef" = panel_coef,
                   "simple" = panel_empty)
        tp_args <- if ("tp_args" %in% names(list(...))) {
            list(...)$tp_args
        } else {
            list(...)[names(list(...)) %in% names(formals(tp_fun))[-1]]
        }
        tp_args$part <- part
        tp_args$gp <- gp
        
        ## inner panel
        inner_panel <- if ("inner_panel" %in% names(list(...))) {
            list(...)$inner_panel
        } else {
            switch(type,
                   "default" = node_inner,
                   "coef" = node_inner,
                   "simple" = node_inner)
        }
        ip_args <- if ("ip_args" %in% names(list(...))) {
            list(...)$ip_args
        } else {
            list(...)[names(list(...)) %in% names(formals(inner_panel))[-1]]
        }
        
        ## edge panel
        edge_panel <- if ("edge_panel" %in% names(list(...))) {
            list(...)$edge_panel
        } else {
            edge_default
        }
        ep_args <- if ("ep_args" %in% names(list(...))) {
            list(...)$ep_args
        } else {
            list(...)[names(list(...)) %in% names(formals(edge_panel))[-1]]
        }
        if ("justmin" %in% names(formals(edge_panel)) &&
            !"justmin" %in% names(ep_args)) ep_args$justmin <- 5
        
        ## other arguments
        dotargs <- list(...)[names(list(...)) %in%
                             names(formals(plot.party))[-1]]
        dotargs <- dotargs[setdiff(names(dotargs),
                                   c("tp_args", "ip_args", "ep_args"))]
        
        ## prepare call
        call <- list(name = as.name("plot.party"),
                     x = quote(x),
                     terminal_panel = quote(tp_fun),
                     tp_args = quote(tp_args),
                     inner_panel = quote(inner_panel),
                     ip_args = quote(ip_args),
                     edge_panel = quote(edge_panel),
                     ep_args = quote(ep_args),
                     drop_terminal = quote(drop_terminal),
                     tnex = quote(tnexCall),
                     newpage = quote(newpage),
                     main = quote(main[pid]),
                     pop = pop, gp = gp)
        call <- append(call, dotargs)    
        mode(call) <- "call"
        
        if (ask) {
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))
        }
        
        ## call
        for (pid in seq_along(part)) {
            tp_args$part <- part[pid]
            x$node <- x$info$node[[part[pid]]]
            tnexCall <- if (missing(tnex)) max(1, depth(x$node)) else tnex
            eval(call)
        }
    }
}
                    

panel_partdep <- function(object, parm = NULL,
                          var = NULL, ask = NULL,
                          prob = NULL, neval = 50, add = FALSE,
                          etalab = c("int", "char", "eta"), ...) {

  ## set and check 'parm'
  allParm <- unique(unlist(tvcm_get_vcparm(object)))
  if (is.null(parm)) parm <- allParm
  if (!all(parm %in% allParm))
    warnings("some 'parm' were not recognized.")
  parm <- intersect(parm, allParm)
  if (length(parm) == 0L) stop("no valid 'parm'.")
  stopifnot(is.null(ask) | is.logical(ask) && length(ask) == 1L)
  stopifnot(is.null(prob) | (is.numeric(prob) && length(prob) == 1L))
  stopifnot(prob > 0 & prob <= 1)
  stopifnot(is.numeric(neval) && length(eval) == 1L && neval > 0)
  stopifnot(is.logical(add) && length(add) == 1L)
  etalab <- match.arg(etalab)
  
  ## labels for coefficients
  termLabs <- parm
  if (object$info$fit == "olmm") 
    termLabs <- olmm_rename(parm, levels(object$info$model$y),
                            object$info$family, etalab)
  
  ## set and check 'var'
  data <- object$data
  allVars <- colnames(data)
  if (is.null(var)) var <- allVars[1L]
  if (!all(var %in% allVars))
    warnings("some variables in 'var' were not recognized.")
  var <- intersect(var, allVars)
  if (length(var) == 0L) stop("no valid variables in 'var'.")
  
  ## set defaults
  if (is.null(ask))
    ask <- ifelse(length(parm) * length(var) == 1L, FALSE, TRUE)
  
  ## dotlist
  dotList <- list(...)

  ## get value spaces
  space <- vcrpart_value_space(data, neval)

  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  for (i in 1:length(var)) {

    z <- space[[var[i]]]

    ## draw a random subsample to save time
    if (is.null(prob)) {
      if (length(z) * nrow(data) > 10000) {
        probi <- 10000 / (length(z) * nrow(data))
      } else {
        probi <- 1.0
      }
    } else {
      probi <- prob
    }
    tmp <- if (probi < 1) { data[tvcm_folds(object, folds_control(K = 1L, type = "subsampling", prob = probi)) == 1L, , drop = FALSE] } else { data }

    newdata <- NULL
    for (j in 1:length(z)) newdata <- rbind(newdata, tmp)
    newdata[, var[i]] <- rep(z, each = nrow(tmp))
    beta <- predict(object, newdata = newdata, type = "coef")
    parm <- intersect(parm, colnames(beta))
    
    ## for each coefficients ...
    for (j in 1:length(parm)) {

      ## .. compute the expectation of beta_i at a value z_j
      ## over the distribution  of z_-j
      betaj <- tapply(beta[, parm[j]], newdata[, var[i]], mean)

      ## define arguments for plot or points call
      ylab <-
        paste("coef(", termLabs[j], "|", var[i], ")", sep = "")

      if (is.numeric(z)) {
        call <- appendDefArgs(dotList, list(type = "s"))
        call <- appendDefArgs(list(x = z, y = betaj), call)
      } else {
        call <- dotList
        call <- appendDefArgs(list(x = as.integer(z), y = betaj,
                                            type = "b"), call)
        if (!add) call <- appendDefArgs(call, list(axes = FALSE))
      }
      
      call <- appendDefArgs(call, list(xlab = var[i], ylab = ylab))
      call <- append(list(name = as.name(ifelse(add, "points", "plot.default"))), call)
      mode(call) <- "call"
      
      ## call plot or points
      eval(call)
      if (!add && is.factor(z) && !"axes" %in% names(list(...)) &&
          !call$axes) {
        box()
        axis(1, 1:nlevels(z), levels(z))
        axis(2)
      }
    }
  }
}


##' -------------------------------------------------------- #
##' Specifies a title for a terminal node
##'
##' @param object a \code{\link{tvcm}} object.
##' @param node   a \code{\link{partynode}} object. The
##'    terminal node for which the title is to be specified.
##' @param id     a logical scalar. Whether the id of the node
##'    is to be included in the title.
##' @param nobs   a logical scalar. Whether the number of
##'    observations in the node is to be included in the
##'    title.
##'
##' @return a character string.
##' -------------------------------------------------------- #

panel_get_main <- function(object, node, id, nobs) {
  rval <- ""
  if (id) rval <-
    paste(rval, paste0(names(object)[id_node(node)]))
  if (id && nobs) rval <- paste0(rval, ":")
  if (nobs) rval <- paste0(rval, "n=", node$info$dims["n"])
  return(rval)
}


panel_coef <- function(object, parm = NULL, 
                       id = TRUE, nobs = TRUE,
                       exp = FALSE,
                       plot_gp = list(),
                       margins, yadj = 0.1,
                       mean = FALSE, mean_gp = list(),
                       conf.int = FALSE, conf.int_gp = list(),
                       abbreviate = TRUE,
                       etalab = c("int", "char", "eta"),
                       ...) {

  ## checks
  stopifnot(is.logical(id) && length(id) == 1L)
  stopifnot(is.logical(nobs) && length(nobs) == 1L)
  stopifnot(is.list(plot_gp))
  if (!missing(margins)) stopifnot(is.numeric(margins) && length(margins) == 4L)
  stopifnot(is.list(mean_gp))
  stopifnot(is.list(conf.int_gp))
  if (conf.int) warning("The shown confidence intervals do not ",
                        "account for the model selection procedure. ",
                        "Therefore interpret them carefully.") 
  
  ## get partition
  part <- which(sapply(object$info$node, identical, object$node))
  
  ## get arguments for setting y tick labels
  etalab <- match.arg(etalab)
  if (object$info$fit == "olmm") {
    yLevs <- levels(object$info$model$y)
  } else {
    yLevs <- all.vars(object$info$formula$original)[1L]
  }
  
  ## extract coefficients
  coef <- extract(object, "coef")$vc[[part[1L]]]
  if (length(na.omit(coef)) < 1L) {
      rval <- function(node) {
          pushViewport(viewport())
          grid.text("(no split)")
          upViewport(2L)
      }
      return(rval)
  }

  ## extract subset of coefficients
  if (is.null(parm)) parm <- colnames(coef)
  if (!is.list(parm)) parm <- list(parm)
  if (is.numeric(unlist(parm))) {
    if (!all(unlist(parm) %in% seq_along(colnames(coef))))
      stop("'parm' must be a subset of ", 1, " to ", ncol(coef), ".")
    parm <- lapply(parm, function(x) colnames(coef)[x])
  } else if (is.character(unlist(parm))) {
    if (!all(unlist(parm) %in% colnames(coef)))
      stop("'parm' must be a subset of ",
           paste(paste("'", colnames(coef), "'", sep = ""), collapse = ", "),
           ".")
  } else {
    stop("'parm' must be either an integer or a character vector.")
  }
  coef <- coef[, unlist(parm), drop = FALSE]
  coefList <- lapply(parm, function(trms) coef[, trms, drop = FALSE])
  
  if (conf.int) {
    sd <- extract(object, "sd")$vc[[part[1L]]]
    sd <- sd[, unlist(parm), drop = FALSE]
    sdList <- lapply(parm, function(trms) sd[, trms, drop = FALSE])
  }

  argsToList <- function(x) {
    if (is.null(names(x))) {
      if (length(x) < length(parm))
      x <- append(x, vector("list", length(parm) - length(x)))
    } else {
      x <- lapply(seq_along(parm), function(i) x)
    }
    return(x)
  }

  ## determine 'margins' if required
  if (missing(margins)) {
      margins <- c(0.5, 0.5, 0, 0)
      if (length(parm) > 1) margins[1L] <- 1.5
  }
  
  ## whether y labels should be the exponential of their actual values
  exp <- rep(exp, length.out = length(parm))
  
  ## plot parameters
  plot_gp <- argsToList(plot_gp)
  plot_gp <-
      lapply(seq_along(coefList),
           function(i) {
             if (conf.int) {
               ylim <- range(c(c(coefList[[i]] - 2 * sdList[[i]]),
                               c(coefList[[i]] + 2 * sdList[[i]])),
                             na.rm = TRUE)
             } else {
               ylim <- range(coefList[[i]], na.rm = TRUE)
             }
             ylim <- range(c(0, ylim))
             ylim <- ylim + c(-1, 1) * 0.1 * diff(ylim)
             xlabel <- colnames(coefList[[i]])
             if (object$info$fit == "olmm")
               xlabel <- olmm_rename(xlabel, yLevs, object$info$family, etalab)
             xlabel <- unlist(lapply(xlabel,
                                     function(x) {
                                       x <- unlist(strsplit(x, ":"))
                                       if (abbreviate) x <- abbreviate(x)
                                       return(paste(x, collapse = ":"))
                                     }))
             rval <- list(xlim = c(0.75, ncol(coefList[[i]]) + 0.25),
                          pch = 1L, ylim = ylim,
                          ylab = if (exp[i]) "exp(coef)" else "coef",
                          type = "b",
                          xlabel = xlabel, 
                          height = 1,
                          width = min(0.6, 0.15 * width(object$node)), 
                          gp = gpar(cex = 1L, fontsize = 8))
             if (length(plot_gp) > 0L)
               rval <- appendDefArgs(plot_gp[[i]], rval)
             rval$yat <- axisTicks(rval$ylim, log = FALSE, nint = 3)
             rval$ylabel <- if (exp[i]) exp(rval$yat) else rval$yat
             rval$ylabel <- format(rval$ylabel, digits = 2L)
             return(rval)
           })

  if (conf.int) {

    conf.int_gp <- argsToList(conf.int_gp)
    conf.int_gp_def <- list(angle = 90,
                            length = unit(1, "mm"),
                            ends = "both",
                            type = "open")
    conf.int_gp <- lapply(seq_along(parm),
                          function(i) appendDefArgs(conf.int_gp[[i]],
                                                    conf.int_gp_def))
  }
  
  ## population mean
  meanCoef <- NULL
  if (mean) {

    mean_gp <- argsToList(mean_gp)
    mean_gp_def <- list(gp = gpar(col = "grey50", lty = 1,
                          lwd = 0.75, fontsize = 8), pch = 1L)
    mean_gp <- lapply(seq_along(parm),
                      function(i) appendDefArgs(mean_gp[[i]], mean_gp_def))

    w <- tapply(weights(object), predict(object, type = "node")[, part], sum) /
      sum(weights(object))
    meanCoef <- colSums(coef * matrix(w, nrow(coef), ncol(coef)))
    meanCoef <- lapply(parm, function(trms) meanCoef[trms, drop = FALSE])
  }

  qN <- qnorm(0.975)
  
  rval <- function(node) {

    strUnit <- if (any(c(id, nobs))) unit(2, "strheight", "A") else unit(0, "npc")
    pushViewport(viewport(layout = grid.layout(3, 1,
                              heights = unit.c(unit(yadj, "npc"),
                                  strUnit,
                                  unit(1, "npc") - unit(yadj, "npc") - strUnit)))) # 1
    
    pushViewport(viewport(layout.pos.row = 2L)) # 2
    grid.rect(gp = gpar(fill = "white", col = 0))
    grid.text(panel_get_main(object, node, id, nobs))
    upViewport()
    
    pushViewport(viewport(layout.pos.row = 3L)) # 2
    grid.rect(gp = gpar(fill = "white", col = 0))
    pushViewport(viewport(layout = grid.layout(length(coefList), 1))) # 3
    
    for (i in 1:length(coefList)) {
      
      pushViewport(viewport(layout.pos.row = i)) # 4
      
      pushViewport(viewport(height = unit(plot_gp[[i]]$height, "npc"),
                            width = unit(plot_gp[[i]]$width, "npc"))) # 5
      pushViewport(plotViewport(margins = margins,
                                xscale = plot_gp[[i]]$xlim,
                                yscale = plot_gp[[i]]$ylim,
                                default.units = "native")) # 6

      ## vertical lines
      grid.segments(unit(1:ncol(coefList[[i]]), "native"),
                    unit(rep(plot_gp[[i]]$ylim[1], ncol(coefList[[i]])), "native"),
                    unit(1:ncol(coefList[[i]]), "native"),
                    unit(rep(plot_gp[[i]]$ylim[2], ncol(coefList[[i]])), "native"),
                    gp = gpar(col = "lightgrey", lty = 3))

      if (plot_gp[[i]]$ylim[1] < 0.0 & plot_gp[[i]]$ylim[2] > 0)
      grid.segments(unit(plot_gp[[i]]$xlim[1], "native"), unit(0, "native"),
                    unit(plot_gp[[i]]$xlim[2], "native"), unit(0, "native"),
                    gp = gpar(col = "black"))

      subs <- coefList[[i]][as.character(id_node(node)),] >= plot_gp[[i]]$ylim[1L] &
          coefList[[i]][as.character(id_node(node)),]<= plot_gp[[i]]$ylim[2L]

      subsMean <- NULL
      if (mean) {
          subsMean <- meanCoef[[i]] > plot_gp[[i]]$ylim[1L] &
              meanCoef[[i]] > plot_gp[[i]]$ylim[1L]
      }

      ## option 'conf.int = TRUE'
      if  (conf.int) {

        ## crop the lines
        nCoef <- length(coefList[[i]][as.character(id_node(node)),])
        endCi <- rep(conf.int_gp[[i]]$ends, length.out = nCoef)
        lenCi <- rep(conf.int_gp[[i]]$length, length.out = nCoef)        
        lwr <- coefList[[i]][as.character(id_node(node)),] -
          qN * sdList[[i]][as.character(id_node(node)),]
        endCi[lwr < plot_gp[[i]]$ylim[1L]] <- "last"
        lwr[lwr < plot_gp[[i]]$ylim[1L]] <- plot_gp[[i]]$ylim[1L]
        lwr[lwr > plot_gp[[i]]$ylim[2L]] <- NA
        upr <- coefList[[i]][as.character(id_node(node)),] +
          qN * sdList[[i]][as.character(id_node(node)),]
        endCi[upr > plot_gp[[i]]$ylim[2L] & endCi == "last"] <- "none"
        endCi[upr > plot_gp[[i]]$ylim[2L] & endCi == "both"] <- "first"
        upr[upr > plot_gp[[i]]$ylim[2L]] <- plot_gp[[i]]$ylim[2]
        upr[upr < plot_gp[[i]]$ylim[1L]] <- NA
        subsCi <- !is.na(lwr) & !is.na(upr)
        lenCi[endCi == "none"] <- 0
        endCi[endCi == "none"] <- "both"
        
        ## plot
        if (any(subsCi))
          grid.segments(unit(which(subsCi), "native"),
                        unit(lwr[subsCi], "native"),
                        unit(which(subsCi), "native"),
                        unit(upr[subsCi], "native"),
                        arrow = arrow(angle = conf.int_gp[[i]]$angle,
                          length = lenCi, 
                          ends = endCi,
                          type = conf.int_gp[[i]]$type),
                        gp = plot_gp[[i]]$gp)
        
      }
      
      ## option 'type = "p"'
      if (plot_gp[[i]]$type %in% c("p", "b")) {
          
        if (mean && any(subsMean))
          grid.points(unit(which(subsMean), "native"),
                      unit(meanCoef[[i]][subsMean], "native"),
                      pch = mean_gp[[i]]$pch, gp = mean_gp[[i]]$gp)
        
        if (any(subs)) 
          grid.points(unit(which(subs), "native"),
                      unit(coefList[[i]][as.character(id_node(node)),][subs],
                           "native"),
                      pch = plot_gp[[i]]$pch, gp = plot_gp[[i]]$gp)
      }
      
      ## option 'type = "l"'
      if (plot_gp[[i]]$type %in% c("l", "b")) {
        
        if (mean && any(subsMean)) 
          grid.lines(unit(which(subsMean), "native"),
                     unit(meanCoef[[i]][subsMean], "native"),
                     gp = mean_gp[[i]]$gp)

        if (any(subs))
            grid.lines(unit(which(subs), "native"),
                       unit(coefList[[i]][as.character(id_node(node)),][subs],
                            "native"),
                       gp = plot_gp[[i]]$gp)
      }
      
      if (id_node(node) == min(nodeids(object, terminal = TRUE))) {
        grid.yaxis(at = plot_gp[[i]]$yat,
                   label = plot_gp[[i]]$ylabel,
                   gp = gpar(lineheight = 0.2))
        grid.text(plot_gp[[i]]$ylab,
                  unit(-0.75 - 0.6 * max(nchar(plot_gp[[i]]$ylabel)), "char"),
                  unit(0.5, "npc"), rot = 90)
      }
      
      grid.xaxis(at = 1:ncol(coefList[[i]]),
                 label = plot_gp[[i]]$xlabel,
                 gp = gpar(lineheight = 0.4))
      grid.rect()
      upViewport(3L)

    }
    
    ## close viewports
    upViewport(3L)
    return(rval)
  }
}
class(panel_coef) <- "grapcon_generator"

panel_empty <- function(object, part = 1L, id = TRUE, nobs = TRUE, ...) {
  
  stopifnot(is.logical(id) && length(id) == 1L)
  stopifnot(is.logical(nobs) && length(nobs) == 1L)
  
  rval <- function(node) {
    
    nid <- id_node(node)
    pushViewport(viewport(height = unit(3, "strheight", "A")))
    grid.rect(gp = gpar(fill = "white", col = 0))
    grid.text(panel_get_main(object, node, id, nobs))
    upViewport()
  }
  return(rval)
}
class(panel_empty) <- "grapcon_generator"

edge_default <- function(obj, digits = 3, abbreviate = FALSE,
                         justmin = Inf,
                         just = c("alternate", "increasing", "decreasing", "equal")) {
  meta <- obj$data
  
  justfun <- function(i, split) {
    myjust <- if(mean(nchar(split)) > justmin) {
      match.arg(just, c("alternate", "increasing", "decreasing", "equal"))
    } else {
      "equal"
    }
    k <- length(split)
    rval <- switch(myjust,
                   "equal" = rep.int(0, k),
                   "alternate" = rep(c(0.5, -0.5), length.out = k),
                   "increasing" = seq(from = -k/2, to =  k/2, by = 1),
                   "decreasing" = seq(from =  k/2, to = -k/2, by = -1)
                   )
    unit(0.5, "npc") + unit(rval[i], "lines")
  }
  
  ## panel function for simple edge labelling
  function(node, i) {
    split <- character_split(split_node(node), meta, digits = digits)$levels
    y <- justfun(i, split)
    split <- split[i]
    ## try() because the following won't work for split = "< 10 Euro", for example.
    if(any(grep(">", split) > 0) | any(grep("<", split) > 0)) {
      tr <- suppressWarnings(try(parse(text = paste("phantom(0)", split)),
                                 silent = TRUE))
      if(!inherits(tr, "try-error")) split <- tr
    }
    grid.rect(y = y, gp = gpar(fill = "white", col = 0),
              width = unit(1, "strwidth", split))
    grid.text(split, y = y, just = "center")
  }
}
class(edge_default) <- "grapcon_generator"
