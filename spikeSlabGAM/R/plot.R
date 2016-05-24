#' @include spikeSlabGAM.R
{}

#' Plot the estimated effect of a model term.
#'
#' Plots the estimated linear predictor for a model term, i.e a main effect or a
#' (2- or 3-way) interaction. Plots for the joint effect of two numerical
#' covariates show an overlay if quantiles were specified: Regions where the
#' pointwise credible intervals do not contain zero are plotted in muted red
#' (\eqn{>0}) and blue (\eqn{< 0}), overlaid by coloured contour lines that show
#' the \code{aggregate} values.  Contour lines are shown only inside the convex
#' hull of the original observations. Plots for
#' \code{\link{srf}}:\code{\link{lin}} terms show the spatially varying
#' coefficient, i.e. the contour lines represent the change in the linear
#' predictor when the \code{lin}-covariate increases by 1 standard deviation.
#' For this reason, a cumulative plot makes no sense and the routine will set
#' \code{cumulative = FALSE} with a warning.
#'
#' Limitations: Plots for 4-way (or higher) interactions are not implemented.
#' Requesting such a plot will return the NULL object with a warning. Plots for
#' \code{\link{mrf}} just treat the grouping variable as a conventional factor
#' i.e. will not incorporate any neighborhood or map information.
#'
#' @param label (character) the label of the effect to be plotted.
#' @param m a fitted \code{spikeSlabGAM} model
#' @param cumulative Defaults to TRUE, in which case the lower order terms that
#'   are associated with the covariates involved in \code{label} are cumulated
#'   and then plotted. (e.g, if \code{label} denotes a smoooth term, the sum of
#'   the linear and smooth effect is displayed if TRUE, if \code{label} is a
#'   factor-factor interaction, the sum of both main effects and their
#'   interaction is displayed etc.) If FALSE, only the marginal effect of
#'   \code{label} is displayed.
#' @param aggregate (function) which summary statistic to use for the posterior
#'   of the effect. Defaults to the mean.
#' @param quantiles which quantiles to use for the borders of credible regions.
#'   Defaults to 10% and 90% percentiles. Set to NULL to omit plotting credible
#'   regions. Cannot deal with more than two quantiles.
#' @param gridlength length of the (univariate) grids on which to evaluate the
#'   posterior.
#' @param contours use how many contour lines for the joint effect of two
#'   numerical covariates. Defaults to 30.
#' @param ggElems a list of plot elements to give to
#'   \code{\link[ggplot2]{ggplot}}. Use this to supply custom themes, colour
#'   schemes etc. Defaults to an empty list, which yields the standard settings.
#' @export
#' @return an object of class \code{ggplot}. Use \code{print} or wrap the call
#'   to \code{plotTerm} in parentheses to directly render on a graphic device.
#' @author Fabian Scheipl
#' @examples
#' #see help for spikeSlabGAM
#' @import ggplot2
#' @import reshape
#' @importFrom scales muted
#' @import grDevices
plotTerm <- function(label, m, cumulative = TRUE,
  aggregate = mean, quantiles = c(.1, .9), gridlength = 40, contours = 30,
  ggElems = list()) {

  nBins <- contours
  alphaRug <- .2

  plotCR <- !is.null(quantiles)

  stopifnot(label %in% names(m$predvars),
    any(!plotCR, length(quantiles)== 2),
    class(m)=="spikeSlabGAM",
    is.function(aggregate),
    all(is.numeric(gridlength), gridlength>0))

  quantiles <- sort(quantiles)


  # find main effect(s) and add interactions
  trms <- allInvolvedTerms(label, m)


  #find variables that are involved...
  vars <- sapply(trms, function(x) gsub("\\)","", gsub("([[:graph:]]+\\()","", x)))
  xu <- lapply(m$predvars[trms], function(x) environment(x$makeNewX)$xu)
  names(xu) <- vars
  vars <- unique(vars)
  xu <- xu[!duplicated(names(xu))]


  if(length(vars)>3) {
    warning("You'll have to figure that one out by yourself.",
      "Plots for 4-way or higher interactions are not implemented, returning NULL.")
    return(NULL)
  }
  isMrf <- grepl("mrf(", trms, fixed = TRUE)
  if(any(isMrf)) {
    warning("Default plots for mrf() display the MRF as a conventional factor.")
  }

  isSrf <- grepl("srf(", trms, fixed = TRUE)& !grepl(":", trms)
  if(any(isSrf)) {
    #numeric-srf interaction only possible for lin()
    if(any(grepl("sm(", trms, fixed = TRUE))) {
      warning("Plots for srf():sm() not implemented, returning NULL.")
      return(NULL)
    }
    isFct <- grepl("fct(", trms, fixed = TRUE) & !grepl(":", trms)
    isLin <- grepl("lin(", trms, fixed = TRUE) & !grepl(":", trms)
    varsLin <- sapply(trms[isLin], function(x) gsub("\\)","",
      gsub("([[:graph:]]+\\()","", x)))
    varsSrf <- sapply(trms[isSrf], function(x) gsub("\\)","",
      gsub("([[:graph:]]+\\()","", x)))

    if(any(isLin)) {
      # TODO: check degree = 1
      if(any(cumulative)) {
        warning("Cumulative plot for a space-varying linear effect ", label,
          " makes no sense. Setting cumulative = FALSE.")
        cumulative <- FALSE
      }
    }



    dfs <- sapply(xu, is.data.frame)
    if(sum(dfs)> 1) {
      warning("Plots for srf():srf() not implemented, returning NULL.")
      return(NULL)
    }
    srfGrid <- {
      tmp <-  xu[dfs][[1]]
      uvals <- apply(tmp, 2, function(x) seq(min(x), max(x), l = gridlength/2))
      tmp <-  expand.grid(uvals[, 1], uvals[, 2])
      rng <- environment(m$predvars[trms[isSrf][1]][[1]]$makeNewX)$rng
      insidePoly <- environment(m$predvars[trms[isSrf][1]][[1]]$makeNewX)$insidePoly
      tmp <- tmp[insidePoly(tmp, rng),]
      names(tmp) <- names(xu[dfs][[1]])
      tmp
    }

    #sort & make grid for factor variables:
    facGrid <- if(any(isFct)) {
      do.call(expand.grid, lapply(xu[!dfs & !(names(xu)%in%varsLin)],
        function(x) return(sort(x))))
    } else NULL

    linGrid <- if(any(isLin)) {
      tmp <- as.list(rep(.5, length(varsLin)))
      names(tmp) <- varsLin
      data.frame(tmp)
    } else NULL

    if(NCOL(linGrid)>1) {
      warning("Plots for 3 way numeric interactions with srf() not implemented.",
        "Returning NULL.")
      return(NULL)
    }

    grid <- if(any(isLin) | any(isFct)) {
      tmp <- do.call(expand.grid, c(as.list(linGrid), as.list(facGrid)))
      reps <- nrow(tmp)
      ind <- rep(1:reps, e = nrow(srfGrid))
      tmp <- data.frame(tmp[ind,])
      names(tmp) <- c(names(linGrid), names(facGrid))
      srfGrid <- do.call(rbind, replicate(reps, list(srfGrid)))
      eval(substitute(tmp$.b <- srfGrid, list(.b = names(xu)[dfs])))
      tmp
    } else {
      tmp <- data.frame(FILLER = rep(NA, nrow(srfGrid)))
      eval(substitute(tmp$.b <- srfGrid, list(.b = names(xu)[dfs])))
      tmp
    }

    # get predictions
    pred <- data.frame(evalTerm(label, model = m, newdata = grid,
      aggregate = aggregate, quantiles = quantiles, returnData = FALSE))
    fcts <- sapply(grid[names(grid) !="FILLER"], is.factor)
    nums <- sapply(grid[names(grid) !="FILLER"],
      function(x) !is.factor(x)&!is.data.frame(x))
    data <- data.frame(pred, grid, srfGrid)

    quantileStr <- colnames(pred)[-1]
    aggregateStr <- colnames(pred)[1]


    titleStr <- if(cumulative) paste(vars, collapse ="*") else label
    plotElems <- list(labs(title = titleStr))

    if(length(vars)== 1)	{
      # srf
      nms <-  names(xu[[1]])
      aesStr <- aes_string(x = nms[1], y = nms[2], z = aggregateStr,
        colour ="..level..")
      plotElems <- c(geom_contour(data = data, mapping = aesStr, bins = nBins),
        list(scale_colour_gradient2(name = expression(eta),
          low =  muted('darkblue'), mid = 'grey80', high = muted('darkred'))),
        geom_point(data = m$data, mapping = aes_string(x = nms[1], y = nms[2]),
          alpha = alphaRug, size =.5),
        plotElems)
      if(plotCR) {
        data$sign <- ifelse(apply(data[, 2:3], 1, prod)>0, sign(data[, 3]), 0)
        plotElems <- c(geom_tile(data = data,
          aes_string(x = nms[1], y = nms[2], fill ="sign"), colour = NA, alpha =.1),
          list(scale_fill_gradient2(name ="CI", low =  muted(alpha('darkblue', .1)),
            mid = 'white', high = muted(alpha('darkred', .1)),
            breaks = c(-1, 1), limits = c(-1.1, 1.1),  labels = c("< 0","> 0"))),
          plotElems)

      }
    }
    if(length(vars)== 2 & any(fcts)) {
      #srf:fct
      nms <- paste(paste(names(xu)[dfs],"$", sep =""), names(srfGrid), sep ="")
      aesStr <- aes_string(x = nms[1], y = nms[2], z = aggregateStr,
        colour ="..level..")
      plotElems <- c(
        facet_grid(eval(substitute(formula( ~ .a),
          list(.a = as.name(names(fcts)[fcts])))), labeller = label_both),
        geom_contour(data = data, mapping = aesStr, bins = nBins),
        list(scale_colour_gradient2(name = expression(eta),
          low =  muted('darkblue'), mid = 'grey80', high = muted('darkred'))),
        geom_point(data = m$data, mapping = aes_string(x = nms[1], y = nms[2]),
          alpha = alphaRug, size =.2),
        plotElems)
      if(plotCR) {
        data$sign <- ifelse(apply(data[, 2:3], 1, prod)>0, sign(data[, 3]), 0)
        plotElems <- c(geom_tile(data = data, aes_string(x = nms[1], y = nms[2],
          fill ="sign"), colour = NA, alpha =.1),
          list(scale_fill_gradient2(name ="CI", low =  muted(alpha('darkblue', .1)),
            mid = 'white', high = muted(alpha('darkred', .1))),
            breaks = c(-1, 1), limits = c(-1.1, 1.1),  labels = c("< 0","> 0")),
          plotElems)
      }
    }
    if(length(vars)== 2 & !any(fcts)) {
      #srf:lin
      nms <-  names(xu[dfs][[1]])
      nmLin <- gsub("[aeiou]","", varsLin[1])
      linLabel <- parse(text = paste("frac(partialdiff * eta, partialdiff *",
        nmLin,")* sd(", nmLin,")"))
      aesStr <- aes_string(x = nms[1], y = nms[2], z = aggregateStr,
        colour ="..level..")
      plotElems <- c(geom_contour(data = data, mapping = aesStr, bins = nBins),
        list(scale_colour_gradient2(name = linLabel, low =  muted('darkblue'),
          mid = 'grey80', high = muted('darkred'))),
        geom_point(data = m$data, mapping = aes_string(x = nms[1], y = nms[2]),
          alpha = alphaRug, size =.5),
        plotElems)
      if(plotCR) {
        data$sign <- ifelse(apply(data[, 2:3], 1, prod)>0, sign(data[, 3]), 0)
        plotElems <- c(geom_tile(data = data, aes_string(x = nms[1], y = nms[2],
          fill ="sign"), colour = NA, alpha =.1),
          list(scale_fill_gradient2(name ="CI",
            low =  muted(alpha('darkblue', .1)), mid = 'white',
            high = muted(alpha('darkred', .1))),
            breaks = c(-1, 1), limits = c(-1.1, 1.1),  labels = c("< 0","> 0")),
          plotElems)
      }
    }
    if(length(vars)== 3 & sum(fcts)== 2) {
      #sort factors by decreasing no. of levels
      lvls <- sapply(grid[, fcts, drop = F], nlevels)
      lvlsO <- order(lvls)
      vars <- c(vars[fcts][lvlsO], vars[!fcts])

      nms <- paste(paste(names(xu)[dfs],"$", sep =""), names(srfGrid), sep ="")
      aesStr <- aes_string(x = nms[1], y = nms[2], z = aggregateStr,
        colour ="..level..")
      plotElems <- c(
        facet_grid(eval(substitute(formula(.b ~ .a),
          list(.a = as.name(names(fcts)[fcts][1]),
            .b = as.name(names(fcts)[fcts][2])))), labeller = label_both),
        geom_contour(data = data, mapping = aesStr, bins = nBins),
        list(scale_colour_gradient2(name = expression(eta),
          low =  muted('darkblue'), mid = 'grey80', high = muted('darkred'))),
        geom_point(data = m$data, mapping = aes_string(x = nms[1], y = nms[2]),
          alpha = alphaRug, size =.2),
        plotElems)
      if(plotCR) {
        data$sign <- ifelse(apply(data[, 2:3], 1, prod)>0, sign(data[, 3]), 0)
        plotElems <- c(geom_tile(data = data, aes_string(x = nms[1], y = nms[2],
          fill ="sign"), colour = NA, alpha =.1),
          list(scale_fill_gradient2(name ="CI",
            low =  muted(alpha('darkblue', .1)), mid = 'white',
            high = muted(alpha('darkred', .1)),
            breaks = c(-1, 1), limits = c(-1.1, 1.1),  labels = c("< 0","> 0"))),
          plotElems)
      }
    }
    if(length(vars)== 3 & sum(fcts)== 1) {
      nms <-  paste(paste(names(xu)[dfs],"$", sep =""), names(srfGrid), sep ="")
      nmLin <- gsub("[aeiou]","", varsLin[1])
      linLabel <- parse(text = paste("frac(partialdiff * eta, partialdiff *",
        nmLin,")"))
      aesStr <- aes_string(x = nms[1], y = nms[2], z = aggregateStr,
        colour ="..level..")
      plotElems <- c(
        facet_grid(eval(substitute(formula(. ~ .a),
          list(.a = as.name(names(fcts)[fcts])))), labeller = label_both),
        geom_contour(data = data, mapping = aesStr, bins = nBins),
        list(scale_colour_gradient2(name = linLabel,
          low =  muted('darkblue'), mid = 'grey80', high = muted('darkred'))),
        geom_point(data = m$data, mapping = aes_string(x = nms[1], y = nms[2]),
          alpha = alphaRug),
        plotElems)
      if(plotCR) {
        data$sign <- ifelse(apply(data[, 2:3], 1, prod)>0, sign(data[, 3]), 0)
        plotElems <- c(
          geom_tile(data = data, aes_string(x = nms[1], y = nms[2], fill ="sign"),
            colour = NA, alpha =.1),
          list(scale_fill_gradient2(name ="CI",
            low =  muted(alpha('darkblue', .1)), mid = 'white',
            high = muted(alpha('darkred', .1)),
            breaks = c(-1, 1), limits = c(-1.1, 1.1),  labels = c("< 0","> 0"))),
          plotElems)
      }
    }
    # end any(isSrf)
  } else {
    #sort / make grid:
    grid <- do.call(expand.grid, lapply(xu, function(x) {
      if(!is.factor(x)) {
        x <- seq(min(x), max(x), length = gridlength)
      } else x <- sort(x)
      return(x)
    }))


    pred <- data.frame(if(cumulative) {
      Reduce("+", lapply(trms, evalTerm, model = m, newdata = grid,
        aggregate = aggregate, quantiles = quantiles, returnData = FALSE))
    } else{
      evalTerm(label, model = m, newdata = grid, aggregate = aggregate,
        quantiles = quantiles, returnData = FALSE)
    })


    quantileStr <- colnames(pred)[-1]
    aggregateStr <- colnames(pred)[1]

    data <- data.frame(grid, pred)


    titleStr <- if(cumulative) paste(vars, collapse ="*") else label
    plotElems <- list(labs(title = titleStr))

    if(length(vars)== 1) {
      if(is.factor(grid[, 1])) {
        #factor
        geom <- if(plotCR) {
          geom_pointrange(aes_string(x = vars, y = aggregateStr,
            ymin = quantileStr[1], ymax = quantileStr[2]))
        } else {
          geom_point(aes_string(x = vars, y = aggregateStr))
        }
        plotElems <- c(geom, list(ylab(expression(eta))),
          plotElems)
      } else {
        #num
        aesAggrStr <- aes_string(x = vars, y = aggregateStr)
        plotElems <- c(	geom_line(data = data, aesAggrStr),
          geom_rug(data = data.frame(m$data), aes_string(x = vars, y = NULL),
            alpha = alphaRug),
          list(ylab(expression(eta))),
          plotElems)
        if(plotCR) {
          aesCRStr <- aes_string(x = vars, ymin = quantileStr[1],
            ymax = quantileStr[2])
          plotElems <- c(geom_ribbon(data = data, aesCRStr, alpha =.3), plotElems)
        }
      }
    } else {
      fcts <- sapply(grid, is.factor)
      nums <- !fcts
      if(length(vars)== 2) {
        if(any(fcts)) {
          if(all(fcts)) {
            #factor:factor
            lvls <- sapply(grid, nlevels)
            f1 <- vars[which.min(lvls)]
            f2 <- vars[-which.min(lvls)]
            geom <- if(plotCR) {
              aesStr <- aes_string(x = f2, y = aggregateStr, colour = f1,
                group = f1, ymin = quantileStr[1], ymax = quantileStr[2])
              geom_pointrange(aesStr, position = position_dodge(width =.66))
            } else {
              aesStr <- aes_string(x = f2, y = aggregateStr, colour = f1,
                group = f1)
              geom_point(aesStr, position = position_dodge(width =.66))
            }
            plotElems <- c(geom, list(ylab(expression(eta))),
              plotElems)
          } else {
            #num:factor
            aesAggrStr <- aes_string(x = vars[nums], y = aggregateStr,
              col = vars[!nums])
            plotElems <- c(geom_line(aesAggrStr),
              geom_rug(data = m$data, aes_string(x = vars[nums],
                col = paste("as.factor(", vars[!nums],")"), y = NULL),
                alpha = alphaRug),
              list(ylab(expression(eta))),
              plotElems)
            if(plotCR) {
              aesCRStr <- aes_string(x = vars[nums], ymin = quantileStr[1],
                ymax = quantileStr[2], fill = vars[!nums])
              plotElems <- c(geom_ribbon(aesCRStr, alpha =.3), plotElems)
            }
          }
        } else {
          #num:num
          meltdata <- melt(data, id.vars = vars)
          #data for geom_rug has to contain the faceting variable (and everything else too for some reason?!?)
          #so we blow up the data a little
          suppressWarnings(rugdata <- data.frame(do.call(cbind,
            c(xu, value = 0))))
          rugdata <- do.call(rbind, apply(t(levels(meltdata$variable)), 2,
            function(x) cbind(rugdata, variable = x)))

          convhull <- m$data[chull(as.matrix(m$data[, vars])), vars]
          keep <- which(insidePoly(data[, vars], convhull))

          aesStr <- aes_string(x = vars[1], y = vars[2], z = aggregateStr,
            colour ="..level..")
          plotElems <- c(geom_contour(data = data[keep, ], mapping = aesStr,
            bins = nBins),
            list(scale_colour_gradient2(name = expression(eta),
              low =  muted('darkblue'), mid = 'grey80', high = muted('darkred'),
              guide ="legend")),
            geom_rug(data = rugdata,
              mapping = aes_string(x = vars[1], y = vars[2]), alpha = alphaRug),
            plotElems)
          if(plotCR) {
            data$sign <- ifelse(apply(data[, 4:5], 1, prod)>0,
              sign(data[, 5]), 0)
            plotElems <- c(geom_tile(data = data, aes_string(x = vars[1],
              y = vars[2], fill ="sign"), colour = NA, alpha =.1),
              list(scale_fill_gradient2(name ="CI",
                low =  muted(alpha('darkblue', .1)), mid = 'white',
                high = muted(alpha('darkred', .1)), guide ="legend",
                breaks = c(-1, 1), limits = c(-1.1, 1.1),
                labels = c("< 0","> 0"))),
              plotElems)

          }
        }
      } else {
        if(!any(fcts)) { stop("You'll have to figure that one out by yourself.",
          "Plots for 3-way interactions of numeric variables are not implemented.")}

        #sort factors by decreasing no. of levels
        lvls <- sapply(grid[, fcts, drop = F], nlevels)
        lvlsO <- order(lvls)
        vars <- c(vars[fcts][lvlsO], vars[nums])
        if(sum(fcts)== 1) {
          #fct:num:num
          meltdata <- melt(data, id.vars = vars)
          #data for geom_rug has to contain the faceting variable (and everything else too for some reason?!?)
          #so we blow up the data a little
          #suppressWarnings(rugdata <- data.frame(do.call(cbind, c(xu, value = 0))))
          #rugdata <- do.call(rbind, apply(t(levels(meltdata$variable)), 2, function(x) cbind(rugdata, variable = x)))

          convhull <- m$data[chull(as.matrix(m$data[, vars[2:3]])), vars[2:3]]
          keep <- which(insidePoly(data[, vars[2:3]], convhull))

          aesStr <- aes_string(x = vars[2], y = vars[3], z = aggregateStr,
            colour ="..level..")
          plotElems <- c(geom_contour(data = data[keep, ], mapping = aesStr,
            bins = nBins),
            list(scale_colour_gradient2(name = expression(eta),
              low =  muted('darkblue'), mid = 'grey80', high = muted('darkred'))),
            plotElems)

          if(plotCR) {
            data$sign <- ifelse(apply(data[, 5:6], 1, prod)>0,
              sign(data[, 5]), 0)
            plotElems <- c(geom_tile(data = data,
              aes_string(x = vars[2], y = vars[3], fill ="sign"), colour = NA,
              alpha =.1),
              list(scale_fill_gradient2(name ="CI",
                low =  alpha(muted('darkblue'), .1), mid = alpha('white', .1),
                high = alpha(muted('darkred'), .1), breaks = c(-1, 1),
                limits = c(-1.1, 1.1),  labels = c("< 0","> 0"))),
              plotElems)
          }
          plotElems <- c(list(facet_grid(eval(substitute(formula( ~ .a),
            list(.a = as.name(vars[1])))), labeller = label_both)), plotElems)



        }
        if(sum(fcts)== 2) {
          #fct:fct:num
          aesAggrStr <- aes_string(x = vars[3], y = aggregateStr, col = vars[2])
          plotElems <- c(geom_line(aesAggrStr),
            geom_rug(data = m$data, aes_string(x = vars[3],
              col = paste("as.factor(", vars[2],")"), y = NULL),
              alpha = alphaRug),
            plotElems)
          if(plotCR) {
            aesCRStr <- aes_string(x = vars[3], ymin = quantileStr[1],
              ymax = quantileStr[2], fill = vars[2])
            plotElems <- c(geom_ribbon(aesCRStr, alpha =.3), plotElems)
          }
          facetFrml <- eval(substitute(formula( ~ .a),
            list(.a = as.name(vars[1]))))
          plotElems <- c(list(facet_grid(facetFrml, labeller = label_both),
            list(ylab(expression(eta)))), plotElems)
        }
        if(all(fcts)) {
          #fct:fct:fct
          geom <- if(plotCR) {
            aesStr <- aes_string(x = vars[3], y = aggregateStr,
              colour = vars[2], group = vars[2], ymin = quantileStr[1],
              ymax = quantileStr[2])
            geom_pointrange(aesStr, position = position_dodge(width =.66))
          } else {
            aesStr <- aes_string(x = vars[3], y = aggregateStr, colour = vars[2],
              group = vars[2])
            geom_point(aesStr, position = position_dodge(width =.66))
          }
          facetFrml <- eval(substitute(formula( ~ .a),
            list(.a = as.name(vars[1]))))
          plotElems <- c(list(facet_grid(facetFrml, labeller = label_both),
            geom, list(ylab(expression(eta)))), plotElems)
        }
      }
    }
  }
  ret <- ggplot(data = data)
  return(p <- ret + plotElems + ggElems)
}

#'Generates graphical summaries of a fitted model
#'
#'This function plots the estimated linear predictors of the terms in a model on
#'a grid of values. By default displays all 3-way, 2-way interactions and main
#'effects present in the model. Starting with ggplot-0.9.2 these are no longer
#'aligned by their axes due to internal changes in grid and ggplot2. Uses
#'\pkg{gridExtra}'s  \code{\link[gridExtra]{marrangeGrob}} to arrange the plots
#'for the terms, also over multiple pages if necessary. This means the graphics
#'device type is temporarily set to the value of \code{interactive.dev} in
#'interactive use in RStudio if necessary since the \code{RStudioGD} graphical
#'device does not support opening multiple pages.
#'
#'@param x a fitted \code{spikeSlabGAM} model
#'@param labels a character vector of names of model terms to be plotted
#'@param cumulative Defaults to TRUE, in which case all lower order terms that
#'  are involved in an interaction are cumulated and then plotted (e.g, if a
#'  model contains 2 smooth effects and their interaction, ONLY the sum of the
#'  marginal smooth and linear terms and all their interactions are plotted.) If
#'  FALSE, a separate plot for every term in the model is produced.
#'@param commonEtaScale use the same limits for all vertical axes of the
#'  different panels? Defaults to FALSE. Can be useful to compare effect sizes
#'  more easily between panels, but tends to mess up the scales.
#'@param aggregate (function) which summary statistic to use for the posterior
#'  of the model terms. Defaults to the mean.
#'@param quantiles which quantiles to use for the borders of credible regions.
#'  Defaults to 10\% and 90\% percentiles. Set to NULL to omit plotting credible
#'  regions.
#'@param gridlength length of the (univariate) grids for the covariates on which
#'  to evaluate the posterior. Defaults to 20.
#'@param base_size default base font size for plot (see e.g.
#'  \code{\link[ggplot2]{theme_gray}})
#'@param ggElems a list of plot elements to give to \code{ggplot}. Use this to
#'  supply custom themes or colors, fortify the plot(s) with partial residuals
#'  etc. Defaults to an empty list. Unless specified differently here, the
#'  default ggplot-theme (\code{\link[ggplot2]{theme_gray}}) is changed to a
#'  white background with major gridlines in gray (\code{'grey95'}), no minor
#'  grid lines, and smaller text for the legends.
#'@param nrow	number of rows per page, defaults to min(sqrt(no. of plots), 3).
#'  See \code{\link[gridExtra]{marrangeGrob}}.
#'@param ncol  number of columns per page, defaults to min((no. of plots)/nrow,
#'  3). See \code{\link[gridExtra]{marrangeGrob}}.
#'@param interactive.dev alternative device to use in interactive mode in
#'  RStudio if output needs to be spread on multiple pages, since the RStudio
#'  graphical device does not support opening multiple displays.
#'@param ... 	arguments passed to \code{\link[gridExtra]{marrangeGrob}}.
#'
#'@seealso \code{\link{plotTerm}} for more details on the specific plots
#'@note Note that \code{cumulative = TRUE} will only find all relevant terms to
#'  accumulate if, for all numeric covariates that have a smooth term, the
#'  smooth term is specified \emph{after} the linear term in the formula.
#'@export
#'@importFrom gridExtra marrangeGrob
#'@return a list of \code{\link[ggplot2]{ggplot}}-objects (invisible)
#'@author Fabian Scheipl
#' @examples
#' #see ?spikeSlabGAM
#FIXME: only wrks if smooth trms are specified AFTER linear terms...
plot.spikeSlabGAM <- function(x, labels = NULL, cumulative = TRUE,
  commonEtaScale = FALSE, aggregate = mean, quantiles = c(.1, .9),
  gridlength = 20, base_size = 12, ggElems = list(),
  nrow = min(ceiling(sqrt(length(plotList))), 3),
  ncol = min(ceiling(length(plotList)/nrow), 3),
  interactive.dev = c("x11", "quartz", "windows"), ...) {

  interactive.dev <- match.arg(interactive.dev)

  # add custom theme / overwrite default theme
  ggElems <- c(list(theme(panel.grid.minor = element_line(colour = NA),
    panel.grid.major = element_line(colour = "grey95"),
    panel.background = element_rect(fill ="white"),
    legend.key = element_rect(fill = "white", colour = "white"),
    legend.text = element_text(size = base_size * 0.6),
    legend.title = element_text(face = "bold", size = base_size * 0.6,
      hjust = 0))),
    ggElems)


  if(is.null(labels)) {
    if(cumulative) {
      #find out which interactions can subsume all lower order effects...
      trms <- terms(x$formula)
      terms <- names(x$predvars)
      order <- attr(trms,"order")
      if(max(order)>3) {
        warning("You'll have to figure that one out by yourself.",
          "Plots for 4-way interactions and higher are not implemented.")
        rm <- which(order>3)
        terms <- terms[-rm]
        order <- order[-rm]
      }


      labels <- terms[order == 1]
      varsList <- sapply(terms,
        function(x) {
          sapply(
            unlist(strsplit(gsub("(:?[[:lower:]]+)(\\([[:alnum:]._]+\\))","\\2",
              x), c("\\)"))),
            function(y) gsub("(", "", y, fixed = T))
        })
      names(varsList) <- terms
      vars <- unique(unlist(varsList))
      #FIXME: this only works if smooth terms are specified in the formula
      # AFTER linear terms...
      combos <-  rev(names(varsList))[!duplicated(rev(varsList))]
      comboOrder <- rev(order)[!duplicated(rev(varsList))]




      plotLabels <- c()
      for(o in max(order):1) {
        plotLabels <- c(plotLabels, combos[comboOrder == o])
        if(any(comboOrder == o)) {
          # find main effect(s)
          done <- sapply(combos[comboOrder == o], allInvolvedTerms, m = x)
          #rm all term combinations that are already involved in the current one
          rm <- unlist(sapply(done, match, table = combos))
          rm <- sort(unique(rm[!is.na(rm)]))
          combos <- combos[-rm]
          comboOrder <- comboOrder[-rm]
        }
      }
    } else {
      plotLabels <- names(x$predvars)
    }
  } else{
    plotLabels <- labels
  }

  #TODO: order plots as in formula or via term importance...?
  if(is.null(labels)) plotLabels <- rev(plotLabels)


  plotList <- lapply(as.list(plotLabels), plotTerm, m = x, cumulative = cumulative,
    aggregate = aggregate, quantiles = quantiles, gridlength = gridlength,
    ggElems = ggElems)
  nPlots <- length(plotList)


  if(commonEtaScale) {
    #get range of eta across all plots...
    rangeEta <- range(lapply(plotList, function(p) {
      useCols <- unlist(mapply(grep,
        pattern = c("^eta$", "percentile$"),
        MoreArgs = list(x = colnames(p$data))))
      range(p$data[, useCols])
    }))
    #.. and apply to all 2-d plots
    plotList <- lapply(plotList, function(p) {
      geoms <- sapply(1:length(p$layers), function(i) {
        p$layers[[i]]$geom$objname
      })
      if("contour" %in% geoms) {
        return(p)
      }
      else{
        return(p + coord_cartesian(ylim = rangeEta))
      }
    })
  }
  if(dev.interactive()) {
    # work around RStudio quirk,
    # see https://github.com/baptiste/gridextra/issues/6
    rstudio <- options()$device == "RStudioGD"
    if(rstudio & nrow * ncol < length(plotList)) {
      message("Temporarily setting <options(device ='", interactive.dev,
        "')> to open multiple graphics displays since dispatch on multiple graphic devices is not supported in RStudio.")
      options(device = interactive.dev)
      on.exit(options(device ="RStudioGD"))
    }
  }
  print(marrangeGrob(plotList, nrow = nrow, ncol = ncol, ...))
  invisible(plotList)
}
