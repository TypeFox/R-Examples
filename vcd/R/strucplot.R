################################################################
### strucplot - generic plot framework for mosaic-like layouts
### 2 core functions are provided: struc_mosaic and struc_assoc
################################################################

strucplot <- function(## main parameters
                      x,
                      residuals = NULL,
                      expected = NULL,
		      condvars = NULL,
                      shade = NULL,
                      type = c("observed", "expected"),
                      residuals_type = NULL,
                      df = NULL,

                      ## layout
                      split_vertical = NULL,
                      spacing = spacing_equal,
                      spacing_args = list(),
                      gp = NULL,
		      gp_args = list(),
                      labeling = labeling_border,
                      labeling_args = list(),
                      core = struc_mosaic,
                      core_args = list(),
                      legend = NULL,
                      legend_args = list(),

                      main = NULL,
                      sub = NULL,
                      margins = unit(3, "lines"),
                      title_margins = NULL,
                      legend_width = NULL,

                      ## control parameters
                      main_gp = gpar(fontsize = 20),
                      sub_gp = gpar(fontsize = 15),
                      newpage = TRUE,
                      pop = TRUE,
                      return_grob = FALSE,
                      keep_aspect_ratio = NULL,
                      prefix = "",
                      ...
                      ) {
  ## default behaviour of shade
  if (is.null(shade)) shade <- !is.null(gp) || !is.null(expected)

  type <- match.arg(type)
  if (is.null(residuals)) {
      residuals_type <- if (is.null(residuals_type))
          "pearson"
       else
           match.arg(tolower(residuals_type),
                     c("pearson", "deviance", "ft"))
  } else {
      if (is.null(residuals_type))
          residuals_type <- ""
  }

  ## convert structable object
  if (is.structable(x)) {
    if (is.null(split_vertical))
      split_vertical <- attr(x, "split_vertical")
    x <- as.table(x)
  }
  if (is.null(split_vertical))
    split_vertical <- FALSE

  ## table characteristics
  d <- dim(x)
  dl <- length(d)
  dn <- dimnames(x)
  if (is.null(dn))
    dn <- dimnames(x) <- lapply(d, seq)
  dnn <- names(dimnames(x))
  if (is.null(dnn))
    dnn <- names(dn) <- names(dimnames(x)) <- LETTERS[1:dl]

  ## replace NAs by 0
  if (any(nas <- is.na(x))) x[nas] <- 0

  ## model fitting:
  ## calculate df and expected if needed
  ## (used for inference in some shading (generating) functions).
  ## note: will *not* be calculated if residuals are given
  if ((is.null(expected) && is.null(residuals)) ||
      !is.numeric(expected)) {
      if (!is.null(df))
          warning("Using calculated degrees of freedom.")
      if (inherits(expected, "formula")) {
          fm <- loglm(expected, x, fitted = TRUE)
          expected <- fitted(fm)
          df <- fm$df
      } else {
          if (is.null(expected))
              expected <- if (is.null(condvars))
                  as.list(1:dl)
              else
                  lapply((condvars + 1):dl, c, seq(condvars))

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

  ## splitting
  if (length(split_vertical) == 1)
    split_vertical <- rep(c(split_vertical, !split_vertical), length.out = dl)

  if (is.null(keep_aspect_ratio))
    keep_aspect_ratio <- dl < 3

  ## spacing
  if (is.function(spacing)) {
    if (inherits(spacing, "grapcon_generator"))
      spacing <- do.call("spacing", spacing_args)
    spacing <- spacing(d, condvars)
  }

  ## gp (color, fill, lty, etc.) argument
  if (shade) {
    if (is.null(gp)) gp <- shading_hcl
    if (is.function(gp)) {
      if (is.null(legend) || (is.logical(legend) && legend))
        legend <- legend_resbased
      gpfun <- if (inherits(gp, "grapcon_generator"))
        do.call("gp", c(list(x, residuals, expected, df), as.list(gp_args))) else gp
      gp <- gpfun(residuals)
    } else if (!is.null(legend) && !(is.logical(legend) && !legend))
      stop("gp argument must be a shading function for drawing a legend")
  } else {
    if(!is.null(gp)) {
      warning("gp parameter ignored since shade = FALSE")
      gp <- NULL
    }
  }

  ## choose gray when no shading is used
  if (is.null(gp)) gp <- gpar(fill = grey(0.8))

  ## recycle gpar values in the last dimension
  size <- prod(d)
  FUN <- function(par) {
    if (is.structable(par))
      par <- as.table(par)
    if (length(par) < size || is.null(dim(par))) aperm(array(par, dim = rev(d))) else par
  }
  gp <- structure(lapply(gp, FUN), class = "gpar")

  ## set up page
  if (newpage)
      grid.newpage()
  if (keep_aspect_ratio)
      pushViewport(viewport(width = 1, height = 1, default.units = "snpc"))

  pushViewport(vcdViewport(mar = margins,
                           oma = title_margins,
                           legend = shade && !(is.null(legend) || is.logical(legend) && !legend),
                           main = !is.null(main), sub = !is.null(sub),
                           keep_aspect_ratio = keep_aspect_ratio,
                           legend_width = legend_width,
                           prefix = prefix))

  ## legend
  if (inherits(legend, "grapcon_generator"))
    legend <- do.call("legend", legend_args)
  if (shade && !is.null(legend) && !(is.logical(legend) && !legend)) {
    seekViewport(paste(prefix, "legend", sep = ""))
    residuals_type <- switch(residuals_type,
                             deviance = "deviance\nresiduals:",
                             ft = "Freeman-Tukey\nresiduals:",
                             pearson = "Pearson\nresiduals:",
                             residuals_type)
    legend(residuals, gpfun, residuals_type)
  }

  ## titles
  if (!is.null(main)) {
    seekViewport(paste(prefix, "main", sep = ""))
    if (is.logical(main) && main)
      main <- deparse(substitute(x))
    grid.text(main, gp = main_gp)
  }

  if (!is.null(sub)) {
    seekViewport(paste(prefix, "sub", sep = ""))
    if (is.logical(sub) && sub && is.null(main))
      sub <- deparse(substitute(x))
    grid.text(sub, gp = sub_gp)
  }

  ## make plot
  seekViewport(paste(prefix, "plot", sep = ""))

  if (inherits(core, "grapcon_generator"))
    core <- do.call("core", core_args)
  core(residuals = residuals,
       observed = if (type == "observed") x else expected,
       expected = if (type == "observed") expected else x,
       spacing = spacing,
       gp = gp,
       split_vertical = split_vertical,
       prefix = prefix)

  upViewport(1)

  ## labels
  if (is.logical(labeling))
    labeling <- if (labeling) labeling_border else NULL
  if (!is.null(labeling)) {
    if (inherits(labeling, "grapcon_generator"))
      labeling <- do.call("labeling", c(labeling_args, list(...)))
    labeling(dn, split_vertical, condvars, prefix)
  }

  ## pop/move up viewport

  seekViewport(paste(prefix, "base", sep = ""))
  ## one more up if sandwich-mode
  if (pop)
      popViewport(1 + keep_aspect_ratio)
  else
      upViewport(1 + keep_aspect_ratio)

  ## return visualized table
  if (return_grob)
      invisible(structure(structable(if (type == "observed") x else expected,
                                     split_vertical = split_vertical),
                          grob = grid.grab()
                          )
                )
  else
      invisible(structable(if (type == "observed") x else expected,
                           split_vertical = split_vertical))
}

vcdViewport <- function(mar = rep.int(2.5, 4),
                        legend_width = unit(5, "lines"),
                        oma = NULL,
                        legend = FALSE, main = FALSE, sub = FALSE,
                        keep_aspect_ratio = TRUE,
                        prefix = "")
{
  ## process parameters
  if (is.null(legend_width))
    legend_width <- unit(5 * legend, "lines")
  if (!is.unit(legend_width))
    legend_width <- unit(legend_width, "lines")

  if (legend && !main && !sub && keep_aspect_ratio) main <- sub <- TRUE
  mar <- if (!is.unit(mar))
    unit(pexpand(mar, 4, rep.int(2.5, 4), c("top","right","bottom","left")), "lines")
  else
    rep(mar, length.out = 4)
  if (is.null(oma)) {
    space <- if (legend && keep_aspect_ratio)
      legend_width + mar[2] + mar[4] - mar[1] - mar[3]
    else unit(0, "lines")
    oma <- if (main && sub)
      max(unit(2, "lines"), 0.5 * space)
    else if (main)
      unit.c(max(unit(2, "lines"), space), unit(0, "lines"))
    else if (sub)
      unit.c(unit(0, "lines"), max(unit(2, "lines"), space))
    else
      0.5 * space
  }
  oma <- if (!is.unit(oma))
    unit(pexpand(oma, 2, rep.int(2, 2), c("top","bottom")), "lines")
  else
    rep(oma, length.out = 2)

  ## set up viewports
  vpPlot <- vpStack(viewport(layout.pos.col = 2, layout.pos.row = 3),
                    viewport(width = 1, height = 1, name = paste(prefix, "plot", sep = ""),
                             default.units = if (keep_aspect_ratio) "snpc" else "npc"))
  vpMarginBottom <- viewport(layout.pos.col = 2, layout.pos.row = 4,
                             name = paste(prefix, "margin_bottom", sep = ""))
  vpMarginLeft <- viewport(layout.pos.col = 1, layout.pos.row = 3,
                           name = paste(prefix, "margin_left", sep = ""))
  vpMarginTop <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                          name = paste(prefix, "margin_top", sep = ""))
  vpMarginRight <- viewport(layout.pos.col = 3, layout.pos.row = 3,
                            name = paste(prefix, "margin_right", sep = ""))
  vpCornerTL <- viewport(layout.pos.col = 1, layout.pos.row = 2,
                         name = paste(prefix, "corner_top_left", sep = ""))
  vpCornerTR <- viewport(layout.pos.col = 3, layout.pos.row = 2,
                         name = paste(prefix, "corner_top_right", sep = ""))
  vpCornerBL <- viewport(layout.pos.col = 1, layout.pos.row = 4,
                         name = paste(prefix, "corner_bottom_left", sep = ""))
  vpCornerBR <- viewport(layout.pos.col = 3, layout.pos.row = 4,
                         name = paste(prefix, "corner_bottom_right", sep = ""))

  vpLegend <- viewport(layout.pos.col = 4, layout.pos.row = 3,
                       name = paste(prefix, "legend", sep = ""))
  vpLegendTop <- viewport(layout.pos.col = 4, layout.pos.row = 2,
                          name = paste(prefix, "legend_top", sep = ""))
  vpLegendSub <- viewport(layout.pos.col = 4, layout.pos.row = 4,
                          name = paste(prefix, "legend_sub", sep = ""))
  vpBase <- viewport(layout = grid.layout(5, 4,
                       widths = unit.c(mar[4], unit(1, "null"), mar[2], legend_width),
                       heights = unit.c(oma[1], mar[1], unit(1, "null"), mar[3], oma[2])),
                     name = paste(prefix, "base", sep = ""))
  vpMain <- viewport(layout.pos.col = 1:4, layout.pos.row = 1,
                     name = paste(prefix, "main", sep = ""))
  vpSub <- viewport(layout.pos.col = 1:4, layout.pos.row = 5,
                    name = paste(prefix, "sub", sep = ""))

  vpTree(vpBase, vpList(vpMain, vpMarginBottom, vpMarginLeft, vpMarginTop,
                        vpMarginRight, vpLegendTop, vpLegend,
                        vpLegendSub, vpCornerTL, vpCornerTR,
                        vpCornerBL, vpCornerBR, vpPlot, vpSub))
}

