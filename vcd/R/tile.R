tile <- function(x, ...)
  UseMethod("tile")

tile.formula <-
function(formula, data = NULL,
         ..., main = NULL, sub = NULL, subset = NULL, na.action = NULL)
{
  if (is.logical(main) && main)
    main <- deparse(substitute(data))
  else if (is.logical(sub) && sub)
    sub <- deparse(substitute(data))

  m <- match.call(expand.dots = FALSE)
  edata <- eval(m$data, parent.frame())

  fstr <- strsplit(paste(deparse(formula), collapse = ""), "~")
  vars <- strsplit(strsplit(gsub(" ", "", fstr[[1]][2]), "\\|")[[1]], "\\+")
  varnames <- vars[[1]]

  dep <- gsub(" ", "", fstr[[1]][1])
  if (!dep %in% c("","Freq")) {
     if (all(varnames == ".")) {
       varnames <- if (is.data.frame(data))
         colnames(data)
       else
         names(dimnames(as.table(data)))
       varnames <- varnames[-which(varnames %in% dep)]
     }

    varnames <- c(varnames, dep)
  }


  if (inherits(edata, "ftable") || inherits(edata, "table") || length(dim(edata)) > 2) {
    dat <- as.table(data)
    if(all(varnames != ".")) {
      ind <- match(varnames, names(dimnames(dat)))
      if (any(is.na(ind)))
        stop(paste("Can't find", paste(varnames[is.na(ind)], collapse=" / "), "in", deparse(substitute(data))))

      dat <- margin.table(dat, ind)
    }
    tile.default(dat, main = main, sub = sub, ...)
  } else {
      m <- m[c(1, match(c("formula", "data", "subset", "na.action"), names(m), 0))]
      m[[1]] <- as.name("xtabs")
      m$formula <-
          formula(paste(if("Freq" %in% colnames(data)) "Freq",
                        "~",
                        paste(varnames, collapse = "+")))
      tab <- eval(m, parent.frame())
      tile.default(tab, main = main, sub = sub, ...)
  }
}

tile.default <-
function(x,
         tile_type = c("area", "squaredarea", "height", "width"),
         halign = c("left", "center", "right"),
         valign = c("bottom", "center", "top"),
         split_vertical = NULL,
         shade = FALSE,
         spacing = spacing_equal(unit(1, "lines")),
         set_labels = NULL,
         margins = unit(3, "lines"),
         keep_aspect_ratio = FALSE,
         legend = NULL,
         legend_width = NULL,
         squared_tiles = TRUE,
         main = NULL, sub = NULL, ...)
{
    ## argument handling
    if (is.logical(main) && main)
        main <- deparse(substitute(x))
    else if (is.logical(sub) && sub)
        sub <- deparse(substitute(x))

    tile_type <- match.arg(tile_type)
    halign <- match.arg(halign)
    valign <- match.arg(valign)

    x <- as.table(x)
    dl <- length(d <- dim(x))

    ## determine starting positions
    xpos <- 1 - (halign == "left") - 0.5 * (halign == "center")
    ypos <- 1 - (valign == "bottom") - 0.5 * (valign == "center")

    ## heuristic to adjust right/bottom margin to obtain squared tiles
    ## FIXME: better push another viewport?
    if (squared_tiles) {
        ## splitting argument
        if (is.structable(x) && is.null(split_vertical))
            split_vertical <- attr(x, "split_vertical")
        if (is.null(split_vertical))
            split_vertical <- FALSE
        if (length(split_vertical) == 1)
            split_vertical <- rep(c(split_vertical, !split_vertical),
                                  length.out = dl)
        if (length(split_vertical) < dl)
            split_vertical <- rep(split_vertical, length.out = dl)

        ## compute resulting dimnension
        dflat <- dim(unclass(structable(x, split_vertical = split_vertical)))

        ## adjust margins
        spacing <- spacing(d)
        delta <- abs(dflat[1] - dflat[2])
        fac <- delta / max(dflat)
        un <-  unit(fac, "npc") - unit(fac * 5 / spacing[[1]][[1]], "lines")
        leg <- if (shade) {
            if (is.null(legend_width))
                unit(5, "lines")
            else legend_width
        } else unit(0, "npc")
        if (dflat[1] < dflat[2])
            margins <- margins + unit.c(unit(0, "npc"), unit(0, "npc"),
                                        un + leg, unit(0, "npc"))
        if (dflat[1] > dflat[2])
            margins <- margins + unit.c(unit(0, "npc"), un - leg,
                                        unit(0, "npc"), unit(0, "npc"))
        if (dflat[1] == dflat[2])
            margins <- margins + unit.c(unit(0, "npc"), unit(0, "npc"),
                                        leg, unit(0, "npc"))
    }

    ## create dummy labels if some are duplicated
    ## and set the labels via set_labels
    dn <- dimnames(x)
    if (any(unlist(lapply(dn, duplicated)))) {
        dimnames(x) <- lapply(dn, seq_along)
        if (is.null(set_labels))
            set_labels <- lapply(dn, function(i) structure(i, names = seq(i)))
    }

    ## workhorse function creating bars
    panelfun <- function(residuals, observed, expected, index, gp, name) {
        xprop <- expected / max(expected)
        if (tile_type == "height")
            grid.rect(x = xpos, y = ypos,
                      height = xprop[t(index)], width = 1,
                      gp = gp, just = c(halign, valign), name = name)
        else if (tile_type == "width")
            grid.rect(x = xpos, y = ypos,
                      width = xprop[t(index)], height = 1,
                      gp = gp, just = c(halign, valign), name = name)
        else if (tile_type == "area")
            grid.rect(x = xpos, y = ypos,
                      width = sqrt(xprop[t(index)]),
                      height = sqrt(xprop[t(index)]),
                      gp = gp, just = c(halign, valign), name = name)
        else
            grid.rect(x = xpos, y = ypos,
                      width = xprop[t(index)], height = xprop[t(index)],
                      gp = gp, just = c(halign, valign), name = name)

    }

    mycore <- function(residuals, observed, expected = NULL,
                       spacing, gp, split_vertical, prefix = "") {
        struc_mosaic(panel = panelfun)(residuals,
                     array(1, dim = d, dimnames = dimnames(observed)),
                     expected = observed,
                     spacing, gp, split_vertical, prefix)
    }

    strucplot(x,
              core = mycore,
              spacing = spacing,
              keep_aspect_ratio = keep_aspect_ratio,
              margins = margins,
              shade = shade,
              legend = legend,
              legend_width = legend_width,
              main = main,
              sub = sub,
              set_labels = set_labels,
              ...)
}
