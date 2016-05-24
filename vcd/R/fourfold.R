## Modifications - MF - 1 Dec 2010
#  -- change default colors to more distinguishable values
#  -- allow to work with >3 dimensional arrays
#  -- modified defaults for mfrow/mfcol to give landscape display, nr <= nc, rather than nr >= nc

# Take a 2+D array and return a 3D array, with dimensions 3+ as a single dimension
# Include as a separate function, since it is useful in other contexts
array3d <- function(x, sep=':') {
	if(length(dim(x)) == 2) {
        x <- if(is.null(dimnames(x)))
            array(x, c(dim(x), 1))
        else
            array(x, c(dim(x), 1), c(dimnames(x), list(NULL)))
        return(x)
    }
  else if(length(dim(x))==3) return(x)
  else {
	  x3d <- array(x, c(dim(x)[1:2], prod(dim(x)[-(1:2)])))
	  if (!is.null(dimnames(x))) {
		  n3d <- paste(names(dimnames(x))[-(1:2)], collapse=sep)
		  d3d <- apply(expand.grid(dimnames(x)[-(1:2)]), 1, paste, collapse=sep)
		  dimnames(x3d) <- c(dimnames(x)[1:2], list(d3d))
		  names(dimnames(x3d))[3] <- n3d
  	}
  	return(x3d)
  }
}

"fourfold" <-
function(x,
#        color = c("#99CCFF","#6699CC","#FF5050","#6060A0", "#FF0000", "#000080"),
         color = c("#99CCFF","#6699CC","#FFA0A0","#A0A0FF", "#FF0000", "#000080"),
         conf_level = 0.95,
         std = c("margins", "ind.max", "all.max"), margin = c(1, 2),
         space = 0.2, main = NULL, sub = NULL, mfrow = NULL, mfcol = NULL, extended = TRUE,
         ticks = 0.15, p_adjust_method = p.adjust.methods, newpage = TRUE,
         fontsize = 12, default_prefix = c("Row", "Col", "Strata"), sep = ": ", varnames = TRUE,
         return_grob = FALSE)

{
    ## Code for producing fourfold displays.
    ## Reference:
    ##   Friendly, M. (1994).
    ##   A fourfold display for 2 by 2 by \eqn{k} tables.
    ##   Technical Report 217, York University, Psychology Department.
    ##   http://datavis.ca/papers/4fold/4fold.pdf
    ##
    ## Implementation notes:
    ##
    ##   We need plots with aspect ratio FIXED to 1 and glued together.
    ##   Hence, even if k > 1 we prefer keeping everything in one plot
    ##   region rather than using a multiple figure layout.
    ##   Each 2 by 2 pie is is drawn into a square with x/y coordinates
    ##   between -1 and 1, with row and column labels in [-1-space, -1]
    ##   and [1, 1+space], respectively.  If k > 1, strata labels are in
    ##   an area with y coordinates in [1+space, 1+(1+gamma)*space],
    ##   where currently gamma=1.25.  The pies are arranged in an nr by
    ##   nc layout, with horizontal and vertical distances between them
    ##   set to space.
    ##
    ##   The drawing code first computes the complete are of the form
    ##     [0, totalWidth] x [0, totalHeight]
    ##   needed and sets the world coordinates using plot.window().
    ##   Then, the strata are looped over, and the corresponding pies
    ##   added by filling rows or columns of the layout as specified by
    ##   the mfrow or mfcol arguments.  The world coordinates are reset
    ##   in each step by shifting the origin so that we can always plot
    ##   as detailed above.

    if(!is.array(x))
        stop("x must be an array")
    dimx <- dim(x)   # save original dimensions for setting default mfrow/mfcol when length(dim(x))>3
    x <- array3d(x)
    if(any(dim(x)[1:2] != 2))
        stop("table for each stratum must be 2 by 2")
    dnx <- dimnames(x)
    if(is.null(dnx))
        dnx <- vector("list", 3)
    for(i in which(sapply(dnx, is.null)))
        dnx[[i]] <- LETTERS[seq(length = dim(x)[i])]
    if(is.null(names(dnx)))
        i <- 1 : 3
    else
        i <- which(is.null(names(dnx)))
    if(any(i > 0))
        names(dnx)[i] <- default_prefix[i]
    dimnames(x) <- dnx
    k <- dim(x)[3]

    if(!((length(conf_level) == 1) && is.finite(conf_level) &&
         (conf_level >= 0) && (conf_level < 1)))
        stop("conf_level must be a single number between 0 and 1")
    if(conf_level == 0)
        conf_level <- FALSE

    std <- match.arg(std)

    findTableWithOAM <- function(or, tab) {
        ## Find a 2x2 table with given odds ratio `or' and the margins
        ## of a given 2x2 table `tab'.
        m <- rowSums(tab)[1]
        n <- rowSums(tab)[2]
        t <- colSums(tab)[1]
        if(or == 1)
            x <- t * n / (m + n)
        else if(or == Inf)
            x <- max(0, t - m)
        else {
            A <- or - 1
            B <- or * (m - t) + (n + t)
            C <- - t * n
            x <- (- B + sqrt(B ^ 2 - 4 * A * C)) / (2 * A)
        }
        matrix(c(t - x, x, m - t + x, n - x), nrow = 2)
    }

    drawPie <- function(r, from, to, n = 500, color = "transparent") {
        p <- 2 * pi * seq(from, to, length = n) / 360
        x <- c(cos(p), 0) * r
        y <- c(sin(p), 0) * r
        grid.polygon(x, y, gp = gpar(fill = color), default.units = "native")
        invisible(NULL)
    }

    stdize <- function(tab, std, x) {
        ## Standardize the 2 x 2 table `tab'.
        if(std == "margins") {
            if(all(sort(margin) == c(1, 2))) {
                ## standardize to equal row and col margins
                u <- sqrt(odds(tab)$or)
                u <- u / (1 + u)
                y <- matrix(c(u, 1 - u, 1 - u, u), nrow = 2)
            }
            else if(margin %in% c(1, 2))
                y <- prop.table(tab, margin)
            else
                stop("incorrect margin specification")
        }
        else if(std == "ind.max")
            y <- tab / max(tab)
        else if(std == "all.max")
            y <- tab / max(x)
        y
    }

    odds <- function(x) {
        ## Given a 2 x 2 or 2 x 2 x k table `x', return a list with
        ## components `or' and `se' giving the odds ratios and standard
        ## deviations of the log odds ratios.
        if(length(dim(x)) == 2) {
            dim(x) <- c(dim(x), 1)
            k <- 1
        }
        else
            k <- dim(x)[3]
        or <- double(k)
        se <- double(k)
        for(i in 1 : k) {
            f <- x[ , , i]
            if(any(f == 0))
                f <- f + 0.5
            or[i] <- (f[1, 1] * f[2, 2]) / (f[1, 2] * f[2, 1])
            se[i] <- sqrt(sum(1 / f))
        }
        list(or = or, se = se)
    }

    gamma <- 1.25                       # Scale factor for strata labels

    angle.f <- c( 90, 180,  0, 270)     # `f' for `from'
    angle.t <- c(180, 270, 90, 360)     # `t' for `to'

    byrow <- FALSE
    if(!is.null(mfrow)) {
        nr <- mfrow[1]
        nc <- mfrow[2]
    }
    else if(!is.null(mfcol)) {
        nr <- mfcol[1]
        nc <- mfcol[2]
        byrow <- TRUE
    }
    else if(length(dimx)>3) {
        nr <- dimx[3]
        nc <- prod(dimx[-(1:3)])
    }
    else {
#       nr <- ceiling(sqrt(k))
        nr <- round(sqrt(k))
        nc <- ceiling(k / nr)
    }
    if(nr * nc < k)
        stop("incorrect geometry specification")
    if(byrow)
        indexMatrix <- expand.grid(1 : nc, 1 : nr)[, c(2, 1)]
    else
        indexMatrix <- expand.grid(1 : nr, 1 : nc)

    totalWidth <- nc * 2 * (1 + space) + (nc - 1) * space
    totalHeight <- if(k == 1)
        2 * (1 + space)
    else
        nr * (2 + (2 + gamma) * space) + (nr - 1) * space
    xlim <- c(0, totalWidth)
    ylim <- c(0, totalHeight)
    if (newpage) grid.newpage()
    if (!is.null(main) || !is.null(sub))
      pushViewport(viewport(height = 1 - 0.1 * sum(!is.null(main), !is.null(sub)),
                            width = 0.9,
                            y = 0.5 - 0.05 * sum(!is.null(main), - !is.null(sub))
                            )
                   )

    pushViewport(viewport(xscale = xlim, yscale = ylim,
                           width = unit(min(totalWidth / totalHeight, 1), "snpc"),
                           height = unit(min(totalHeight / totalWidth, 1), "snpc")))
    o <- odds(x)

    ## perform logoddsratio-test for each stratum (H0: lor = 0) and adjust p-values
    if(is.numeric(conf_level) && extended)
      p.lor.test <- p.adjust(sapply(1 : k, function(i) {
                               u <- abs(log(o$or[i])) / o$se[i]
                               2 * (1 - pnorm(u))
                             }),
                             method = p_adjust_method
                             )

    scale <- space / (2 * convertY(unit(1, "strheight", "Ag"), "native", valueOnly = TRUE) )
    v <- 0.95 - max(convertX(unit(1, "strwidth", as.character(c(x))), "native", valueOnly = TRUE) ) / 2

    fontsize = fontsize * scale

    for(i in 1 : k) {

        tab <- x[ , , i]

        fit <- stdize(tab, std, x)

        xInd <- indexMatrix[i, 2]
        xOrig <- 2 * xInd - 1 + (3 * xInd - 2) * space
        yInd <- indexMatrix[i, 1]
        yOrig <- if(k == 1)
            (1 + space)
        else
            (totalHeight
             - (2 * yInd - 1 + ((3 + gamma) * yInd - 2) * space))
        pushViewport(viewport(xscale = xlim - xOrig, yscale = ylim - yOrig))

        ## drawLabels()
        u <- 1 + space / 2
        adjCorr <- 0.2
        grid.text(
                  paste(names(dimnames(x))[1],
                        dimnames(x)[[1]][1],
                        sep = sep),
                  0, u,
                  gp = gpar(fontsize = fontsize),
                  default.units = "native"
             )
        grid.text(
                  paste(names(dimnames(x))[2],
                        dimnames(x)[[2]][1],
                        sep = sep),
                  -u, 0,
                  default.units = "native",
                  gp = gpar(fontsize = fontsize),
                  rot = 90)
        grid.text(
                  paste(names(dimnames(x))[1],
                        dimnames(x)[[1]][2],
                        sep = sep),
                  0, -u,
                  gp = gpar(fontsize = fontsize),
                  default.units = "native"
             )
        grid.text(
                  paste(names(dimnames(x))[2],
                        dimnames(x)[[2]][2],
                        sep = sep),
                  u, 0,
                  default.units = "native",
                  gp = gpar(fontsize = fontsize),
                  rot = 90)
        if (k > 1) {
            grid.text(if (!varnames)
                          dimnames(x)[[3]][i]
                      else
                          paste(names(dimnames(x))[3],
                                dimnames(x)[[3]][i],
                                sep = sep),
                      0, 1 + (1 + gamma / 2) * space,
                      gp = gpar(fontsize = fontsize * gamma),
                      default.units = "native"
                      )
          }

        ## drawFrequencies()

        ### in extended plots, emphasize charts with significant logoddsratios
        emphasize <- if(extended && is.numeric(conf_level))
          2 * extended * (1 + (p.lor.test[i] < 1 - conf_level))
        else 0

        d <- odds(tab)$or
        drawPie(sqrt(fit[1,1]),  90, 180, col = color[1 + (d > 1) + emphasize])
        drawPie(sqrt(fit[2,1]), 180, 270, col = color[2 - (d > 1) + emphasize])
        drawPie(sqrt(fit[1,2]),   0,  90, col = color[2 - (d > 1) + emphasize])
        drawPie(sqrt(fit[2,2]), 270, 360, col = color[1 + (d > 1) + emphasize])

        u <- 1 - space / 2
        grid.text(as.character(c(tab))[1],
                  -v, u,
                  just = c("left", "top"),
                  gp = gpar(fontsize = fontsize),
                  default.units = "native")
        grid.text(as.character(c(tab))[2],
                  -v, -u,
                  just = c("left", "bottom"),
                  gp = gpar(fontsize = fontsize),
                  default.units = "native")
        grid.text(as.character(c(tab))[3],
                  v, u,
                  just = c("right", "top"),
                  gp = gpar(fontsize = fontsize),
                  default.units = "native")
        grid.text(as.character(c(tab))[4],
                  v, -u,
                  just = c("right", "bottom"),
                  gp = gpar(fontsize = fontsize),
                  default.units = "native")

        ## draw ticks
        if(extended && ticks)
          if(d > 1) {
            grid.lines(c(sqrt(fit[1,1])           * cos(3*pi/4),
                         (sqrt(fit[1,1]) + ticks) * cos(3*pi/4)),
                       c(sqrt(fit[1,1])           * sin(3*pi/4),
                         (sqrt(fit[1,1]) + ticks) * sin(3*pi/4)),
                       gp = gpar(lwd = 1),
                       default.units = "native"
                       )
            grid.lines(c(sqrt(fit[2,2])           * cos(-pi/4),
                         (sqrt(fit[2,2]) + ticks) * cos(-pi/4)),
                       c(sqrt(fit[2,2])           * sin(-pi/4),
                         (sqrt(fit[2,2]) + ticks) * sin(-pi/4)),
                       gp = gpar(lwd = 1),
                       default.units = "native"
                       )
          } else {
            grid.lines(c(sqrt(fit[1,2])           * cos(pi/4),
                         (sqrt(fit[1,2]) + ticks) * cos(pi/4)),
                       c(sqrt(fit[1,2])           * sin(pi/4),
                         (sqrt(fit[1,2]) + ticks) * sin(pi/4)),
                       gp = gpar(lwd = 1),
                       default.units = "native"
                       )
            grid.lines(c(sqrt(fit[2,1])           * cos(-3*pi/4),
                         (sqrt(fit[2,1]) + ticks) * cos(-3*pi/4)),
                       c(sqrt(fit[2,1])           * sin(-3*pi/4),
                         (sqrt(fit[2,1]) + ticks) * sin(-3*pi/4)),
                       gp = gpar(lwd = 1),
                       default.units = "native"
                       )
          }

        ## drawConfBands()
        if(is.numeric(conf_level)) {
            or <- o$or[i]
            se <- o$se[i]
            ## lower
            theta <- or * exp(qnorm((1 - conf_level) / 2) * se)
            tau <- findTableWithOAM(theta, tab)
            r <- sqrt(c(stdize(tau, std, x)))
            for(j in 1 : 4)
                drawPie(r[j], angle.f[j], angle.t[j])
            ## upper
            theta <- or * exp(qnorm((1 + conf_level) / 2) * se)
            tau <- findTableWithOAM(theta, tab)
            r <- sqrt(c(stdize(tau, std, x)))
            for(j in 1 : 4)
                drawPie(r[j], angle.f[j], angle.t[j])
        }

        ## drawBoxes()
        grid.polygon(c(-1,  1, 1, -1),
                     c(-1, -1, 1,  1),
                     default.units = "native",
                     gp = gpar(fill = "transparent")
                     )
        grid.lines(c(-1, 1), c(0, 0), default.units = "native")
        for(j in seq(from = -0.8, to = 0.8, by = 0.2))
            grid.lines(c(j, j), c(-0.02, 0.02), default.units = "native")
        for(j in seq(from = -0.9, to = 0.9, by = 0.2))
            grid.lines(c(j, j), c(-0.01, 0.01), default.units = "native")
        grid.lines(c(0, 0), c(-1, 1), default.units = "native")
        for(j in seq(from = -0.8, to = 0.8, by = 0.2))
            grid.lines(c(-0.02, 0.02), c(j, j), default.units = "native")
        for(j in seq(from = -0.9, to = 0.9, by = 0.2))
            grid.lines(c(-0.01, 0.01), c(j, j), default.units = "native")

        popViewport(1)

    }

    if(!is.null(main) || !is.null(sub)) {
        if (!is.null(main))
          grid.text(main,
                    y = unit(1, "npc") + unit(1, "lines"),
                    gp = gpar(fontsize = 20, fontface = 2))
        if (!is.null(sub))
          grid.text(sub,
                    y = unit(0, "npc") - unit(1, "lines"),
                    gp = gpar(fontsize = 20, fontface = 2))
        popViewport(1)
      }

    popViewport(1)

    if (return_grob)
        return(invisible(grid.grab()))
    else
        return(invisible(NULL))
}
