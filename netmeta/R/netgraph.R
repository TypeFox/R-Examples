netgraph <- function(x, seq = x$seq,
                     labels = rownames(x$TE.fixed),
                     cex = 1, col = "slateblue", offset = 0.0175,
                     scale = 1.10,
                     plastic, thickness, lwd = 5, lwd.min = lwd / 2.5, lwd.max = lwd * 4,
                     dim = "2d",
                     ##
                     highlight = NULL, col.highlight = "red2",
                     lwd.highlight = lwd, highlight.split = ":",
                     ##
                     multiarm = any(x$narms > 2),
                     col.multiarm = NULL,
                     alpha.transparency = 0.5,
                     ##
                     points = FALSE, col.points = "red",
                     cex.points = 1, pch.points = 20,
                     ##
                     start.layout = ifelse(dim == "2d", "circle", "eigen"),
                     eig1 = 2, eig2 = 3, eig3 = 4,
                     iterate,
                     tol = 0.0001, maxit = 500, allfigures = FALSE,
                     A.matrix = x$A.matrix,
                     N.matrix = sign(A.matrix),
                     ##
                     xpos = NULL, ypos = NULL, zpos = NULL,
                     ...) {
  
  
  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")
  
  dim <- meta:::setchar(dim, c("2d", "3d"))
  is_2d <- dim == "2d"
  is_3d <- !is_2d
  ##
  start.layout <- meta:::setchar(start.layout, c("eigen", "prcomp", "circle", "random"))
  ##
  if (!missing(seq) & is.null(seq))
    stop("Argument 'seq' must be not NULL.")
  ##
  if (!missing(labels) & is.null(labels))
    stop("Argument 'labels' must be not NULL.")
  
  if (missing(iterate))
    iterate <- ifelse(start.layout == "circle", FALSE, TRUE)
  
  
  if (missing(plastic))
    if (start.layout == "circle" & iterate == FALSE & is_2d)
      plastic <- TRUE
    else
      plastic <- FALSE
  
  
  if (missing(thickness)) {
    if (start.layout == "circle" & iterate == FALSE & plastic == TRUE) {
      thick <- "se.fixed"
      thickness <- "se.fixed"
    }
    else {
      thick <- "equal"
      thickness <- "equal"
    }
  }
  else {
    if (!is.matrix(thickness)) {
      if (length(thickness) == 1 & is.character(thickness))
        thick <- meta:::setchar(thickness,
                                c("equal", "number.of.studies",
                                  "se.fixed", "se.random", "w.fixed", "w.random"))
      ##
      else if (length(thickness) == 1 & is.logical(thickness)) {
        if (thickness)
          thick <- "se.fixed"
        else
          thick <- "equal"
      }
    }
    else {
      if ((dim(thickness)[1] != dim(A.matrix)[1]) |
          (dim(thickness)[2] != dim(A.matrix)[2]))
        stop("Dimension of argument 'A.matrix' and 'thickness' are different.")
      if (is.null(dimnames(thickness)))
        stop("Matrix 'thickness' must have row and column names identical to argument 'A.matrix'.")
      else {
        if (any(rownames(thickness) != rownames(A.matrix)))
          stop("Row names of matrix 'thickness' must be identical to argument 'A.matrix'.")
        if (any(colnames(thickness) != colnames(A.matrix)))
          stop("Column names of matrix 'thickness' must be identical to argument 'A.matrix'.")
      }
      ##
      W.matrix <- thickness
      thick <- "matrix"
    }
  }
  
  
  if (allfigures & is_3d) {
    warning("Argument 'allfigures' set to FALSE for 3-D network plot.")
    allfigures <- FALSE
  }
  
  
  if (is.null(seq) | !(start.layout == "circle" & iterate == FALSE)) {
    seq1 <- 1:length(labels)
    if (!missing(seq) & !is.null(seq) & (is.null(xpos) & is.null(ypos)))
      warning("Argument 'seq' only considered if start.layout=\"circle\" and iterate=FALSE.")
  }
  else {
    rn <- rownames(x$TE.fixed)
    seq1 <- charmatch(setseq(seq, rn), rn)
  }
  ##
  A.matrix <- A.matrix[seq1, seq1]
  N.matrix <- N.matrix[seq1, seq1]
  ##
  if (thick == "matrix")
    W.matrix <- W.matrix[seq1, seq1]
  ##
  labels <- labels[seq1]
  
  
  A.sign <- sign(A.matrix)
  
  
  if ((is_2d & (is.null(xpos) & is.null(ypos))) |
      (is_3d & (is.null(xpos) & is.null(ypos) & is.null(zpos)))) {
    stressdata <- stress(x,
                         A.matrix = A.matrix,
                         N.matrix = N.matrix,
                         ##
                         dim = dim,
                         start.layout = start.layout,
                         iterate = iterate,
                         eig1 = eig1, eig2 = eig2, eig3 = eig3,
                         tol = tol,
                         maxit = maxit,
                         ##
                         allfigures = allfigures,
                         ##
                         seq = seq,
                         ##
                         labels = labels,
                         cex = cex,
                         col = col,
                         offset = offset,
                         scale = scale,
                         ##
                         plastic = plastic,
                         thickness = thickness,
                         lwd = lwd,
                         lwd.min = lwd.min,
                         lwd.max = lwd.max,
                         ##
                         highlight = highlight,
                         col.highlight = col.highlight,
                         lwd.highlight = lwd.highlight,
                         highlight.split = highlight.split,
                         ## multiarm
                         col.multiarm = col.multiarm,
                         alpha.transparency = alpha.transparency,
                         ##
                         points = points, col.points = col.points,
                         cex.points = cex.points, pch.points = pch.points,
                         ##
                         ...)
    ##
    xpos <- stressdata$x
    ypos <- stressdata$y
    if (is_3d)
      zpos <- stressdata$z
  }
  
  
  if (allfigures)
    return(invisible(NULL))
  
  
  n <- dim(A.matrix)[1]
  d <- scale * max(abs(c(min(c(xpos, ypos), na.rm = TRUE),
                         max(c(xpos, ypos), na.rm = TRUE))))
  
  
  ## Generate dataset for plotting
  ##
  if (is_2d)
    pd <- data.frame(xpos, ypos, labels, seq)
  else
    pd <- data.frame(xpos, ypos, zpos, labels, seq)
  ##
  pd$adj1 <- NA
  pd$adj2 <- NA
  pd$adj3 <- NA
  ##
  pd$adj1[pd$xpos >= 0] <- 0
  pd$adj1[pd$xpos <  0] <- 1
  ##
  pd$adj2[pd$ypos >  0] <- 0
  pd$adj2[pd$ypos <= 0] <- 1
  ##
  if (!is_2d) {
    pd$adj3[pd$zpos >  0] <- 0
    pd$adj3[pd$zpos <= 0] <- 1
  }
  ##
  offset <- offset * 2 * d
  ##
  if (is_2d) {
    pd$xpos.labels <- pd$xpos - offset + 2 * (pd$adj1 == 0) * offset
    pd$ypos.labels <- pd$ypos - offset + 2 * (pd$adj2 == 0) * offset
  }
  else {
    pd$xpos.labels <- pd$xpos
    pd$ypos.labels <- pd$ypos
    pd$zpos.labels <- pd$zpos
  }
  
  
  ##
  ## Define coloured regions for multi-arm studies
  ##
  if (multiarm) {
    td1 <- data.frame(studies = x$studies, narms = x$narms)
    td1 <- td1[rev(order(td1$narms)), ]
    td1 <- td1[td1$narms > 2, ]
    multiarm.studies <- td1$studies
    ##
    n.multi <- length(multiarm.studies)
    ##
    missing.col.multiarm <- missing(col.multiarm)
    ##
    if (missing.col.multiarm | is.null(col.multiarm)) {
      ## Check for R package colorspace & use various gray values if
      ## not installed packages
      if (!any(as.data.frame(installed.packages())$Package == "colorspace"))
        col.polygon <- grDevices::rainbow(n.multi, alpha = alpha.transparency)
      else
        col.polygon <- colorspace::sequential_hcl(n.multi, alpha = alpha.transparency)
    }
    else {
      ##
      if (is.function(col.multiarm)) {
        mcname <- deparse(substitute(col.multiarm))
        ##
        csfun <- function(fcall, fname) {
          is.cs <- length(grep(fname, fcall)) > 0
          if (is.cs)
            meta:::is.installed.package("colorspace")
          is.cs
        }
        ##
        if (csfun(mcname, "rainbow_hcl"))
          col.polygon <- colorspace::rainbow_hcl(n.multi, start = 240, end = 60, alpha = alpha.transparency)
        else if (csfun(mcname, "sequential_hcl"))
          col.polygon <- colorspace::sequential_hcl(n.multi, alpha = alpha.transparency)
        else if (csfun(mcname, "diverge_hcl"))
          col.polygon <- colorspace::diverge_hcl(n.multi, alpha = alpha.transparency)
        else if (csfun(mcname, "heat_hcl"))
          col.polygon <- colorspace::heat_hcl(n.multi, alpha = alpha.transparency)
        else if (csfun(mcname, "terrain_hcl"))
          col.polygon <- colorspace::terrain_hcl(n.multi, alpha = alpha.transparency)
        else if (csfun(mcname, "diverge_hsv"))
          col.polygon <- colorspace::diverge_hsv(n.multi, alpha = alpha.transparency)
        else if (csfun(mcname, "choose_palette")) {
          fcolm <- colorspace::choose_palette(n = n.multi)
          col.polygon <- fcolm(n = n.multi)
        }
        else
          col.polygon <- sapply(n.multi, col.multiarm, alpha = alpha.transparency)
        ##
        if (csfun(mcname, "sequential_hcl") |
            csfun(mcname, "diverge_hcl") |
            csfun(mcname, "heat_hcl"))
          col.polygon <- rev(col.polygon)
      }
    }
    ##
    if (!missing.col.multiarm & is.character(col.multiarm)) {
      if (length(col.multiarm) > 1 & length(col.multiarm) != n.multi)
        stop("Length of argument 'col.multiarm' must be equal to one or the number of multi-arm studies: ", n.multi)
      col.polygon <- col.multiarm
    }
  }
  ##
  ## Define line width
  ##
  if (thick == "number.of.studies") {
    W.matrix <- lwd.max * A.matrix / max(A.matrix)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "equal") {
    W.matrix <- lwd * A.sign
  }
  else if (thick == "se.fixed") {
    IV.matrix <- x$seTE.direct.fixed[seq1, seq1]
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * min(IV.matrix, na.rm = TRUE) / IV.matrix
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "se.random") {
    IV.matrix <- x$seTE.direct.random[seq1, seq1]
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * min(IV.matrix, na.rm = TRUE) / IV.matrix
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "w.fixed") {
    IV.matrix <- 1 / x$seTE.direct.fixed[seq1, seq1]^2
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * IV.matrix / max(IV.matrix, na.rm = TRUE)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "w.random") {
    IV.matrix <- 1 / x$seTE.direct.random[seq1, seq1]^2
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * IV.matrix / max(IV.matrix, na.rm = TRUE)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "matrix") {
    W.matrix[is.infinite(W.matrix)] <- NA
    if (min(W.matrix[W.matrix != 0], na.rm = TRUE) == max(W.matrix[W.matrix != 0], na.rm = TRUE))
      W.matrix <- lwd * W.matrix
    else
      W.matrix <- lwd.max * W.matrix / max(W.matrix, na.rm = TRUE)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  
  
  
  
  
  ##
  ##
  ## Plot graph
  ##
  ##
  range <- c(-d, d)
  ##
  if (is_2d) {
    oldpar <- par(xpd = TRUE, pty = "s")
    on.exit(par(oldpar))
    ##
    plot(xpos, ypos,
         xlim = range, ylim = range,
         type = "n", axes = FALSE, bty = "n",
         xlab = "", ylab = "",
         ...)
    ##
    ## Add coloured regions for multi-arm studies
    ##
    if (multiarm) {
      ##
      if (n.multi > 0) {
        multiarm.labels <- vector("list", n.multi)
        if (length(col.polygon) == 1)
          col.polygon <- rep(col.polygon, n.multi)
        for (i in 1:n.multi) {
          treat1 <- x$treat1[x$studlab %in% multiarm.studies[i]]
          treat2 <- x$treat2[x$studlab %in% multiarm.studies[i]]
          multiarm.labels[[i]] <- sort(unique(c(treat2, treat1)))
          ##
          pdm <- pd[pd$seq %in% multiarm.labels[[i]], ]
          if (nrow(pdm) == 0)
            pdm <- pd[pd$labels %in% multiarm.labels[[i]], ]
          ##
          ## Clockwise ordering of polygon coordinates
          ##
          polysort <- function(x, y) {
            xnorm <- (x - mean(x)) / sd(x) # Normalise coordinate x
            ynorm <- (y - mean(y)) / sd(y) # Normalise coordinate y
            r <- sqrt(xnorm^2 + ynorm^2)   # Calculate polar coordinates
            cosphi <- xnorm / r
            sinphi <- ynorm / r
            s <- as.numeric(sinphi > 0) # Define angles to lie in [0, 2 * pi]
            phi <- acos(cosphi)
            alpha <- s * phi + (1 - s) * (2 * pi - phi)
            ##
            res <- order(alpha)
            res
          }
          ##
          pdm <- pdm[polysort(pdm$xpos, pdm$ypos), ]
          ##
          polygon(pdm$xpos, pdm$ypos,
                  col = col.polygon[i], border = NA)
        }
      }
    }
    ##
    ## Draw lines
    ##
    if (plastic) {
      n.plastic <- 30
      lwd.multiply <- rep(NA, n.plastic)
      cols <- rep("", n.plastic)
      j <- 0
      for (i in n.plastic:1) {
        j <- j + 1
        lwd.multiply[j] <- sin(pi * i / 2 / n.plastic)
        cols[j] <- paste("gray", round(100 * (1 - i / n.plastic)), sep = "")
      }
    }
    else {
      lwd.multiply <- 1
      cols <- col
    }
    ##
    for (n.plines in 1:length(lwd.multiply)) {
      for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
          if (A.sign[i, j] > 0) {
            lines(c(xpos[i], xpos[j]), c(ypos[i], ypos[j]),
                  lwd = W.matrix[i, j] * lwd.multiply[n.plines],
                  col = cols[n.plines])
          }
        }
      }
    }
    ##
    ## Add highlighted comparisons
    ##
    if (!is.null(highlight)) {
      for (high in highlight) {
        highs <- unlist(strsplit(high, split = highlight.split))
        if (length(highs) != 2)
          stop("Wrong format for argument 'highlight' (see helpfile of plotgraph command).")
        ##
        if (sum(pd$labels %in% highs) != 2)
          stop(paste("Argument 'highlight' must contain two of the following values (separated by \":\"):\n  ",
                     paste(paste("'", pd$labels, "'", sep = ""),
                           collapse = " - "), sep = ""))
        ##
        pdh <- pd[pd$labels %in% highs, ]
        ##
        if (is_2d) {
          lines(pdh$xpos, pdh$ypos,
                lwd = W.matrix[labels == highs[1], labels == highs[2]],
                col = col.highlight)
        }
      }
    }
    ##
    ## Add points for labels
    ##
    if (points)
      points(xpos, ypos,
             pch = pch.points, cex = cex.points, col = col.points)
    ##
    ## Print treatment labels
    ##
    if (!is.null(labels))
      for (i in 1:n)
        text(pd$xpos.labels[i], pd$ypos.labels[i],
             labels = pd$labels[i],
             cex = cex,
             adj = c(pd$adj1[i], pd$adj2[i]))
  }
  else {
    plot3d(xpos, ypos, zpos,
           size = 10, col = col.points, cex = cex.points,
           axes = FALSE, box = FALSE,
           xlab = "", ylab = "", zlab = "")
    ##
    ## Add points for labels
    ##
    if (points)
      points3d(xpos, ypos, zpos,
               pch = pch.points, cex = cex.points, col = col.points)
    ##
    ## Print treatment labels
    ##
    if (!is.null(labels))
      for (i in 1:n)
        text3d(pd$xpos.labels[i], pd$ypos.labels[i], pd$zpos.labels[i],
               texts = pd$labels[i],
               cex = cex,
               adj = c(pd$adj1[i], pd$adj2[i]))
    ##
    ## Add highlighted comparisons
    ##
    if (!is.null(highlight)) {
      for (high in highlight) {
        highs <- unlist(strsplit(high, split = highlight.split))
        if (length(highs) != 2)
          stop("Wrong format for argument 'highlight' (see helpfile of plotgraph command).")
        ##
        if (sum(pd$labels %in% highs) != 2)
          stop(paste("Argument 'highlight' must contain two of the following values (separated by \":\"):\n  ",
                     paste(paste("'", pd$labels, "'", sep = ""),
                           collapse = " - "), sep = ""))
        ##
        pdh <- pd[pd$labels %in% highs, ]
        ##
        lines3d(pdh$xpos*(1+1e-4), pdh$ypos*(1+1e-4), pdh$zpos*(1+1e-4),
                lwd = W.matrix[labels == highs[1], labels == highs[2]],
                col = col.highlight)
      }
    }
    ##
    ## Add coloured regions for multi-arm studies
    ##
    if (multiarm) {
      ##
      morethan3 <- FALSE
      ##
      if (n.multi > 0) {
        multiarm.labels <- vector("list", n.multi)
        if (length(col.polygon) == 1)
          col.polygon <- rep(col.polygon, n.multi)
        for (i in 1:n.multi) {
          treat1 <- x$treat1[x$studlab %in% multiarm.studies[i]]
          treat2 <- x$treat2[x$studlab %in% multiarm.studies[i]]
          multiarm.labels[[i]] <- sort(unique(c(treat2, treat1)))
          ##
          pdm <- pd[pd$seq %in% multiarm.labels[[i]], ]
          if (nrow(pdm) == 0)
            pdm <- pd[pd$labels %in% multiarm.labels[[i]], ]
          if (nrow(pdm) == 3)
            triangles3d(pdm$xpos, pdm$ypos, pdm$zpos,
                        col = col.polygon[i])
          else
            morethan3 <- TRUE
        }
      }
      if (morethan3)
        warning("Multi-arm studies with more than three treatments not shown in 3-D plot.")
    }
    ##
    ## Draw lines
    ##
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (A.sign[i, j] > 0) {
          lines3d(c(xpos[i], xpos[j]), c(ypos[i], ypos[j]), c(zpos[i], zpos[j]),
                  lwd = W.matrix[i, j],
                  col = col)
        }
      }
    }
  }
  
  
  invisible(NULL)
}
