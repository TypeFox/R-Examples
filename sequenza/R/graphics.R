cp.plot <- function (cp.table, xlab = "Ploidy", ylab = "Cellularity", zlab = "Scaled rank LPP", 
                     colFn = colorRampPalette(c('white', 'lightblue')), ...) {
  z <- matrix(rank(cp.table$lpp), nrow = nrow(cp.table$lpp)) / length(cp.table$lpp)
  map <- makecmap(c(0, 1), colFn = colFn, include.lowest = TRUE)
  colorgram(x = cp.table$ploidy, y = cp.table$cellularity, z = z, 
            map = map, las = 1, 
            xlab = xlab, ylab = ylab, zlab = zlab, ...)
}

cp.plot.contours <- function(cp.table, likThresh = c(0.95), alternative = TRUE,
                             col = palette(), legend.pos = 'bottomright', pch = 18, alt.pch = 3, ...) {
   znormsort <- sort(cp.table$lpp, decreasing = TRUE)
   znormcumLik <- cumsum(znormsort)
   n <- sapply(likThresh, function(x) sum(znormcumLik < x) + 1)
   LikThresh <- znormsort[n]
   names(LikThresh) <- paste0(likThresh * 100, '%')
   contour(x = cp.table$ploidy, y = cp.table$cellularity, z = cp.table$lpp,
           levels = znormsort[n], col = col, drawlabels = FALSE,
           xlab = "Ploidy", ylab = "Cellularity", ...)
   max.xy <- which(cp.table$lpp == max(cp.table$lpp), arr.ind = TRUE)
   points(x = cp.table$ploidy[max.xy[, "row"]],
          y = cp.table$cellularity[max.xy[, "col"]], pch = pch)
   if (alternative == TRUE){
      alt.sol <- alternative.cp.solutions(cp.table)
      alt.sol <- alt.sol[-1, ]
      points(x = alt.sol$ploidy,
             y = alt.sol$cellularity, pch = alt.pch)
   }
   if(!is.na(legend.pos)) {
      if (alternative == FALSE) {
         legend(legend.pos, legend = c(paste("C.R.", names(LikThresh), sep = " "), "Point estimate"),
                col = c(col[1:length(LikThresh)], "black"), lty = c(rep(1, length(LikThresh)), NA),
                pch = c(rep(NA, length(LikThresh)), pch), border = NA, bty = "n")
      } else {
         legend(legend.pos, legend = c(paste("C.R.", names(LikThresh), sep = " "),
                                       "Point estimate", "Alternative solutions"),
                col = c(col[1:length(LikThresh)], "black", "black"), lty = c(rep(1, length(LikThresh)), NA, NA),
                pch = c(rep(NA, length(LikThresh)), pch, alt.pch), border = NA, bty = "n")         
      }
   }
   invisible(LikThresh)
}

# plot.fit.model <- function(mufreq.tab, cellularity, ploidy, chr23 = "XY",
#                            cn.ratio.range = c(0.5:2), avg.depth.ratio = avg.depth.ratio,
#                            cex.m = 1, cex.d = 1, ...) {
#    xy.index   <- mufreq.tab$chr == "chrX" | mufreq.tab$chr == "chrY"
#    plot(x = mufreq.tab$F, y = mufreq.tab$adjusted.ratio,
#         xlab = "mutation frequency", ylab = "depth.ratio",
#         las = 1, type="n", ...)
#    points(x = mufreq.tab$F[!xy.index], y = mufreq.tab$adjusted.ratio[!xy.index],
#           pch = 19, col = "blue", cex = cex.d)
#    types      <- types.matrix(cn.ratio.range = cn.ratio.range, chr23 = chr23)
#    if (length(which(xy.index)) >= 1) {
#       if (chr23 == "XY") {
#          points(x = mufreq.tab$F[xy.index], y = mufreq.tab$adjusted.ratio[xy.index], pch = 19, col = "green")
#       } else {
#          points(x = mufreq.tab$F[xy.index], y = mufreq.tab$adjusted.ratio[xy.index], pch = 19, col = "blue")
#       }
#    }
#    points.fit <-model.points(cellularity = cellularity, ploidy = ploidy,
#                              types = types, avg.depth.t = avg.depth.t,
#                              avg.depth.r = avg.depth.r)
#    points(points.fit,pch = 19, col = "red", cex = cex.m)
#    if (length(which(xy.index)) >= 1 & chr23 == "XY") {
#       legend("bottomright", c("Male chr X/Y",paste(paste("C", cellularity,sep = " : "), paste("P", ploidy,sep = " : "), sep = "; ")), pch = 19, col = c("green","red"))
#    } else {
#       legend("bottomright", paste(paste("C", cellularity,sep = " : "), paste("P", ploidy,sep = " : "), sep = "; "), pch = 19, col = "red")
#    }
# }

# plot.gc.depth <- function(d.values, gc.contents, ...) {
#    gc.values <- sort(unique(gc.contents))
#    d.list <- list()
#    #sizevect <- rep(0,length(gc.values))
#    for ( i in 1:length(gc.values)) {
#       d.list[[i]] <- d.values[gc.contents == gc.values[i]]
#       #sizevect[i]    <- length(d.list[[i]])
#    }
#    bxplot(d.list, names = gc.values, ...)
# }

plotWindows <- function(seqz.window, m.lty = 1, m.lwd = 3,
                         m.col = "black", q.bg = "lightblue", log2.plot = FALSE,
                         n.min = 1, xlim, ylim, add = FALSE, ...) {
   if (log2.plot) {
      seqz.window[, c(3, 4, 5)] <- log2(seqz.window[, c(3, 4, 5)])
   }
   if(!add) {
      if(missing(xlim))
         xlim <- c(seqz.window$start[1], seqz.window$end[nrow(seqz.window)])
      if(missing(ylim))
         ylim <- c(min(seqz.window$q0, na.rm = TRUE), max(seqz.window$q1, na.rm = TRUE))
      plot(xlim, ylim, type = "n", ...)
   }
   seqz.window <- seqz.window[seqz.window$N >= n.min, ]
   rect(xleft = seqz.window$start, ybottom = seqz.window$q0,
        xright = seqz.window$end, ytop = seqz.window$q1,
        col = q.bg, border = NA)
   segments(y0 = seqz.window$mean, x0 = seqz.window$start, x1 = seqz.window$end, 
            lty = m.lty, lwd = m.lwd, col = m.col)

}

chromosome.view <- function(baf.windows, ratio.windows, mut.tab = NULL, segments = NULL,  min.N.baf = 1, min.N.ratio = 1e4,
                            main = "", vlines = FALSE, legend.inset = c(-20 * strwidth("a", units = 'figure'), 0), BAF.style = "lines",
                            CNn = 2, cellularity = NULL, ploidy = NULL, avg.depth.ratio = NULL, model.lwd = 1,
                            model.lty = "24", model.col = 1, x.chr.space = 10) {
   make.polygons <- function(segments, model.baf) {
      max.B      <- max(model.baf$B[model.baf$CNt == max(segments$CNt)])
      mat.polygs <- matrix(ncol = max.B+1, nrow = nrow(segments))
      colnames(mat.polygs) <- 0:max.B
      get.B <- function (CNt, B) model.baf$BAF[model.baf$CNt == CNt & model.baf$B == B]
      polyg.coords <- sapply( 0:max.B, FUN = function (k) as.numeric(sapply(segments$CNt, FUN = function(i) get.B(i, k))))
      polyg.coords[is.na(polyg.coords)] <- 1
      polyg.coords <- cbind(polyg.coords, 1)
      polyg.pos    <- segments[, c("start.pos", "end.pos")]
      edge1        <- polyg.pos$end.pos[-nrow(polyg.pos)]
      edge2        <- polyg.pos$start.pos[-1]
      no.dat       <- c(1:nrow(polyg.pos))[edge2 - edge1 >= 1e6]
      v.gaps       <- apply(rbind(edge1[no.dat], edge2[no.dat]), 2, mean)
      v.gaps       <- cbind(start.pos = v.gaps, end.pos = v.gaps)
      polyg.p.new  <- rbind(polyg.pos, v.gaps)
      polyg.c.new  <- polyg.coords
      for (i in no.dat) {
         polyg.c.new  <- rbind(polyg.c.new, 0)
      }
      new_idx      <- c(seq_along(polyg.pos$start.pos), no.dat+0.5)
      polyg.pos    <- polyg.p.new[order(new_idx), ]
      polyg.coords <- polyg.c.new[order(new_idx), ]
      Xs      <- unlist(lapply(1:nrow(polyg.coords), function(k) polyg.pos[k, ]))
      #color <- colorRampPalette(c("grey99", "grey20"))( max.B + 1 )
      color <- gray.colors((max.B + 1), start = 0.5, end = 0.9, alpha = 0.3)
      extra.x <- c(max(Xs), min(Xs))
      extra.y = c(0 ,0)
      for (k in 1:ncol(polyg.coords)) {
         Ys <- unlist(lapply(polyg.coords[, k], function (i) rep(i, 2)))
         polygon(x = c(Xs,extra.x), y = c(Ys,extra.y), col = color[k], border = NA)
         extra.y <- rev(Ys)
         extra.x <- rev(Xs)
      }
   }
   if (is.null(segments)) {
      data.model <- NULL
   } else {
      if ("CNt" %in% colnames(segments)) {
         if (length(c(cellularity, ploidy, avg.depth.ratio)) != 3) {
            data.model <- NULL
         } else {
            data.model     <- list()
            CNt.max        <- max(segments$CNt, na.rm = TRUE) + 1
            CNt.min        <- 0
            data.model$baf <- expected.baf(sd = mean(segments$sd.BAF, na.rm = TRUE), CNn = CNn, CNt = CNt.max, cellularity = cellularity)
            if (CNn == 2) {
               data.model$baf <- rbind(c(0,0,max(data.model$baf$BAF),0), data.model$baf)
            } else {
               data.model$baf <- rbind(c(0,0,1,0), data.model$baf)
            }       
            types          <- types.matrix(CNt.min = CNt.min, CNt.max = CNt.max, CNn = CNn)
            data.model$muf <- cbind(types, model.points(cellularity = cellularity, ploidy = ploidy,
                                                   types = types, avg.depth.ratio = avg.depth.ratio))
         }
      } else {
         data.model <- NULL
      }
   }
   if (is.null(mut.tab)) {
      par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(2,1), xaxt='n')
      min.x <- min(c(min(baf.windows$start), min(ratio.windows$start)))
      max.x <- max(c(max(baf.windows$end), max(ratio.windows$end)))
      xlim <- c(min.x, max.x)
   } else {
      min.x <- min(c(min(baf.windows$start), min(ratio.windows$start), min(mut.tab$position)))
      max.x <- max(c(max(baf.windows$end), max(ratio.windows$end), max(mut.tab$position)))
      xlim <- c(min.x, max.x)
      par(mar = c(0, 4, 0, 10), oma = c(5, 0, 4, 0), mfcol = c(3,1), xaxt='n', xpd = TRUE)
      mutation.colors <- c(
         'A>C' = rgb(red =   0, green = 178, blue = 238, alpha = 120, maxColorValue = 255),
         'T>G' = rgb(red =   0, green = 178, blue = 238, alpha = 120, maxColorValue = 255),
         'A>G' = rgb(red = 255, green =  64, blue =  64, alpha = 120, maxColorValue = 255),
         'T>C' = rgb(red = 255, green =  64, blue =  64, alpha = 120, maxColorValue = 255),
         'A>T' = rgb(red =  34, green = 139, blue =  34, alpha = 120, maxColorValue = 255),
         'T>A' = rgb(red =  34, green = 139, blue =  34, alpha = 120, maxColorValue = 255),
         'C>A' = rgb(red = 139, green =  90, blue =   0, alpha = 120, maxColorValue = 255),
         'G>T' = rgb(red = 139, green =  90, blue =   0, alpha = 120, maxColorValue = 255),
         'C>G' = rgb(red = 127, green =   0, blue = 255, alpha = 120, maxColorValue = 255),
         'G>C' = rgb(red = 127, green =   0, blue = 255, alpha = 120, maxColorValue = 255),
         'C>T' = rgb(red = 255, green = 215, blue =   0, alpha = 120, maxColorValue = 255),
         'G>A' = rgb(red = 255, green = 215, blue =   0, alpha = 120, maxColorValue = 255)
      )
      plot(x = mut.tab$position, y = mut.tab$F,
           ylab = "Mutant allele frequency", las = 1, pch = 19,
           col = c(mutation.colors, 'NA' = NA)[as.character(mut.tab$mutation)],
           ylim = c(min(mut.tab$F, na.rm = TRUE), 1), xlim = xlim)
      unique.colors <- unique(mutation.colors)
      labels <- sapply(unique.colors, function(a) paste(names(mutation.colors)[mutation.colors == a], collapse = ", "))
      #legend("topleft", legend = labels, fill = unique.colors, border = NA, bty = "n")
      legend(y = "center", x  = "right", legend = labels,
             inset = legend.inset, pch = 19, col = unique.colors,
             pt.bg = unique.colors, border = NA, bty = "n")
      if (!is.null(segments)){
         if (vlines) {
            abline(v = segments$end.pos, lwd = 1, lty = 2)
         }
         if (!is.null(data.model)) {
            for (i in 1:nrow(segments)) {
                segments(x0 = segments$start.pos[i], x1 = segments$end.pos[i],
                         y0 = unique(data.model$muf$mufreqs[data.model$muf$CNt == segments$CNt[i]]), lwd = model.lwd, lty = model.lty, col = model.col)
            }
         }
      }

   }
   if (!is.null(segments)){
      plot(ylab = "B allele frequency", type = "n",
           x = xlim, y = c(0, 0.5), las = 1)
      if (BAF.style == "blocks") {
         if (!is.null(data.model)) {
            make.polygons(segments, data.model$baf)
            axis(side = 4, line = 0, las = 1,
                 labels = data.model$baf$B[data.model$baf$CNt == segments$CNt[nrow(segments)]],
                 at = data.model$baf$BAF[data.model$baf$CNt == segments$CNt[nrow(segments)]])
            mtext(text = "Number of B alleles", side = 4, line = 2, cex = par("cex.lab")*par("cex"))
         }
      }
      plotWindows(baf.windows, ylab = "B allele frequency",
                  xlim = xlim, ylim = c(0, 0.5), las = 1,
                  n.min = min.N.baf, add = TRUE)
      if (vlines) {
         abline(v = segments$end.pos, lwd = 1, lty = 2)
      }
      segments(x0 = segments$start.pos, y0 = segments$Bf, x1=segments$end.pos, y1 = segments$Bf, col = "red", lwd = 3)
      if (BAF.style == "lines") {
         if (!is.null(data.model)) {
            for (i in 1:nrow(segments)) {
               segments(x0 = segments$start.pos[i], x1 = segments$end.pos[i],
                        y0 = unique(data.model$baf$BAF[data.model$baf$CNt == segments$CNt[i]]), lwd = model.lwd, lty = model.lty, col = model.col)
            }
         }
      }
   }
   else {
      plotWindows(baf.windows, ylab = "B allele frequency",
                  xlim = xlim, ylim = c(0, 0.5), las = 1,
                  n.min = min.N.baf)
   }
   plotWindows(ratio.windows, ylab = "Depth ratio",
               las = 1, n.min = min.N.ratio, ylim = c(0, 2.5))
   if (!is.null(segments)){
      if (vlines) {
         abline(v = segments$end.pos, lwd = 1, lty = 2)
      }
      segments(x0 = segments$start.pos, y0 = segments$depth.ratio, x1=segments$end.pos, y1 = segments$depth.ratio, col = "red", lwd = 3)
      if (!is.null(data.model)) {
         ratios.theoric <- unique(data.model$muf[,c('CNt', 'depth.ratio')])

         segments(x0 = rep(min(segments$start.pos, na.rm =TRUE), times = nrow(ratios.theoric)),
                  x1 = rep(max(segments$end.pos, na.rm = TRUE), times = nrow(ratios.theoric)),
                  y0 = ratios.theoric$depth.ratio, lwd = model.lwd, lty = model.lty, col = model.col)
         #text(x = rep(min(segments$start.pos, na.rm =TRUE), times = nrow(ratios.theoric)),
         #     y = ratios.theoric$depth.ratio, labels = ratios.theoric$CNt, pos = 2, offset = 0.5, cex = 0.8)
         axis(labels = as.character(ratios.theoric$CNt), side = 4, line = 0, las = 1,
              at = ratios.theoric$depth.ratio)
         mtext(text = "Copy number", side = 4, line = 2, cex = par("cex.lab")*par("cex"))
      }
   }
   par(xaxt='s')
   axis(labels = as.character(round(seq(xlim[1]/1e6, xlim[2]/1e6, by = x.chr.space), 0)), side = 1 , line = 0,
        at = seq(xlim[1], xlim[2], by = 1e6 * x.chr.space), outer = FALSE, cex = par("cex.axis")*par("cex"))
   mtext("Position (Mb)", side = 1, line = 3, outer = FALSE, cex = par("cex.lab")*par("cex"))
   mtext(main, 3, outer = TRUE, cex = par("cex.main")*par("cex"), line = 2)
}

#genome.view <- function(baf.windows, ratio.windows, segments = NULL, main = "",
#                            min.N.baf = 1, min.N.ratio = 1e4, CNn = rep(2, length(ratio.windows)),
#                            cellularity = NULL, ploidy = NULL, avg.depth.ratio = NULL) {
#   chr.metrics <- list()
#   for (i in 1:length(ratio.windows)) {
#      chr.metrics[[i]] <- range(ratio.windows[[i]]$mean, na.rm = TRUE)
#   }
#   chr.metrics <- do.call(rbind, chr.metrics)
#   x0 <- chr.metrics[1,1]
#
#}

genome.view <- function(seg.cn, info.type = "AB", ...) {
   chr.order <- unique(seg.cn$chromosome)
   seg.list  <- split(x = seg.cn[,c("chromosome", "start.pos", "end.pos", "A", "B", "CNt")],
                      f = seg.cn$chromosome)
   seg.list  <- seg.list[order(order(chr.order))]
   seg.max   <- lapply(X = seg.list, FUN = function(x) x[nrow(x), "end.pos" ])
   seg.pos   <- lapply(seg.list, "[", TRUE, c("start.pos", "end.pos"))
   seg.max   <- cumsum(as.numeric(do.call(rbind, seg.max)))
   chr.offset <- 0
   for (i in 1:length(seg.pos)){
      seg.pos[[i]] <- seg.pos[[i]] + chr.offset
      colnames(seg.pos[[i]]) <- c("abs.start","abs.end")
      chr.offset   <- seg.max[i]
   }
   seg.max      <- sapply(X = seg.pos, FUN = function(x) x[nrow(x), "abs.end" ])
   abs.list     <- mapply(cbind, seg.list, seg.pos, SIMPLIFY = FALSE)
   abs.segments <- do.call(rbind, abs.list)
   if (info.type == "AB") {
      abs.segments <- na.exclude(abs.segments)
      plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
           y = c(-0.1, (max(abs.segments$A)+0.1)), type = "n",
           ylab = "Copy number", xlab = "Position (Mb)",
           xaxt='n',  yaxt='n', xaxs = "i", ...)
      axis(labels = 0:max(abs.segments$A),
           at = 0:max(abs.segments$A),
           side = 2, line = 0, las = 1)
      #abline(h = c(0:max(abs.segments$A)), lty = 2)
      segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
               y0 = (abs.segments$B-0.1), y1 = (abs.segments$B-0.1), col="blue", lwd = 5, lend = 1)
      segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
               y0 = (abs.segments$A+0.1), y1 = (abs.segments$A+0.1), col="red", lwd = 5, lend = 1)
   } else {
      abs.segments <- abs.segments[!is.na(abs.segments$CNt), ]
      plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
           y = c(min(abs.segments$CNt), max(abs.segments$CNt)), type = "n",
           ylab = "Copy number", xlab = "Position (Mb)",
           xaxt='n', yaxt = 'n', xaxs = "i", ...)
      axis(labels = min(abs.segments$CNt):max(abs.segments$CNt),
           at = min(abs.segments$CNt):max(abs.segments$CNt),
           side = 2, line = 0, las = 1)
      #abline(h = c(min(abs.segments$CNt):max(abs.segments$CNt)), lty = 2)
      segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
               y0 = abs.segments$CNt, y1= abs.segments$CNt, col="red", lwd = 5, lend = 1)
   }
   abline(v = c(0, seg.max), lty = 3)
   for (i in 1:length(abs.list)){
      max.pos <- nrow(abs.list[[i]])
      mtext(chr.order[i], side = 3, line = 0,
            at = sum(abs.list[[i]]$abs.start[1], abs.list[[i]]$abs.end[max.pos])/2)
      #axis(labels = as.character(round(seq(abs.list[[i]]$start.pos[1]/1e6, abs.list[[i]]$end.pos[max.pos]/1e6, by = 20), 0)),
      #     at = seq(abs.list[[i]]$abs.start[1], abs.list[[i]]$abs.end[max.pos], by = 2e7), outer = FALSE, cex = par("cex.axis")*par("cex"),
      #     side = 1 , line = 0)
   }
   axis(labels = as.character(round(seq(abs.list[[1]]$start.pos[1]/1e6, abs.list[[1]]$end.pos[nrow(abs.list[[1]])]/1e6, by = 50), 0)),
        at = seq(abs.list[[1]]$abs.start[1], abs.list[[1]]$abs.end[nrow(abs.list[[1]])], by = 5e7), outer = FALSE, cex = par("cex.axis")*par("cex"),
        side = 1 , line = 1)
}

baf.ratio.model.fit <- function(cellularity, ploidy, segs, BAF.space = seq(0.001, 0.5, 0.005),
                                ratio.space = seq(0.01, 2.5, 0.05), avg.depth.ratio = 1, CNt.max = 7,
                                segment.filter = 3e6) {
   s.b   <- mean(segs$sd.BAF, na.rm = TRUE)
   s.r   <- mean(segs$sd.ratio, na.rm = TRUE)
   l.s   <- segs$end.pos - segs$start.pos
   s.big <- l.s >= segment.filter
   test.values <- expand.grid(Bf = BAF.space, ratio = ratio.space,
                              KEEP.OUT.ATTRS = FALSE)
   both.space  <- baf.bayes(Bf = test.values$Bf, CNt.max = CNt.max, CNt.min = 0,
                           depth.ratio = test.values$ratio,
                           cellularity = cellularity, ploidy = ploidy,
                           avg.depth.ratio = avg.depth.ratio,
                           sd.Bf = s.b, weight.Bf = 10,
                           sd.ratio = s.r, weight.ratio = 10, ratio.priority = F,
                           CNn = 2) 
   both.space <- as.data.frame(both.space)
   z <- tapply(both.space$LPP, list(test.values$Bf, test.values$ratio), mean)
   x <- as.numeric(rownames(z))
   y <- as.numeric(colnames(z))
   t <- types.matrix(CNt.min = 0, CNt.max = CNt.max, CNn = 2)
   mpts <- cbind(t,
                 model.points(cellularity = cellularity,
                              ploidy = ploidy, types = t,
                              avg.depth.ratio = avg.depth.ratio)
   )
   mpts <- unique(mpts[, c("CNt", "depth.ratio")])
   par(mar = c(5.1, 4.1, 4.1, 4.1))
   rev.heat <- function(...){rev(heat.colors(...))}
   suppressWarnings(colorgram(x, y, z, key = NA, nz = 1000, xlab = "B allele frequency", ylab = "Depth ratio",
             main = paste("cellularity:", cellularity, "ploidy:", ploidy, "sd.BAF:", round(s.b,2), sep = " "),
             map = makecmap(z, breaks = unique(quantile(z, seq(.25,1,0.0001))), right = TRUE, n = 1000, colFn = rev.heat),
             outlier = "white", las = 1, xlim = c(0, 0.5)))
   axis(side = 4, at = mpts$depth.ratio, labels = mpts$CNt, las = 1)
   mtext(text = "Copy number", side = 4, line = 2)
   points(x = segs$Bf[s.big], y = segs$depth.ratio[s.big], pch = 1, cex = 1)
   points(x = segs$Bf[!s.big], y = segs$depth.ratio[!s.big], pch = ".", cex = 1)
}

plotRawGenome <- function(sequenza.extract, cellularity, ploidy, CNt.max = 7, main = "", ...){
   max.end <- sapply(sequenza.extract$ratio, FUN = function(x) max(x$end, na.rm = T))
   max.end <- c(0, cumsum(as.numeric(max.end)))
   chrs <- names(sequenza.extract$ratio)
   coords.names <- (max.end + c(diff(max.end)/2,0))[1:length(chrs)]
   new.coords <- function(win.list, max.end){
      lapply(1:length(win.list), FUN = function(x) {
         y <- win.list[[x]]
         y$start <- y$start + max.end[x]
         y$end <- y$end + max.end[x]
         y
      })}
   new.coords.segs <- function(segs, max.end){
      lapply(1:length(segs), FUN = function(x) {
         y <- segs[[x]]
         y$start.pos <- y$start.pos + max.end[x]
         y$end.pos <- y$end.pos + max.end[x]
         y
      })}
   ratio.new <- new.coords(sequenza.extract$ratio,max.end)
   BAF.new   <- new.coords(sequenza.extract$BAF,max.end)
   segs.new  <- do.call(rbind, new.coords.segs(sequenza.extract$segments,max.end))
   avg.depth.ratio <- 1
   #avg.depth.ratio = mean(sequenza.extract$gc$adj[,2])
   #avg.depth.ratio <- center.ratio(segs.new)
   par(mar = c(1, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(2,1), ...)
   plot(x = c(min(max.end), max(max.end)), y = c(0,0.5), main = main, xlab = NA,
        ylab = "B allele frequency", type = "n", las = 1, xaxs = "i", yaxs = "i", xaxt = "n" )
   plotWindows(seqz.window = do.call(rbind, BAF.new), q.bg = "lightblue", m.col = "black", add = T)
   segments(x0 = segs.new$start.pos, x1 = segs.new$end.pos,
            y0 = (segs.new$Bf), y1 = (segs.new$Bf), col = "red", lwd = 2, lend = 1)
   abline(v = max.end, lty = 1)
   plot(x=c(min(max.end), max(max.end)), y = c(0,2.5), main = "", xlab = NA,
        ylab = "Depth ratio", type = "n", las = 1, xaxs = "i", yaxs = "i", xaxt = "n")
   plotWindows(seqz.window = do.call(rbind, ratio.new), q.bg = "lightblue", m.col = "black", add = T)
   segments(x0 = segs.new$start.pos, x1 = segs.new$end.pos,
            y0 = (segs.new$depth.ratio), y1 = (segs.new$depth.ratio), col = "red", lwd = 2, lend = 1)
   if (!missing(ploidy) & !missing(cellularity)){
      types        <- types.matrix(CNt.min = 0, CNt.max = CNt.max, CNn = 2)
      depth.ratios <- model.points(cellularity = cellularity, ploidy = ploidy,
                                   avg.depth.ratio = avg.depth.ratio,
                                   types = types)[, "depth.ratio"]
      depth.ratios <- unique(data.frame(CNt = types$CNt, ratio = depth.ratios))
      abline(h = depth.ratios$ratio, lty = 2)
      axis(labels = as.character(depth.ratios$CNt), side = 4, line = 0, las = 1,
           at = depth.ratios$ratio)
      mtext(text = "Copy number", side = 4, line = 2, cex = par("cex.lab")*par("cex"))
   }
   abline(v = max.end, lty = 1)
   axis(labels = chrs, at = coords.names, side = 1, cex.axis = 1)
}
