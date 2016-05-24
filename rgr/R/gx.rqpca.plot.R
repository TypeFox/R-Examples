gx.rqpca.plot <-
function(save, v1 = 1, v2 = 2, rplot = TRUE, qplot = TRUE, rowids = NULL, 
         ifrot = TRUE, main = "", cex.lab = 0.9, cex.main = 0.9, rcex = 1,
         qcex = 0.8, rcol = 1, qcol = 1, ifray = TRUE, if34 = TRUE, ...)
{
     # Function to prepare plots scores on component v2 vs scores on 
     # component v1, loadings of v2 vs v1, or biplots using the saved 
     # output from gx.mva, gx,mva.closed. gx.robmva or gx.robmva.closed.
     # The default is to plot PC-2 vs. PC-1, with data points plotted as
     # default crosses, and the variable names identified.  If other than
     # a cross is required, pch is set to the required value, or the input
     # matrix row number or sample ID may be displayed.
     #
     # Other component combinations may be plotted by setting v1 and v2
     # appropriately.  The plots are always scaled so that both variables
     # and cases (samples) will plot in the space.  Setting ifrot = FALSE
     # causes the unrotated loadings and scores to be plotted from the 
     # gx.rotate saved object.
     #
     # Note: Default is rplot = T & qplot = T & rowids = NULL.
     # The following combinations result in the following plots:
     # rplot = T & qplot = T & rowids = NULL, crosses and variable names
     # rplot = T & qplot = F & rowids = NULL, variable names only
     # rplot = F & qplot = T & rowids = NULL, crosses only
     # rplot = F & qplot = T & rowids = F, sample IDs only
     # rplot = F & qplot = T & rowids = T, input matrix row numbers only
     # rplot = T & qplot = T & rowids = F, sample IDs and variable names
     # rplot = T & qplot = T & rowids = T, input matrix row numbers and variable names
     #
     # By default variable names and sample/case markers are plotted in black
     # and at 100% and 80% font sizes, respectively.  Setting rcol = 2 displays
     # variable names in red, and setting qcol = 4 displays sample/case markers,
     # either current 'pch', row number or sample ID, in blue.  Similarly, changing
     # rcex and qcex provides control of the font sizes.
     #
     frame()
     if(main == "") banner <- paste("PC biplot for",
         deparse(substitute(save)), "\ndata source:", save$input)
     else banner <- main
     nr <- save$nr
     if(is.null(nr)) {
         rload <- save$rload
         rqscore <- save$rqscore
     }
     else {
         if(ifrot) {
             rload <- save$vload
             rqscore <- save$vscore
         }
     else {
         rload <- save$rload
         rqscore <- save$rqscore
         nr <- NULL
         }
     }
     rnames <- save$matnames[[2]]
     qnames <- save$matnames[[1]]
     if(save$proc == "cov") {
         if(is.null(nr)) {
             lv1 <- paste("PC-", as.character(v1), ", ",
                 round(save$econtrib[v1],1), "% of total variability", sep = "")
             lv2 <- paste("PC-", as.character(v2), ", ",
                 round(save$econtrib[v2],1), "% of total variability", sep = "")
             }
         else {
             lv1 <- paste("Varimax Rotated PC-", as.character(v1), ", ",
                 round(save$econtrib[v1],1), "% of total variability", sep = "")
             lv2 <- paste("Varimax Rotated PC-", as.character(v2), ", ",
                 round(save$econtrib[v2],1), "% of total variability", sep = "")
         }
     }
     else {
         if(is.null(nr)) {
             lv1 <- paste("Robust (", save$proc, ") PC-", as.character(v1), ", ",
                 round(save$econtrib[v1],1), "% of total variability", sep = "")
             lv2 <- paste("Robust (", save$proc, ") PC-", as.character(v2), ", ",
                 round(save$econtrib[v2],1), "% of total variability", sep = "")
         }
         else {
             lv1 <- paste("Robust (", save$proc, ") Varimax Rotated PC-",
                 as.character(v1), "\n", round(save$econtrib[v1],1),
                 "% of total variability", sep = "") 
             lv2 <- paste("Robust (", save$proc, ") Varimax Rotated PC-",
                 as.character(v2), "\n", round(save$econtrib[v2],1),
                "% of total variability", sep = "")
         }
     }
     if(!is.null(rowids)) qplot <- TRUE
     if(rplot & !qplot) {
         x1 <- min(rload[, v1])
         x2 <- max(rload[, v1])
         y1 <- min(rload[, v2])
         y2 <- max(rload[, v2])
     }
     if(!rplot & qplot) {
         x1 <- min(rqscore[, v1])
         x2 <- max(rqscore[, v1])
         y1 <- min(rqscore[, v2])
         y2 <- max(rqscore[, v2])
     }
     if(rplot & qplot) {
         x1 <- min(min(rload[, v1]), min(rqscore[, v1]))
         x2 <- max(max(rload[, v1]), max(rqscore[, v1]))
         y1 <- min(min(rload[, v2]), min(rqscore[, v2]))
         y2 <- max(max(rload[, v2]), max(rqscore[, v2]))
     }
     plot(rqscore[, v1], rqscore[, v2], xlab = lv1, ylab = lv2, 
         xlim = c(x1, x2), ylim = c(y1, y2), type = "n", main = banner,
         cex.main = cex.main, cex.lab = cex.lab, ...)
     if((x1 < 0) & (x2 > 0)) abline(v = 0, lty = 2)
     if((y1 < 0) & (y2 > 0)) abline(h = 0, lty = 2)
     if(qplot) {
         if(is.null(rowids))
             points(rqscore[, v1], rqscore[, v2], cex = qcex, col = qcol, ...)
         else if(rowids) text(rqscore[, v1], rqscore[, v2], cex = qcex, col = qcol, ...)
         else text(rqscore[, v1], rqscore[, v2], qnames, cex = qcex, col = qcol, ...)
     }
     if(rplot) {
         if(save$proc == "cov" || !(qplot)) {
             text(rload[, v1], rload[, v2], rnames, cex = rcex, col = rcol, ...)
         }
     else {
         oldpar <- par
         on.exit(par(oldpar))
         par(new = TRUE)
         plot(rload[, v1], rload[, v2], type= "n", xaxt = "n", yaxt = "n",
             xlab = "", ylab = "")
         text(rload[, v1], rload[, v2], rnames, cex = rcex, col = rcol, ...)
         if(if34) {
             axis(3)
             mtext(paste("Robust (", save$proc, ") PC-", v1, " loading", sep = ""),
                 side = 3, line = 2.5, col = rcol)
             axis(4)
             mtext(paste("Robust (", save$proc, ") PC-", v2, " loading", sep = ""),
                 side = 4, line = 2.5, col = rcol)
         } 
     }
         if(ifray) {
             for (i in 1:save$p) {
                 lines(c(0, rload[i, v1]), c(0, rload[i, v2]), col = rcol, lty = 3)
             }
         }
     }
     invisible()
}
