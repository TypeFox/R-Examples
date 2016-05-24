## this version is based on latticeExtra:::c.trellis
## it doesn't work yet with conditioning
print26.AEdotplot <- function(x, ...,
                            main=attr(x, "main"),
                            sub=attr(x,"sub"),
                            ae.key=attr(x, "ae.key"),
                            panel.widths=attr(x,"panel.widths"),
                            AEtable=TRUE) {

  x.in <- x
  if (!AEtable && missing(panel.widths)) panel.widths <- c(.7, .3, 0)
  title.centers <- ifelse (panel.widths[3] > 0, .67, .90)

  ## title.adjust <- function(x, just) {
  ##   if (is.null(x) ||
  ##       is.na(x) ||
  ##       (is.logical(x) && (x == FALSE)) ||
  ##       (!is.list(x) && nchar(x) == 0) ||
  ##       (is.list(x) && nchar(x[[1]]) == 0)) {
  ##     x <- NULL
  ##     x.blank <- x
  ##   }
  ##   else {
  ##     if (!is.list(x)) x <- list(x, cex=1)
  ##     if (!any(names(x) == "label")) {
  ##       empty <- which(names(x) == "")
  ##       if (length(empty) > 0) names(x)[empty[1]] <- "label"
  ##     }
  ##     x.blank <- x
  ##     x.blank[["label"]] <- " "
  ##     if (length(grep("\n", x[["label"]])) > 0)
  ##       x.blank[["label"]] <-
  ##         paste(" ", rep("\n", length(gregexpr("\n", x[["label"]])[[1]])),
  ##               collapse="", sep="")
  ##     x$just <- just
  ##   }
  ##   list(x, x.blank)
  ## }

  ## main.both <- title.adjust(main, just=title.centers)
  ## x$left.plot$main    <- main.both[[2]]
  ## x$right.plot$main   <- main.both[[1]]
  ## x$text.plot$main    <- main.both[[2]]

  ## sub.both <- title.adjust(sub, just=title.centers)
  ## x$left.plot$sub     <- sub.both[[2]]
  ## x$right.plot$sub    <- sub.both[[1]]
  ## x$text.plot$sub     <- sub.both[[2]]

  ## key.blank <- ae.key
  ## key.blank$points$col <- 0
  ## key.blank$text[[1]] <- c(" ", " ")

  ## x$left.plot$legend  <- list(bottom=list(
  ##                               fun="draw.key",
  ##                               args=list(key=ae.key)))
  ## x$right.plot$legend <- list(bottom=list(
  ##                               fun="draw.key",
  ##                               args=list(key=key.blank)))
  ## x$text.plot$legend  <- list(bottom=list(
  ##                               fun="draw.key",
  ##                               args=list(key=key.blank)))

  ## pw <- cumsum(c(0, panel.widths))
  ## pos1 <- c(pw[1], 0, pw[2], 1)
  ## pos2 <- c(pw[2], 0, pw[3], 1)
  ## pos3 <- c(pw[3], 0, pw[4], 1)


if (AEtable || panel.widths[3]==0) {
  ALL3 <- cbind(x$left.plot, x$right.plot, x$text.plot)
  ALL3 <- update(ALL3, scales=list(x=list(relation="free"), alternating=FALSE), between=list(x=.5),
                 main=main, sub=sub)
  ALL3$x.limits <- sapply(x, `[[`, "x.limits")
  ALL3$legend  <- list(bottom=list(fun="draw.key", args=list(key=ae.key)))
  panel.widths.adj <- attr(x,"panel.widths")[1]/3
  panel.widths <- panel.widths + c(-1, .5, .5)*panel.widths.adj
  ALL3 <- resizePanels(ALL3, w=panel.widths)
  ALL3$x.scales$at <- as.list(ALL3$x.scales$at)
  ALL3$x.scales$at[[2]] <- x$right.plot$x.scales$at
  ALL3$x.scales$labels <- list(FALSE, FALSE, FALSE)
  ALL3$x.scales$labels[[2]] <- x$right.plot$x.scales$labels
  print(ALL3)
}
  else {
  ALL2 <- cbind(x$left.plot, x$right.plot)
  ALL2 <- update(ALL2, scales=list(x=list(relation="free"), alternating=FALSE), between=list(x=.5),
                 main=main, sub=sub)
  ALL2$x.limits <- sapply(x, `[[`, "x.limits")[1:2]
  ALL2$legend  <- list(bottom=list(fun="draw.key", args=list(key=ae.key)))
  panel.widths.adj <- attr(x,"panel.widths")[1]/3
  panel.widths <- panel.widths + c(-1, .5, .5)*panel.widths.adj
  ALL2 <- resizePanels(ALL2, w=panel.widths[1:2])
  ALL2$x.scales$at <- as.list(ALL2$x.scales$at)
  ALL2$x.scales$at[[2]] <- x$right.plot$x.scales$at
  ALL2$x.scales$labels <- list(FALSE, FALSE)
  ALL2$x.scales$labels[[2]] <- x$right.plot$x.scales$labels
  print(ALL2)
}

  ## if ((pos1[3] - pos1[1]) > 0)
  ##   print(       x$left.plot                                 ,
  ##         position=pos1, more=TRUE)

  ## if ((pos3[3] - pos3[1]) > 0)
  ##   print(update(x$text.plot, scales=list(y=list(draw=FALSE)), strip.left=FALSE),
  ##         position=pos3, more=TRUE)

  ## ## print right.plot (in the middle) last because it holds main, sub, and key
  ## if ((pos2[3] - pos2[1]) > 0)
  ##   print(update(x$right.plot, scales=list(y=list(draw=FALSE)), strip.left=FALSE),
  ##         position=pos2, more=TRUE)

  ## ##  lattice.setStatus() ## needed because all three panels are set to TRUE
  ## ## lattice:::lattice.setStatus(print.more = FALSE)
  ## lattice.lattice.setStatus(print.more = FALSE)

  invisible(x.in)
}
