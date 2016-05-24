# Diana Hall
# 12-28-2015
# raster plotting function
.plot.mm.s<-function (s, whichcells = NULL, beg = min(unlist(s$spikes), na.rm = TRUE), 
                     end = max(unlist(s$spikes), na.rm = TRUE), label.cells = FALSE,
                     show.burst.number=F ,
                     use.names = TRUE, show.bursts = FALSE, main = NULL, ylab = "Unit", 
                     xlab = "Time (s)", for.figure = FALSE, show.episodes, episode.y = -0.01, 
                     ...) 
{
  if (length(whichcells) > 0 && is.numeric(whichcells[1])) {
  }
  else {
    whichcells = .names.to.indexes(names(s$spikes), whichcells, 
                                  allow.na = TRUE)
  }
  if (is.null(main)) {
    main <- basename(s$file)
  }
  N <- length(whichcells)
  ticpercell <- 1/N
  deltay <- ticpercell * 0.8
  yminadd <- ticpercell
  if (show.bursts) 
    spikes <- s$spikes
  else spikes <- s$spikes
  if (for.figure) {
    plot(c(beg, end), c(0, 1), type = "n", yaxt = "n", bty = "n", 
         main = "", xaxt = "n", xaxs = "i", yaxs = "i", xlab = "", 
         ylab = "", ...)
    mtext(main, side = 3, adj = 0, line = 0.5)
  }
  else {
    plot(c(beg, end), c(0, 1), type = "n", bty = "n", yaxt = "n", 
         main = main, xlab = xlab, ylab = ylab, ...)
  }
  ymin <- 0
  have.bursts <- ((length(s$allb) > 0) && show.bursts)
  for (cell in whichcells) {
    ts <- spikes[[cell]]
    n <- length(ts)
    if (n > 0) {
      ys <- numeric(n) + ymin
      segments(ts, ys, ts, ys + deltay, lwd = 0.2)
      if (have.bursts) {
        burst.times <- s$allb[[cell]]
        if (!is.na(burst.times[1])) {
          nbursts <- nrow(burst.times)
          ys <- rep(ymin + deltay/2, nbursts)
          shimmy <- deltay * 0.25
          odd <- (1:nbursts)%%2 == 1
          ys[odd] <- ys[odd] + shimmy
          start.burst <- ts[burst.times[, "beg"]]
          end.burst <- ts[burst.times[, "beg"] + burst.times[, 
                                                             "len"] - 1]
          segments(start.burst, ys, end.burst, ys, col = "red", 
                   lwd = 2)
          if (show.burst.number){
            text(start.burst, rep(ymin + deltay * 1.1, 
                                  nbursts), labels = burst.times[, "len"], 
                 col = "blue")
          }
        }
      }
    }
    ymin <- ymin + yminadd
  }
  if (label.cells) {
    allys <- seq(from = yminadd/2, by = yminadd, length = N)
    if (use.names) {
      labels <- names(spikes)[whichcells]
    }
    else {
      labels <- whichcells
    }
    axis(2, at = allys, labels = labels, las = 1, tick = F)
  }
  if (missing(show.episodes)) {
    show.episodes <- ("episodes" %in% names(s))
  }
  if (show.episodes) {
    segments(s$episodes[, "beg"], episode.y, s$episodes[, 
                                                        "end"], episode.y, col = "purple", xpd = NA)
  }
}

