## point symbols:
## 2 = pch(3), +
## 4 = pch(4), x
## 8 = pch(8), *

tipom.car <-
  function(lengths, widths, thicknesses, ic=FALSE, bubble=FALSE, ...) {
    lwmin <- pmin(widths, lengths)
    thck <- thicknesses
    ftab <- as.data.frame(table(lwmin, thck))
    ftab <- ftab[ftab$Freq>0, ]
    ftab$lwmin <- as.numeric(ftab$lwmin)
    ftab$thck <- as.numeric(ftab$thck)

    if (bubble == TRUE) {
      tipom.plot <- with(ftab,
                         plot(thck, lwmin,
                              ylim=c(0, max(lwmin)+max(lwmin)*0.1),
                              xlim=c(0, max(thck)+max(thck)*0.1),
                              asp=1,
                              xlab="Thickness",
                              cex=Freq*0.2,
                              pch=20,
                              ylab="Length or Width",
                              ...)
                         )
    } else {
      syms <- c(20, 3, 3, rep(4, 4), rep(8, max(ftab$Freq)))
      cexs <- c(0.2, 0.7, 0.7, rep(1, max(ftab$Freq)))
      tipom.plot <- with(ftab,
                         plot(thck, lwmin,
                              ylim=c(0, max(lwmin)+max(lwmin)*0.1),
                              xlim=c(0, max(thck)+max(thck)*0.1),
                              asp=1,
                              xlab="Thickness",
                              pch=syms[Freq],
                              cex=cexs[Freq],
                              ylab="Length or Width",
                              xaxt="n",
                              yaxt="n",
                              frame=FALSE,
                              ...)
                         )
    }
    if (ic == TRUE) {
      slopes <- c(1, 1.5, 2, 2.5, 4, 8)
      for (s in 1:6)
        abline(a=0, b=slopes[s])
    }
    axis(2, yaxp=c(0,max(lwmin),max(lwmin)), pos=0)
    axis(1, xaxp=c(0,max(thck),max(thck)), pos=0)
    tipom.plot
  }
