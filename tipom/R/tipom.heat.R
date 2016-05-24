
tipom.heat <-
  function(lengths, widths, thicknesses, ...) {
    lwmin <- pmin(widths, lengths)
    thck <- thicknesses
    ftab <- table(lwmin, thck)

    x <- as.numeric(dimnames(ftab)$lwmin)
    y <- as.numeric(dimnames(ftab)$thck)

    tipom.plot <- image(x, y, ftab, col=gray((32:0)/32))
    tipom.plot
  }
