IC <-
  function(length, width, thickness) {
    lw <- pmin(length, width)
    ic <- lw / thickness
    ic
  }
