pal <- palettes$palettes; dim(pal) <- c(5,4); pal
palperc <- 100 * row.perc(pal); palperc
palettes$palperc <- as.vector(palperc)
