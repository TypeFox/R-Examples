primary.colors <- function(n, steps = 3, no.white = TRUE)
  {
    i <- round(seq(0, 255, length.out = steps))
    if(is.R()) {
        res <- rgb(expand.grid(i, i, i), maxColorValue = 255)
    } else {
        cmat <- expand.grid(i, i, i)
        res <- rgb(cmat[,1], cmat[,2], cmat[,3], maxColorValue = 255)
    }
    if ( no.white ) res <- res[-length(res)]
    if ( missing(n) )
      res
    else
      res[round(seq(1, length(res), length.out = n))]
  }

table.ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1)
  {
    x <- seq(0, 1, length.out = n)
    y <- rep(0, length(x))
    sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
    sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
    y[sill.min:sill.max] <- 1
    base.min <- round((n - 1) * (mid - base / 2)) + 1
    base.max <- round((n - 1) * (mid + base / 2)) + 1
    xi <- base.min:sill.min
    yi <- seq(0, 1, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    xi <- sill.max:base.max
    yi <- seq(1, 0, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    height * y
  }

rgb.tables <- function(n,
                       red = c(0.75, 0.25, 1),
                       green = c(0.5, 0.25, 1),
                       blue = c(0.25, 0.25, 1))
  {
    rr <- do.call("table.ramp", as.list(c(n, red)))
    gr <- do.call("table.ramp", as.list(c(n, green)))
    br <- do.call("table.ramp", as.list(c(n, blue)))
    rgb(rr, gr, br)
  }

matlab.like <- function(n) rgb.tables(n)

matlab.like2 <- function(n)
  rgb.tables(n,
             red = c(0.8, 0.2, 1),
             green = c(0.5, 0.4, 0.8),
             blue = c(0.2, 0.2, 1))

blue2green2red <- matlab.like2

blue2red <- function(n)
  rgb.tables(n,
             red = c(1, 1, 1),
             green = c(0.5, 0, 1),
             blue = c(0, 1, 1))

green2red <- function(n)
  rgb.tables(n,
             red = c(1, 0, 2),
             green = c(0, 0, 2),
             blue = c(0, 0, 0, 0))

blue2green <- function(n)
  rgb.tables(n,
             red = c(0, 0, 0, 0),
             green = c(1, 0, 2),
             blue = c(0, 0, 2))

magenta2green <- function(n)
  rgb.tables(n,
             red = c(0, 0, 2),
             green = c(1, 0, 2),
             blue = c(0, 0, 2))

blue2yellow <- function(n)
  rgb.tables(n,
             red = c(1, 0, 2),
             green = c(1, 0, 2),
             blue = c(0, 0, 2))

cyan2yellow <- function(n)
  rgb.tables(n,
             red = c(1, 0, 2),
             green = c(0.5, 1, 2),
             blue = c(0, 0, 2))

ygobb <- function(n)
{
        rg <- approx(c(0, 0.25, 0.5, 0.75, 1),
                     c(1, 2/3, 2/3, 1/3, 0), n = n)$y
        b <- approx(c(0, 0.25, 0.5, 0.75, 1),
                    c(2/3, 2/3, 1/3, 2/3, 0), n = n)$y
        rgb(rg, rg, b)
}


