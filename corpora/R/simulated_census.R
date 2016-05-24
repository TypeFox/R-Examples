simulated.census <- function (N=502202, p.male=0.55, seed.rng=42) {
  # default N = population of Luxembourg, as estimate on 1 Jan 2010 (from http://en.wikipedia.org/wiki/Luxembourg)

  # save current RNG state and set new seed (for reproducible data set)
  runif(1) # make sure there is a .Random.seed
  save.seed <- .Random.seed
  set.seed(seed.rng, kind="default")

  height2weight <- function (h, h0=170, w0=65, slope=1, curvature=0.01, sd.log=.2) {
    x <- h - h0
    w <- w0 + (slope * (1 + curvature * x)) * x
    scale <- 2 ^ rnorm(length(w), 0, sd.log)
    w * scale
  }

  height2shoe <- function (h, h0=170, s0=42, slope=0.2, sd=1) {
    x <- h - h0
    sizes <- s0 + slope * x + rnorm(length(s0), 0, sd)
    round(2 * sizes) / 2
  }

  N.male <- round(N * p.male)
  h.male <- rnorm(N.male, 180, 10)
  w.male <- height2weight(h.male, h0=180, w0=75)
  s.male <- height2shoe(h.male, h0=180, s0=44.5)

  N.female <- N - N.male
  h.female <- rnorm(N.female, 160, 10)
  w.female <- height2weight(h.female, h0=160, w0=50, slope=0.5, curvature=0.02, sd.log=.15)
  s.female <- height2shoe(h.female, h0=160, s0=38, slope=0.15)

  males <- data.frame(height=h.male, weight=w.male, shoe.size=s.male, sex="m")
  females <- data.frame(height=h.female, weight=w.female, shoe.size=s.female, sex="f")
  FakeCensus <- rbind(males, females)

  FakeCensus <- FakeCensus[sample(1:N), ]
  rownames(FakeCensus) <- 1:N

  # restore RNG state, then return generated data set
  .Random.seed <<- save.seed
  FakeCensus 
}