#' Calculates some initial parameters
#' \code{parini} calculates some initial parameters
#' @keywords internal
#' @export
parini <- function(d, x, k, n, guess, lapses, psyfun) {
  calculate_parini <- function(d, x, k, n, guess, lapses, psyfun) {
    ntrials <- unique(d[[n]])
    yq <- d[[k]] / d[[n]]
    if (is.numeric(guess) && is.numeric(lapses)) {
      gue <- guess
      lap  <- lapses
    }
    if (is.logical(guess) && is.logical(lapses)) {
      if (guess && lapses) {
        gue <- min(yq)
        lap  <- 1 - max(yq)
      }
    }
    if (is.logical(guess) && is.numeric(lapses)) {
      lap <- lapses
      if (guess) gue <- min(yq)
      if (!guess) gue <- 0
    }
    if (is.numeric(guess) && is.logical(lapses)) {
      gue <- guess
      if (lapses) lap <- 1 - max(yq)
      if (!lapses) lap <- 0
    }

    ### Transforming y values to be closer to the range (0,1)
    y01 <- (yq - gue) / (1 - gue - lap)
    datp <- data.frame(x = d[[x]], y01)

    ### Replacing 0s and/or 1s by 1 / (2 * n) and 1 - 1 / (2 * n) where n is the number of trials
    datp <- datp %>%
      mutate(y01 = ifelse(y01 == 1, 1 - 1 / (2 * ntrials), y01)) %>%
      mutate(y01 = ifelse(y01 == 0, 1 / (2 * ntrials), y01))

    ### Eliminating probabilities outside (0,1)
    dat <- filter(datp, y01 > 0, y01 <1)

    ### Linear fit
    dat$z <- qnorm(dat$y01)
    coef <- lm(z~x, data = dat)$coefficients

    if (coef[[2]] == 0) { # checking that the slope is not zero
      p1 <- median(dat$x)
      p2 <- (1 - gue - lap) / (max(dat$x)-min(dat$x))
    }
    else {
      p1 <- -coef[[1]] / coef[[2]]
      p2 <- 1 / coef[[2]]
    }

    if (psyfun == 'logistic_fun') p2 <- 1 / p2
    if (psyfun == 'weibull_fun') p2 <- 1 / p2

    if (is.numeric(guess) && is.numeric(lapses)) para <- c(p1, p2)
    if (is.logical(guess) && is.logical(lapses)) {
      if (guess && lapses) para <- c(p1, p2, gue, lap)
      if (!guess && !lapses) para <- c(p1, p2)
    }
    if (is.logical(guess) && is.numeric(lapses)) {
      if (guess) para <- c(p1, p2, gue)
      if (!guess) para <- c(p1, p2)
    }
    if (is.numeric(guess) && is.logical(lapses)) {
      if (lapses) para <- c(p1, p2, lap)
      if (!lapses) para <- c(p1, p2)
    }
    data.frame(paran = paste0('p', seq(1, length(para))), par = para)
  }
  d %>% do(calculate_parini(., x, k, n, guess, lapses, psyfun))
}

