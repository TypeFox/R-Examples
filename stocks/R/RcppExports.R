diffs <- function(x, lag = 1) {
  .Call('stocks_diffs_c', PACKAGE = 'stocks', x, lag)
}

pdiffs <- function(x, lag = 1, percent = FALSE) {
  
  # Check that percent is a logical
  if (!is.logical(percent)) {
    stop("For percent input, please enter TRUE or FALSE")
  }

  # Call pdiffs_c with multiplier depending on percent input
  if (percent) {
    .Call('stocks_pdiffs_c', PACKAGE = 'stocks', x, lag, 100)
  } else {
    .Call('stocks_pdiffs_c', PACKAGE = 'stocks', x, lag, 1)
  }

}

pchanges <- function(x, lag = 1, percent = FALSE) {
  
  # Check that percent is a logical
  if (!is.logical(percent)) {
    stop("For percent input, please enter TRUE or FALSE")
  }
  
  # Call pchanges_c with multiplier depending on percent input
  if (percent) {
    .Call('stocks_pchanges_c', PACKAGE = 'stocks', x, lag, 100)
  } else {
    .Call('stocks_pchanges_c', PACKAGE = 'stocks', x, lag, 1)
  }
  
}

ratios <- function(x) {
  .Call('stocks_ratios_c', PACKAGE = 'stocks', x)
}

convert.rate <- function(rate, days.in = 1, days.out = 1) {
  ((rate + 1)^(days.out/days.in)) - 1
}

daily.yearly <- function(daily.gain, years = 1) {
  (1 + daily.gain)^(251*years) - 1
}

yearly.daily <- function(total.gain, years = 1) {
  (total.gain + 1)^(1/(251*years)) - 1
}

balances <- function(ratios, initial = 10000) {
  c(initial, initial * cumprod(ratios))
}

final.balance <- function(ratios, initial = 10000) {
  initial * prod(ratios)
}

prices.rate <- function(prices, xday.rate = NULL) {
  
  # Get the overall growth rate
  prices.length <- length(prices)
  rate1 <- prices[prices.length] / prices[1] - 1
  
  # Convert to x-day growth rate if xday.rate is specified
  if (! is.null(xday.rate) && xday.rate != prices.length - 1) {
    rate1 <- convert.rate(rate = rate1, days.in = prices.length - 1, days.out = xday.rate)
  }
  
  # Return the rate
  return(rate1)
  
}

gains.rate <- function(gains, xday.rate = NULL) {
  
  # Get the overall growth rate
  gains.length <- length(gains)
  rate1 <- prod(gains + 1) - 1
  
  # Convert to x-day growth rate if xday.rate is specified
  if (! is.null(xday.rate) && ! (xday.rate == gains.length)) {
    rate1 <- convert.rate(rate = rate1, days.in = gains.length, days.out = xday.rate)
  }
  
  # Return the rate
  return(rate1)
  
}

positives <- function(x, include.zero = FALSE) {
  
  # Check that include.zero is a logical
  if (!is.logical(include.zero)) {
    stop("For include.zero input, please enter TRUE or FALSE")
  }
  
  # Get positive values
  if (include.zero) {
    x[which(x >= 0)]
  } else {
    x[which(x > 0)]
  }
  
}

negatives <- function(x, include.zero = FALSE) {
  
  # Check that include.zero is a logical
  if (!is.logical(include.zero)) {
    stop("For include.zero input, please enter TRUE or FALSE")
  }
  
  # Get negative values
  if (include.zero) {
    x[which(x <- 0)]
  } else {
    x[which(x < 0)]
  }
  
}

nonpositives <- function(x) {
  x[which(x <= 0)]
}

nonnegatives <- function(x) {
  x[which(x >= 0)]
}

mdd <- function(prices = NULL, gains = NULL, highs = NULL, lows = NULL, indices = FALSE) {
  
  # Check that indices is a logical
  if (!is.logical(indices)) {
    stop("For indices input, please enter TRUE or FALSE")
  }
  
  # If gains specified rather than prices, convert to prices
  if (!is.null(gains)) {
    prices <- balances(ratios = gains + 1, initial = 1)
  }
  
  # Call C++ function depending on indices and whether prices or highs and lows specified
  if (!is.null(prices)) {
    if (!is.null(highs) | !is.null(lows)) {
      stop("Please input prices OR highs and lows, but not both. If both are available, use prices.")
    }
    if (indices) {
      mdd.out <- .Call('stocks_mdd_p_c2', PACKAGE = 'stocks', prices)
      names(mdd.out) <- c("mdd", "start.index", "end.index")
    } else {
      mdd.out <- .Call('stocks_mdd_p_c1', PACKAGE = 'stocks', prices)
    }
  } else {
    if (!is.null(prices)) {
      stop("Please input prices OR highs and lows, but not both. If both are available, use prices.")
    }
    if (indices) {
      mdd.out <- .Call('stocks_mdd_hl_c2', PACKAGE = 'stocks', highs, lows)
      names(mdd.out) <- c("mdd", "start.index", "end.index")
    } else {
      mdd.out <- .Call('stocks_mdd_hl_c1', PACKAGE = 'stocks', highs, lows)
    }
  }
    
  # Return mdd.out
  return(mdd.out)  
  
}

sharpe <- function(gains = NULL, prices = NULL, rf = 0) {
  
  # Check that either gains or prices is specified
  if (is.null(gains) & is.null(prices)) {
    stop("Please enter a gains vector or a prices vector")
  }
  
  # Convert from prices to gains if necessary
  if (!is.null(prices)) {
    gains <- pchanges(prices)
  }
  
  # Calculate Sharpe ratio
  (mean(gains) - rf) / sd(gains)
  
}

sortino <- function(gains = NULL, prices = NULL, rf = 0) {
  
  # Check that either gains or prices is specified
  if (is.null(gains) & is.null(prices)) {
    stop("Please enter a gains vector or a prices vector")
  }
  
  # Convert from prices to gains if necessary
  if (!is.null(prices)) {
    gains <- pchanges(prices)
  }
  
  # Calculate Sortino ratio
  (mean(gains) - rf) / sd(negatives(gains))
  
}

rrr <- function(prices = NULL, gains = NULL) {
  
  # Check that either prices or gains is specified
  if (is.null(prices) & is.null(gains)) {
    stop("Please enter a prices vector or a gains vector")
  }
  
  # Calculate overall growth rate
  if (!is.null(prices)) {
    ret <- prices.rate(prices)
    max.dd <- mdd(prices = prices)
  } else {
    ret <- gains.rate(gains)
    max.dd <- mdd(gains = gains)
  }
  
  # Calculate risk-return ratio
  ret / max.dd
  
}