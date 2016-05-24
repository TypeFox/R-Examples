## A method called brinson to apply the Brinson model analysis
## Yang Lu Yang.Lu@williams.edu

brinson <- function(x, 
                    date.var = "date",
                    cat.var = "sector",
                    bench.weight = "benchmark",
                    portfolio.weight = "portfolio",
                    ret.var = "return"){
  ## various checks
  ## x must be a data frame.
  stopifnot(is.data.frame(x))

  ## date.var must have length 1.
  stopifnot(length(date.var) == 1)
  
  ## weight.var must have length 1.
  stopifnot(length(cat.var) == 1)

  ## bench.weight must have length 1.
  stopifnot(length(bench.weight) == 1)

  ## portfolio.weight must have length 1.
  stopifnot(length(portfolio.weight) == 1)
  
  ## ret.var must have length 1.
  stopifnot(length(ret.var) == 1)

  dates <- unique(x[[date.var]])
  
  if (length(dates) > 1){
    .fun <- function(i){brinson(x[x[[date.var]] %in% i, ],
                                date.var,
                                cat.var,
                                bench.weight,
                                portfolio.weight,
                                ret.var)}
    
    multiples <- lapply(dates, .fun)
    
    portfolio.multi <- new("brinsonMulti",
                           date.var = as.character(dates),
                           cat.var = cat.var,
                           bench.weight = bench.weight,
                           portfolio.weight = portfolio.weight,
                           ret.var = ret.var,
                           universe = multiples)
    
    ## more slots including sector weights for portfolio, sector
    ## weights for benchmark, ret.port, ret.bench

    ##
    weight.port.mat <- NULL
    for (j in 1:length(dates)){
      weight.port.mat <- cbind(weight.port.mat, multiples[[j]]@weight.port)
    }
    colnames(weight.port.mat) <- as.character(dates)

    ##
    weight.bench.mat <- NULL
    for (j in 1:length(dates)){
      weight.bench.mat <- cbind(weight.bench.mat, multiples[[j]]@weight.bench)
    }
    colnames(weight.bench.mat) <- as.character(dates)

    ##
    ret.port.mat <- NULL
    for (j in 1:length(dates)){
      ret.port.mat <- cbind(ret.port.mat, multiples[[j]]@ret.port)
    }
    colnames(ret.port.mat) <- as.character(dates)

    ##
    ret.bench.mat <- NULL
    for (j in 1:length(dates)){
      ret.bench.mat <- cbind(ret.bench.mat, multiples[[j]]@ret.bench)
    }
    colnames(ret.bench.mat) <- as.character(dates)
    
    portfolio.multi@ret.port <- ret.port.mat
    portfolio.multi@ret.bench <- ret.bench.mat
    portfolio.multi@weight.port <- weight.port.mat
    portfolio.multi@weight.bench <- weight.bench.mat

    brinson.mat <- matrix(NA, nrow = 4, ncol = length(dates))

    ## q4, portfolio returns
    brinson.mat[1, ] <- sapply(1:length(dates),
                               function(i){ret.port.mat[, i] %*% weight.port.mat[, i]})
    ## q3 
    brinson.mat[2, ] <- sapply(1:length(dates),
                               function(i){ret.port.mat[, i] %*% weight.bench.mat[, i]})
    ## q2
    brinson.mat[3, ] <- sapply(1:length(dates),
                               function(i){ret.bench.mat[, i] %*% weight.port.mat[, i]})
    
    ## q1, benchmark returns
    brinson.mat[4, ] <- sapply(1:length(dates),
                               function(i){ret.bench.mat[, i] %*% weight.bench.mat[, i]})
    
    colnames(brinson.mat) <- as.character(dates)
    rownames(brinson.mat) <- c("q4", "q3", "q2", "q1")
  
    portfolio.multi@brinson.mat <- brinson.mat
    
    return(portfolio.multi)
    
  } else {
    ## single period 
    portfolio <- new("brinson",
                     date.var = date.var,
                     cat.var = cat.var,
                     bench.weight = bench.weight,
                     portfolio.weight = portfolio.weight,
                     ret.var = ret.var,
                     universe = x)

    ## portfolio weight by category
    weight.port <- tapply(x[[portfolio.weight]], x[[cat.var]], sum)

    ## benchmark weight by category
    weight.bench <- tapply(x[[bench.weight]], x[[cat.var]], sum)

    ## benchmark returns
    ret.bench <- .cat.ret(x = x,
                          cat.var = cat.var,
                          ret.var = ret.var,
                          var = bench.weight)
    
    ## portfolio returns
    ret.port <- .cat.ret(x = x,
                         cat.var = cat.var,
                         ret.var = ret.var,
                         var = portfolio.weight)

    ## benchmark returns by category
    ret.bench <- ret.bench / weight.bench

    ## portfolio returns by category
    ret.port <- ret.port / weight.port
    
    ret.bench[ret.bench == "NaN"] <- 0
    ret.port[ret.port == "NaN"] <- 0
    
    portfolio@weight.port <- weight.port
    portfolio@weight.bench <- weight.bench
    portfolio@ret.port <- ret.port
    portfolio@ret.bench <- ret.bench

    q4 <- ret.port %*% weight.port ## portfolio returns
    q3 <- ret.port %*% weight.bench
    q2 <- ret.bench %*% weight.port
    q1 <- ret.bench %*% weight.bench ## benchmark returns

    portfolio@q4 <- q4[1, 1]
    portfolio@q3 <- q3[1, 1]
    portfolio@q2 <- q2[1, 1]
    portfolio@q1 <- q1[1, 1]
    
    return(portfolio)
  }
}
