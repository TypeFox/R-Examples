test.simSummary <-
structure(function()
{

  ## Create simple input from a rather silly simulation - fixed to ease testing
  simFun <- function(x)
  {
    ret <- list()
    ret$s <- 1 + x
    ret$v <- 1:5 + x
    ret$m <- matrix(1:25, nrow=5, ncol=5) + x
    ret$a <- array(1:(4*3*2), dim=c(4, 3, 2)) + x
    ret
  }
  sim <- list()
  sim$sim1 <- simFun(x=-1)
  sim$sim2 <- simFun(x=0)
  sim$sim3 <- simFun(x=1)

  ## ... with NA
  simNA <- sim
  simNA$sim3$s          <- NA
  simNA$sim3$v[5]       <- NA
  simNA$sim3$m[5, 5]    <- NA
  simNA$sim3$a[4, 3, 2] <- NA

  ## Test input - 'x' must be a list
  (checkException(simSummary(x=1)))
  (checkException(simSummary(x=data.frame(1:3, 1:3))))

  ## Test input - "all elements of 'x' must be a list"
  tmp <- sim; tmp$sim4 <- 1
  (checkException(simSummary(x=tmp)))

  ## Test input - "elements of inner lists must be either numeric vector, matrix, or array"
  tmp <- sim; tmp$sim1$l <- list(1:10)
  (checkException(simSummary(x=tmp)))
  tmp <- sim; tmp$sim1$d <- data.frame(a=1:10, b=1:10)
  (checkException(simSummary(x=tmp)))
  tmp <- sim; tmp$sim1$s <- "a"
  (checkException(simSummary(x=tmp)))

  ## Test input - all elements of 'x' must have exactly the same structure
  tmp <- sim; tmp$sim1$s <- 1:3
  (checkException(simSummary(x=tmp)))

  ## Test input - "summary functions must return a single (scalar) value"
  (checkException(simSummary(x=tmp, FUN=range)))

  ## Test summary values - length
  ret <- simSummary(x=sim, FUN="length")
  (checkEquals(ret$length$s, 3, checkNames=FALSE))
  (checkEquals(ret$length$v, rep(3, 5), checkNames=FALSE))
  tmp <- ret$length$m; tmp[] <- 3 # to avoid attributes issues
  (checkEquals(ret$length$m, tmp, checkNames=FALSE))
  tmp <- ret$length$a; tmp[] <- 3 # to avoid attributes issues
  (checkEquals(ret$length$a, tmp, checkNames=FALSE))
  
  ## Test summary values - nobs
  ret <- simSummary(x=sim, FUN="nobs")
  (checkEquals(ret$nobs$s, 3, checkNames=FALSE))
  (checkEquals(ret$nobs$v, rep(3, 5), checkNames=FALSE))
  tmp <- ret$nobs$m; tmp[] <- 3 # to avoid attributes issues
  (checkEquals(ret$nobs$m, tmp, checkNames=FALSE))
  tmp <- ret$nobs$a; tmp[] <- 3 # to avoid attributes issues
  (checkEquals(ret$nobs$a, tmp, checkNames=FALSE))

  ## Test summary values - nobs (in real action)
  ret <- simSummary(x=simNA, FUN="nobs")
  (checkEquals(ret$nobs$s, 2, checkNames=FALSE))
  (checkEquals(ret$nobs$v, c(rep(3, 4), 2), checkNames=FALSE))
  tmp <- ret$nobs$m; tmp[] <- 3; tmp[5, 5] <- 2 # to avoid attributes issues
  (checkEquals(ret$nobs$m, tmp, checkNames=FALSE))
  tmp <- ret$nobs$a; tmp[] <- 3; tmp[4, 3, 2] <- 2 # to avoid attributes issues
  (checkEquals(ret$nobs$a, tmp, checkNames=FALSE))

  ## Test summary values - mean
  ret <- simSummary(x=sim, FUN="mean")
  (checkEquals(ret$mean$s, 1, checkNames=FALSE))
  (checkEquals(ret$mean$v, 1:5, checkNames=FALSE))
  tmp <- ret$mean$m; tmp[] <- 1:25 # to avoid attributes issues
  (checkEquals(ret$mean$m, tmp, checkNames=FALSE))
  tmp <- ret$mean$a; tmp[] <- 1:24 # to avoid attributes issues
  (checkEquals(ret$mean$a, tmp, checkNames=FALSE))

  ## Test summary values - mean & NA
  ret <- simSummary(x=simNA, FUN="mean")
  (checkEquals(ret$mean$s, as.numeric(NA), checkNames=FALSE))
  (checkEquals(ret$mean$v, c(1:4, NA), checkNames=FALSE))
  tmp <- ret$mean$m; tmp[] <- c(1:24, NA) # to avoid attributes issues
  (checkEquals(ret$mean$m, tmp, checkNames=FALSE))
  tmp <- ret$mean$a; tmp[] <- c(1:23, NA) # to avoid attributes issues
  (checkEquals(ret$mean$a, tmp, checkNames=FALSE))

  ## Test summary values - mean & NA but with na.rm=TRUE
  ret <- simSummary(x=simNA, FUN="mean", na.rm=TRUE)
  (checkEquals(ret$mean$s, 0.5, checkNames=FALSE))
  (checkEquals(ret$mean$v, c(1:4, 4.5), checkNames=FALSE))
  tmp <- ret$mean$m; tmp[] <- c(1:24, 24.5) # to avoid attributes issues
  (checkEquals(ret$mean$m, tmp, checkNames=FALSE))
  tmp <- ret$mean$a; tmp[] <- c(1:23, 23.5) # to avoid attributes issues
  (checkEquals(ret$mean$a, tmp, checkNames=FALSE))

  ## Test summary values - min
  ret <- simSummary(x=sim, FUN="min")
  (checkEquals(ret$min$s, 0, checkNames=FALSE))
  (checkEquals(ret$min$v, 0:4, checkNames=FALSE))
  tmp <- ret$min$m; tmp[] <- 0:24 # to avoid attributes issues
  (checkEquals(ret$min$m, tmp, checkNames=FALSE))
  tmp <- ret$min$a; tmp[] <- 0:23 # to avoid attributes issues
  (checkEquals(ret$min$a, tmp, checkNames=FALSE))

  ## Test summary values - max
  ret <- simSummary(x=sim, FUN="max")
  (checkEquals(ret$max$s, 2, checkNames=FALSE))
  (checkEquals(ret$max$v, 2:6, checkNames=FALSE))
  tmp <- ret$max$m; tmp[] <- 2:26 # to avoid attributes issues
  (checkEquals(ret$max$m, tmp, checkNames=FALSE))
  tmp <- ret$max$a; tmp[] <- 2:25 # to avoid attributes issues
  (checkEquals(ret$max$a, tmp, checkNames=FALSE))

}, class = c("svTest", "function"))
