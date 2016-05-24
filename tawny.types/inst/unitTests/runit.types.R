require(futile.any)
require(xts)

get_portfolio <- function()
{
  idx <- as.Date("2011-01-01") + seq(1,20)
  rets <- xts(matrix(rnorm(80), ncol=4), idx)
  colnames(rets) <- c('A','B','C','D')
  TawnyPortfolio(rets, 15)
}

test.window_at.portfolio <- function()
{
  p <- get_portfolio()
  #cat("\ndim(p$returns):",dim(p$returns),"\n")

  p1 <- window_at(p,1)
  checkTrue(nrow(p1$returns) == 15)
  checkTrue(all(p1$returns[1,] == p$returns[1,]))

  p2 <- window_at(p,2)
  checkTrue(nrow(p2$returns) == 15)
  checkTrue(all(p2$returns[1,] == p$returns[2,]))

  p3 <- window_at(p,3)
  checkTrue(nrow(p3$returns) == 15)
  checkTrue(all(p3$returns[1,] == p$returns[3,]))

  p4 <- window_at(p,4)
  checkTrue(nrow(p4$returns) == 15)
  checkTrue(all(p4$returns[1,] == p$returns[4,]))

  p5 <- window_at(p,5)
  checkTrue(nrow(p5$returns) == 15)
  checkTrue(all(p5$returns[1,] == p$returns[5,]))

  p6 <- window_at(p,6)
  checkTrue(nrow(p6$returns) == 15)
  checkTrue(all(p6$returns[1,] == p$returns[6,]))
}

test.start.TawnyPortfolio <- function()
{
  p <- get_portfolio()
  row <- start(p)
  #cat("\nRow:",row,"\n")
  #cat("\nRow:",p$rets[5,],"\n")
  checkTrue(all(row == p$rets[5,]))
}

test.end.TawnyPortfolio <- function()
{
  p <- get_portfolio()
  row <- end(p)
  checkTrue(all(row == p$rets[nrow(p$rets),]))
}

# TODO: Add test to verify rownames
test.rollapply.TawnyPortfolio <- function()
{
  p <- get_portfolio()
  out <- rollapply(p, function(x) colSums(x$returns))
  cat("\nout:\n")
  print(out)
  checkTrue(all(out[1,] == colSums(p$returns[1:15,] )))
  checkTrue(all(out[2,] == colSums(p$returns[2:16,] )))
  checkTrue(all(out[3,] == colSums(p$returns[3:17,] )))
  checkTrue(all(out[4,] == colSums(p$returns[4:18,] )))
  checkTrue(all(out[5,] == colSums(p$returns[5:19,] )))
  checkTrue(all(out[6,] == colSums(p$returns[6:20,] )))
}
