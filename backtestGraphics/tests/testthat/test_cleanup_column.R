context("Test for the function 'cleanup_column'")

load("test_cleanup_column.RData")

test_that("cleanup_column function", {
  result.1 <- backtestGraphics:::cleanup_column(x = x,
                                                name.var = "name",   
                                                id.var = "id",
                                                date.var = "date",
                                                nmv.var = "start.nmv",
                                                gmv.var = "gmv",
                                                pnl.var = "pnl.adj",
                                                contract.var = "num.contract.start",
                                                sector.var = "sector",
                                                strategy.var = "strategy",
                                                substrategy.var = "substrategy",
                                                portfolio.var = "port")
  result.2 <- backtestGraphics:::cleanup_column(x = x,
                                                name.var = "wrong.name",   
                                                id.var = "id",
                                                date.var = "date",
                                                nmv.var = "start.nmv",
                                                gmv.var = "gmv",
                                                pnl.var = "pnl.adj",
                                                contract.var = "wrong.contract",
                                                sector.var = "wrong.sector",
                                                strategy.var = "wrong.strategy",
                                                substrategy.var = "wrong.substrategy",
                                                portfolio.var = "wrong.port")
  
  expect_equal(result.1, truth.1, label = "Failed the test for cleaning up the input data set")
  expect_equal(result.2, truth.2, label = "Failed the test for cleaning up the input data set when some columns are missing")
  
  ## Expect an error when both instrument names & ID's are missing"
  
  expect_error(backtestGraphics:::cleanup_column(x = x,
                                                 name.var = "wrong.name",   
                                                 id.var = "bad.id",
                                                 date.var = "date",
                                                 nmv.var = "start.nmv",
                                                 gmv.var = "gmv",
                                                 pnl.var = "pnl.adj",
                                                 contract.var = "wrong.contract",
                                                 sector.var = "wrong.sector",
                                                 strategy.var = "wrong.strategy",
                                                 substrategy.var = "wrong.substrategy",
                                                 portfolio.var = "wrong.port"), 
               label = "Failed to return error when both instrument names & ID's are missing")
  
  ## Expect an error when date is missing
  
  expect_error(backtestGraphics:::cleanup_column(x = x,
                                                 name.var = "wrong.name",   
                                                 id.var = "id",
                                                 date.var = "wrong.date",
                                                 nmv.var = "start.nmv",
                                                 gmv.var = "gmv",
                                                 pnl.var = "pnl.adj",
                                                 contract.var = "wrong.contract",
                                                 sector.var = "wrong.sector",
                                                 strategy.var = "wrong.strategy",
                                                 substrategy.var = "wrong.substrategy",
                                                 portfolio.var = "wrong.port"), 
               label = "Failed to return error when date is missing")
  
  ## Expect an error when both nmv and number of contracts are missing
  
  expect_error(backtestGraphics:::cleanup_column(x = x,
                                                 name.var = "wrong.name",   
                                                 id.var = "id",
                                                 date.var = "date",
                                                 nmv.var = "wrong.nmv",
                                                 gmv.var = "gmv",
                                                 pnl.var = "pnl.adj",
                                                 contract.var = "wrong.contract",
                                                 sector.var = "wrong.sector",
                                                 strategy.var = "wrong.strategy",
                                                 substrategy.var = "wrong.substrategy",
                                                 portfolio.var = "wrong.port"), 
               label = "Failed to return error when both nmv and number of contracts are missing")
})