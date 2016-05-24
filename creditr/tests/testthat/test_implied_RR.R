context("Test implied_RR")

## test case for implied recovery rate using a data frame containing ids, spreads and pds of 
## seven different CDSs

library(creditr)

test_that("test for implied recovery rate", {
  
  ids <- c("CaesarsEntCorp", "ElectroluxAB", "Chorus", "NorskeSkogindustrier", 
           "TokyoElectricPower", "ToysRUs", "Xerox")
  
  pd <- c(0.99998, 0.0827, 0.1915, 0.9128, 0.1830, 0.7813, 0.0880)
  
  spread <- c(12354.53, 99, 243.28, 2785.8889, 250.00, 1737.7289, 105.8)
  
  date <- c(as.Date("2014-04-15"), as.Date("2014-04-22"), as.Date("2014-04-15"), 
            as.Date("2014-04-15"), as.Date("2014-04-15"), as.Date("2014-04-15"), 
            as.Date("2014-04-22"))
  
  endDate <- c(as.Date("2019-06-20"), as.Date("2019-06-20"), as.Date("2019-06-20"), 
               as.Date("2019-06-20"), as.Date("2019-06-20"), as.Date("2019-06-20"), 
               as.Date("2019-06-20"))   
  
  df <- data.frame(ids, pd, spread, date, endDate)
  
  # 3 is the column number for spread and 2 is the column number for pd, 1 for id
  
  result <- implied_RR(x = df, 
                       date.var = "date",
                       spread.var = "spread", 
                       pd.var = "pd", 
                       maturity.var = "endDate")
  
  truth <- c(40, 40, 40, 40, 35, 40, 40)
  
  ## the implied recovery rate is marginally different (less than 0.02%) from the recovery rates
  ## provided by bloomberg, so we round off the results
  
  expect_equal(round(result), truth)
  
})