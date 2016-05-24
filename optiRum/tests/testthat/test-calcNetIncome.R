context("calcNetIncome")

taxyear<-2014

test_that("TaxOwed correct results", {
  incomeTaxable <- c(0,500,30000,32000,1e5,2e5)
  taxBrackets <- fread(system.file("extdata", "annualtaxthresholds.csv", 
                                   package = "optiRum"))[TaxYear==2014&Type=="Income Tax"]
  expect_identical(TaxOwed(incomeTaxable,taxBrackets), c(0, 100, 6000, 6427, 33627, 76127), "tax rates calc")
}) 

test_that("calcNetIncome on a single household", {
  persons <- data.table(householdID                    = c(1L,1L),
                          personID                = c(1L,2L),
                          employedIncome            = c(100000,50000),
                          investmentIncome          = c(2000,500),
                          nonTaxableIncome          = c(0,0),
                          selfEmployedProfits       = c(0,0),
                          salarySacrificePercentage = c(0.05,0.05),
                          taxCode                   = c("",""),
                          studentLoan               = c(FALSE,FALSE),
                          numberOfChildren          = c(1,1) )
  expected_r <- data.table(householdID = c(1L, 1L), # valid for reference data as of 24/03/2015
                           personID = 1:2, 
                           personalAllowance = c(10000, 10000),
                           totalIncome = c(8588.83333333333, 4208.33333333333),
                           netIncome = c(5703.47666666667, 3124.31), 
                           householdNetIncome = c(8827.78666666667, 8827.78666666667), 
                           incomeTax = c(2368.91666666667, 735.583333333333), 
                           class1NI = c(427.606666666667, 348.44), 
                           class4NI = c(0, 0), 
                           childBenefits = c(88.8333333333333, 0), 
                           childBenefitTax = c(88.8333333333333, 0), 
                           studentLoanRepayment = c(0, 0), key=c("householdID","personID"))
  expect_equal(calcNetIncome(persons,financialYear=taxyear),expected_r)
})

test_that("calcNetIncome on multiple households", {
  persons <- data.table(householdID                    = c(1L,1L,2L,2L,3L),
                          personID                = c(1L,2L,1L,2L,1L),
                          employedIncome            = c(45000,14000,220000,190000,140000),
                          investmentIncome          = c(500,0,0,1500,10000),
                          nonTaxableIncome          = c(0,0,3000,0,1500),
                          selfEmployedProfits       = c(0,0,0,5000,2000),
                          salarySacrificePercentage = c(0.05,0.05,0.05,0.05,0.05),
                          taxCode                   = c("","","","",""),
                          studentLoan               = c(FALSE,FALSE,TRUE,FALSE,TRUE),
                          numberOfChildren          = c(1,1,3,3,2))
  expected_r <- data.table(householdID = c(1L, 1L, 2L, 2L, 3L), # valid for reference data as of 24/03/2015
                           personID = c(1L, 2L, 1L, 2L, 1L),
                           personalAllowance = c(10000, 10000, 0, 0, 0), 
                           totalIncome = c(3880.5, 1166.66666666667, 18789.6, 16375, 12939.2166666667), 
                           netIncome = c(2962.72666666667, 1058.22666666667, 9843.635, 10135.9766666667, 7119.46833333333), 
                           householdNetIncome = c(4020.95333333333, 4020.95333333333, 19979.6116666667, 19979.6116666667, 7119.46833333333), 
                           incomeTax = c(577.25, 55, 6681.41666666667, 5668.91666666667, 4235.58333333333), 
                           class1NI = c(340.523333333333, 53.44, 617.606666666667, 570.106666666667, 490.94),
                           class4NI = c(0, 0, 0, 0, 0), 
                           childBenefits = c(88.8333333333333, 0, 206.266666666667, 0, 147.55), 
                           childBenefitTax = c(0, 0, 206.266666666667, 0, 147.55), 
                           studentLoanRepayment = c(0, 0, 1440.675, 0, 945.675), key=c("householdID","personID"))
  expect_equal(calcNetIncome(persons,financialYear=taxyear),expected_r)
})
