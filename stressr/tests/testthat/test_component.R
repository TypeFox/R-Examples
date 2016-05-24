# unit tests for stressr stress component data
require(testthat)
require(stressr)
require(lattice)

context("Stress Component Data")

# this setup mimics the network call function
# reading our abbreviated test file
x <- readHTMLTable("component-data.xml",
                   header=TRUE,
                   skip.rows=1,
                   colClasses=c("character",rep("numeric",22)),
                   stringsAsFactors=FALSE)

components <- x[[1]]
colnames(components) <- c(
  "DATE",
  "ASSET-BACKED SECURITY SPREAD",
  "BANK BOND SPREAD",
  "COMMERCIAL MBS SPREAD",
  "COMMERCIAL PAPER - T-BILL SPREAD",
  "COMMERCIAL REAL ESTATE SPREAD",
  "CORPORATE BOND SPREAD",
  "COVERED INTEREST SPREAD",
  "CREDIT MARKETS",
  "EQUITY MARKETS",
  "FINANCIAL BETA",
  "FOREIGN EXCHANGE MARKETS",
  "INTERBANK COST OF BORROWING",
  "INTERBANK LIQUIDITY SPREAD",
  "INTERBANK MARKETS",
  "LIQUIDITY SPREAD",
  "REAL ESTATE MARKET",
  "RESIDENTIAL MBS SPREAD",
  "RESIDENTIAL REAL ESTATE SPREAD",
  "SECURITIZATION MARKET",
  "STOCK MARKET CRASHES",
  "TREASURY YIELD CURVE SPREAD",
  "WEIGHTED DOLLAR CRASHES"
)

fedColors <- c(
  "#c1ad97", # "ASSET-BACKED SECURITY SPREAD",
  "#e4aa71", # "BANK BOND SPREAD",
  "#a8927a", # "COMMERCIAL MBS SPREAD",
  "#537f9a", # "COMMERCIAL PAPER - T-BILL SPREAD",
  "#597b6a", # "COMMERCIAL REAL ESTATE SPREAD",
  "#d0e7f7",# "CORPORATE BOND SPREAD",
  "#40657c", # "COVERED INTEREST SPREAD",
  "#5e8aa5", # "CREDIT MARKETS",
  "#dad38d", # "EQUITY MARKETS",
  "#bc8146", # "FINANCIAL BETA",
  "#999999", # "FOREIGN EXCHANGE MARKETS",
  "#d0945a", # "INTERBANK COST OF BORROWING",
  "#fac38e", #"INTERBANK LIQUIDITY SPREAD",
  "#a0c81d", # "INTERBANK MARKETS",
  "#324754", # "LIQUIDITY SPREAD",
  "#2ba4a5", # "REAL ESTATE MARKET",
  "#97846e", # "RESIDENTIAL MBS SPREAD",
  "#6e9180", # "RESIDENTIAL REAL ESTATE SPREAD",
  "#df6968", #" SECURITIZATION MARKET",
  "#dbd48e", # "STOCK MARKET CRASHES",
  "#90b7cf", # "TREASURY YIELD CURVE SPREAD",
  "#999999" # WEIGHTED DOLLAR CRASHES"
)


# convert date posted to date class
components$DATE <- as.Date(components$DATE,format="%m/%d/%Y")
components <- xts(components[,-1],order.by=components[,1])
rv <- list(df=components,
           colors=fedColors,
           columns=1:ncol(components),
           main="Financial Stress Index",
           ylab="")
class(rv) <- "stress"
rv

# reporting tests

test_that("correctly loaded data frame", {
  expect_equal(nrow(rv$df),601)
  expect_equal(ncol(rv$df),22)
  expect_equal(colnames(rv$df)[2],c("BANK BOND SPREAD"))
  expect_equal(as.numeric(last(rv$df[,2])),0.99)
  expect_true("xts" %in% class(rv$df))
})

test_that("correctly loaded list", {
  expect_equal(names(rv),c("df","colors","columns","main","ylab"))
  expect_equal(length(rv$colors),22)
  expect_equal(rv$main,"Financial Stress Index")
  expect_equal(rv$ylab,"")
})

test_that("performs xyplot",{
  expect_equal(class(xyplot(rv)),"trellis")
})

test_that("performs stress line chart",{
  expect_equal(class(stressLineChart(rv)),"trellis")
  expect_equal(class(stressLineChart(rv,range="2013")),"trellis")
  expect_equal(class(stressLineChart(rv,range="2012/2013")),"trellis")
})

test_that("performs stress area chart",{
  expect_equal(class(stressAreaChart(rv)),"trellis")
  expect_equal(class(stressAreaChart(rv,range="2013")),"trellis")
  expect_equal(class(stressAreaChart(rv,range="2012/2013")),"trellis")
})

test_that("performs component summary",{
  cs = getComponentSummary(rv)
  expect_equal(class(cs),"stress")
  expect_equal(length(cs$colors),6)
  expect_equal(length(cs$columns),6)
  expect_equal(cs$columns,c(11,8,14,9,16,19))
  expect_equal(class(xyplot(cs)),"trellis")
  expect_equal(class(stressLineChart(cs,range="2013")),"trellis")
  expect_equal(class(stressAreaChart(cs,range="2012/2013")),"trellis")
})
  
test_that("performs equity markets",{
  cs = getEquityMarkets(rv)
  expect_equal(class(cs),"stress")
  expect_equal(length(cs$colors),1)
  expect_equal(length(cs$columns),1)
  expect_equal(cs$columns,c(9))
  expect_equal(cs$main,"Equity Markets")
  expect_equal(class(xyplot(cs)),"trellis")
  expect_equal(class(stressLineChart(cs,range="2013")),"trellis")
  expect_equal(class(stressAreaChart(cs,range="2012/2013")),"trellis")
})

test_that("performs funding markets",{
  cs = getFundingMarkets(rv)
  expect_equal(class(cs),"stress")
  expect_equal(length(cs$colors),4)
  expect_equal(length(cs$columns),4)
  expect_equal(cs$columns,c(10,12,2,13))
  expect_equal(cs$main,"Funding Markets")
  expect_equal(class(xyplot(cs)),"trellis")
  expect_equal(class(stressLineChart(cs,range="2013")),"trellis")
  expect_equal(class(stressAreaChart(cs,range="2012/2013")),"trellis")
})

test_that("performs credit markets",{
  cs = getCreditMarkets(rv)
  expect_equal(class(cs),"stress")
  expect_equal(length(cs$colors),5)
  expect_equal(length(cs$columns),5)
  expect_equal(cs$columns,c(15,7,4,21,6))
  expect_equal(cs$main,"Credit Markets")
  expect_equal(class(xyplot(cs)),"trellis")
  expect_equal(class(stressLineChart(cs,range="2013")),"trellis")
  expect_equal(class(stressAreaChart(cs,range="2012/2013")),"trellis")
})

test_that("performs foreign exchange markets",{
  cs = getForeignExchangeMarkets(rv)
  expect_equal(class(cs),"stress")
  expect_equal(length(cs$colors),1)
  expect_equal(length(cs$columns),1)
  expect_equal(cs$columns,c(22))
  expect_equal(cs$main,"Foreign Exchange Markets")
  expect_equal(class(xyplot(cs)),"trellis")
  expect_equal(class(stressLineChart(cs,range="2013")),"trellis")
  expect_equal(class(stressAreaChart(cs,range="2012/2013")),"trellis")
})

test_that("performs real estate markets",{
  cs = getRealEstateMarkets(rv)
  expect_equal(class(cs),"stress")
  expect_equal(length(cs$colors),2)
  expect_equal(length(cs$columns),2)
  expect_equal(cs$columns,c(5,18))
  expect_equal(cs$main,"Real Estate Markets")
  expect_equal(class(xyplot(cs)),"trellis")
  expect_equal(class(stressLineChart(cs,range="2013")),"trellis")
  expect_equal(class(stressAreaChart(cs,range="2012/2013")),"trellis")
})

test_that("performs securitization markets",{
  cs = getSecuritizationMarkets(rv)
  expect_equal(class(cs),"stress")
  expect_equal(length(cs$colors),3)
  expect_equal(length(cs$columns),3)
  expect_equal(cs$columns,c(17,3,1))
  expect_equal(cs$main,"Securitization Markets")
  expect_equal(class(xyplot(cs)),"trellis")
  expect_equal(class(stressLineChart(cs,range="2013")),"trellis")
  expect_equal(class(stressAreaChart(cs,range="2012/2013")),"trellis")
})
