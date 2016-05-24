test.CAPMList <- function()
{
    
    x <- CAPMList(monthlyReturns, marketIndex = sp500Returns, riskFree = US13wTB)

    .unitTestPath <- BLCOPOptions("unitTestPath")

    expected <- read.csv(file.path(.unitTestPath, "CAPMRes_lm.csv"), row.names = 1 )
    checkEquals( x, expected)
    x <- CAPMList(monthlyReturns, marketIndex = sp500Returns, riskFree = US13wTB, regFunc = "rlm")

    expected <- read.csv(file.path(.unitTestPath, "CAPMRes_rlm.csv"), row.names = 1)
    checkEquals( x, expected)
}