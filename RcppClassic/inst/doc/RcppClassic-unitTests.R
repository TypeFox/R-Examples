### R code from vignette source 'RcppClassic-unitTests.Rnw'

###################################################
### code chunk number 1: RcppClassic-unitTests.Rnw:9-12
###################################################
require( RcppClassic )
prettyVersion <- packageDescription("RcppClassic")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")


###################################################
### code chunk number 2: RcppClassic-unitTests.Rnw:23-29
###################################################
results <- "unitTests-results/RcppClassic-unitTests.txt"
if( file.exists( results ) ){
	writeLines( readLines( results ) )
} else{
	writeLines( "unit test results not available" )
}


