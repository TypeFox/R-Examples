###
###
### Create a plot of income inequality in the US
### merging US Census Table F-1
### Income Limits for Each Fifth and Top 5 Percent of Families
###      (All Races)
###
### with Piketty & Saez of "tax units" = families
###    There are slight systematic differences
### between these two sources, but the differences are modest
### relative to the overall image esp. of the top earners
###
###

##
## 1.  Census F-1:
##     Go to
##http://www.census.gov/hhes/www/income/data/historical/families/
## Select "Table F-1 All Races".  This should download something like
##     'F01AR_2011.xls' # put in the working directory
(F01_ <- dir(pattern='^F01AR_2011'))
(F01.xls <- grep('\\.xls$', F01_, value=TRUE))
# Confirm that this is what you want.

library(gdata)
dim(F01 <- read.xls(F01.xls[1], stringsAsFactors=FALSE))

(id <- grep('Dollars', F01[[1]]))
(nyrs <- diff(id)-1)
F01.xls <- F01.xls[1]

dim(F01adj <- read.xls(F01.xls, stringsAsFactors=FALSE,
                   skip=id[2], nrows=nyrs))
str(F01adj)
names(F01adj) <- c('Year', 'Number.thousands',
                   paste('quintile', 1:4, sep=''), 'p95')

library(Ecdat)

iyrs <- (nyrs:1)
F01r <- F01adj[iyrs, 1:7]
str(F01Adj <- asNumericDF(F01adj))
# check
matplot(F01Adj$Year, F01Adj[3:7], log='y', type='l')

#save(F01Adj, file='F01Adj.rda')

##
## 2.  Piketty & Saez 2010 update
##     Go to http://elsa.berkeley.edu/~saez/
##     Download
## (Tables and Figures Updated to 2010 in Excel format, March 2012)
##  with "Income and Wealth Inequality"
## "Income Inequality in the United States, 1913-1998"
## with Thomas Piketty, Quarterly Journal of Economics, 118(1), 2003, 1-39
## It should download as 'TabFig2010.xls';  put in the working directory
##
(PikSaez_ <- dir(pattern='^TabFig20'))
(PikSaez.xls <- grep('\\.xls', PikSaez_, value=TRUE))

dim(PikS <- read.xls(PikSaez.xls, skip=3, sheet='Table A6',
                stringsAsFactors=FALSE))
# 98 28

#  should start with 1913
# which is on row 6 of the file, but read.xls thinks it's row 4
PikS[[1]] # check
# 1913:2010

(PikHeaders <- read.xls(PikSaez.xls, nrows=4, sheet='Table A6',
                       stringsAsFactors=FALSE))
(PikH <- as.character(PikHeaders[2, ]))

# check
(PikSna <- sapply(PikS, function(x)sum(!is.na(x))))
# Columns with name = NA are all NA
# All other columns are officially full
table(PikSna, PikH=="NA")

names(PikS) <- c('year', PikH[-1])

# However, the first 4 rows of P90 and P95 are NA;
# those elements of PikS == '', which asNumericChar converted to NA.
str(PikSaez <- asNumericDF(PikS[, c(1, 16:21)]))

# save(PikSaez, file='PikSaez.rda')

table(sel47. <- (PikSaez$year>1946))
(NAsNeeded <- rep(NA, nyrs-sum(sel47.)))

str(PikS. <- rbind(PikSaez[sel47., ], NAsNeeded))

##
## 3.  GDP
##     Go to
##     http://www.measuringworth.org/usgdp/
##     Select all variables from 1790 to the present,
##     click "Retrieve", then "Download the Results in a Spreadsheet Format"
##
##

(usgdpCSV <- dir(pattern='^USGDP'))

str(usgdp <- read.csv(usgdpCSV, colClasses='character'))

(ngdp <- nrow(usgdp))

# check Year
usgdp$Year[-(3:(ngdp-3))]
# 1790:2011, "Note: ...", "Citation: ..."
Ngdp <- ngdp-2

str(USGDP <- asNumericDF(usgdp[1:Ngdp,]))

names(USGDP) <- c('Year', 'nominalGDP.M', 'realGDP.M',
                  'GDPdeflator', 'PopulationK',
                  'nominalGDPperCap', 'realGDPperCap')

# select years,
plot(GDPdeflator~Year, USGDP, log='y')
abline(h=100)
abline(v=2005) # 2005 = reference year.

# save(usgdp, file='usgdp.rda')

(selGDP <- which(USGDP$Year>1946))

str(GDP <- USGDP[selGDP,])
# check:  goes to 2011., followed by an NA
GDP$Year

n. <- nrow(GDP)
(GDPadj <- GDP$GDPdeflator[n.]/100)
# to convert from 2005 dollars to 2011 dollars
str(GDP. <- with(GDP, data.frame(Year=Year,
                       realGDP.M=GDPadj*realGDP.M,
                        GDP.Deflator=GDPdeflator,
                        PopulationK=PopulationK,
                        realGDPperCap=GDPadj*realGDPperCap)))
# save(GDP., file='GDP.rda')

##
## 4.  Merge
##
# 4.1.  Median
# median approximated by the geometric mean of quintiles 3 and 4.
# This would be correct for a lognormal distribution,
# which is a reasonable approximation for many quantities
# of the distributions of many statistics involving money
names(F01Adj)
F01. <- with(F01Adj, data.frame(Year, Number.thousands, quintile1,
                                quintile2,
                                median=sqrt(quintile2*quintile3),
                                quintile3, quintile4, p95))
str(F01.)

# 4.2.  Adjust PikSaez from 2010 dollars
#       to 2011 dollars used in the Census table F-1
#       and GDP.

F01[1:6,]
# Income in current and 2011 CPI-U-RS adjusted dollars
# PikS. Incomes are expressed in 2010 dollars

table(PSbad <- is.na(PikS.[2]), exclude=NULL)
(PSlast <- max(PikS.[!PSbad, 1]))

(PSlastRow <- which(GDP.[[1]]==PSlast)[1])
(p95.last2 <- GDP.[PSlastRow:nyrs, 'GDP.Deflator'])

(PSinflat <- p95.last2[length(p95.last2)]/p95.last2[1])
# = 1 when both PS and F01 are in the same units., i.e., 2010 dollars

# 4.3.  Combine
str(GDP.)
str(F01.)
str(PikS.)
# all have 64 rows, years starting 1947

str(incomeInequality <- cbind(F01., PSinflat*PikS.[-1], GDP.[-1]))

# 95th percentile:  ratio of IRS to census numbers
incomeInequality$P95IRSvsCensus <- with(incomeInequality, P95/p95)

sum(is.na(incomeInequality$P95IRSvsCensus)) # 1
quantile(incomeInequality$P95IRSvsCensus, na.rm=TRUE)

incomeInequality$personsPerFamily <- with(incomeInequality,
                                      PopulationK/Number.thousands)
incomeInequality$realGDPperFamily <- with(incomeInequality,
                                      personsPerFamily*realGDPperCap)
incomeInequality$mean.median <- with(incomeInequality,
                                  realGDPperFamily/median)
incomeInequality[n.-2:0,]

save(incomeInequality, file='incomeInequality.rda')

