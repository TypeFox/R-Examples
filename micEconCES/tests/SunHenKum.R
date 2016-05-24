## This script reproduces some of the results of
## Sun, K., Henderson, D.J. and Kumbhakar, S.C. (2011),
## Biases in approximating log production.
## Journal of Applied Econometrics, forthcoming.
## doi: 10.1002/jae.1229
##
## This paper is in fact a replication study of
## Masanjala, W.H. and Papageorgiou, C. (2004), 
## The Solow model with CES technology: nonlinearities and parameter 
## heterogeneity. Journal of Applied Econometrics, 19: 171–201.
## doi: 10.1002/jae.722
##
## The data used in the two above-mentioned papers are actually from
## Durlauf, S.N. and Johnson, P.A. (1995),
## Multiple Regimes and Cross-Country Growth Behavior. 
## Journal of Applied Econometrics, 10: 365–384.

# load the "micEconCES" package
library( "micEconCES" )
options( digits = 3 )

# load data (included in the "AER" package)
data( "GrowthDJ", package = "AER" )

# remove data from oil producing countries
# as this has been done by Masanjala and Papageorgiou (2004)
# and hence, also by Sun, Henderson, and Kumbhakar (2011)
GrowthDJ <- subset( GrowthDJ, oil == "no" )

# calculate "input" variables for the Solow growth model
GrowthDJ$x1 <- 1
GrowthDJ$x2 <- ( GrowthDJ$popgrowth + 5 ) / GrowthDJ$invest

# CES: non-linear least-squares estimation (NLLS)
cesNls <- cesEst( "gdp85", c( "x1", "x2"), data = GrowthDJ )
summary(cesNls)
cesNlsCoef <- coef( cesNls )
delta <- cesNlsCoef[ "delta" ]
rho <- cesNlsCoef[ "rho" ]
print( alpha <- ( delta - 1 ) / delta )
print( sigma <- 1 / ( 1 - rho ) )
cesNlsVar <- vcov(cesNls)
# deltamethod(~-x2/(1-x2), cesNlsCoef, cesNlsVar)
# deltamethod(~1/(1-x3), cesNlsCoef, cesNlsVar)

# Cobb-Douglas: non-linear least-squares estimation (NLLS)
cdNls <- cesEst( "gdp85", c( "x1", "x2"), data = GrowthDJ, rho = 0 )
summary(cdNls)
cdNlsCoef <- coef( cdNls )
delta <- cdNlsCoef[ "delta" ]
print( alpha <- ( delta - 1 ) / delta )

# Cobb-Douglas: estimation with logs
cdLog <- cesEst( "gdp85", c( "x1", "x2"), data = GrowthDJ, rho = 0, multErr = TRUE )
summary(cdLog)
cdLogCoef <- coef( cdLog )
delta <- cdLogCoef[ "delta" ]
print( alpha <- ( delta - 1 ) / delta )

