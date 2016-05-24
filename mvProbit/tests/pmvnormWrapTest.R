library( "mvProbit" )
library( "miscTools" )

# covariance matrix
sigma <- miscTools::symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

######## only upper ##########
upper <- c( -0.3, 0.7, -0.5 )
# Genz + Bretz (default)
pug <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = GenzBretz(), random.seed = 123 )
print( pug )

# Miwa (as function)
pum <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = Miwa, random.seed = 123 )
print( pum )
all.equal( pug, pum, check.attributes = FALSE, tol = 1e-4 )

# Miwa (as object returned from function Miwa())
pum1 <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = Miwa(), random.seed = 123 )
all.equal( pum, pum1 )

# Miwa (as character string)
pum2 <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = "Miwa", random.seed = 123 )
all.equal( pum, pum2 )

# TVPACK
put <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = TVPACK, random.seed = 123 )
print( put )
all.equal( pug, put, check.attributes = FALSE, tol = 1e-4 )
all.equal( pum, put, check.attributes = FALSE, tol = 1e-6 )

# GHK
pughk <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = "ghk", random.seed = 123, nGHK = 1000 )
print( pughk )
all.equal( pug, pughk, tol = 1e-3, check.attributes = FALSE )
all.equal( pum, pughk, tol = 1e-3, check.attributes = FALSE )

# GHK, lower precision
pughk1 <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = "ghk", random.seed = 123, nGHK = 100 )
print( pughk1 )
all.equal( pughk, pughk1, tol = 1e-2, check.attributes = FALSE )
all.equal( pug, pughk1, tol = 1e-2, check.attributes = FALSE )


######## only lower ##########
lower <- c( -0.7, 0.3, -0.9 )
# Genz + Bretz (default)
plg <- mvProbit:::pmvnormWrap( lower = lower, sigma = sigma, 
   algorithm = GenzBretz, random.seed = 123 )
print( plg )

# Miwa
plm <- mvProbit:::pmvnormWrap( lower = lower, sigma = sigma, 
   algorithm = Miwa, random.seed = 123 )
print( plm )
all.equal( plg, plm, tol = 1e-4, check.attributes = FALSE )

# TVPACK
plt <- mvProbit:::pmvnormWrap( lower = lower, sigma = sigma, 
   algorithm = TVPACK, random.seed = 123 )
print( plt )
all.equal( plg, plt, tol =1e-4, check.attributes = FALSE )
all.equal( plm, plt, check.attributes = FALSE, tolerance = 1e-8 )

# GHK
plghk <- mvProbit:::pmvnormWrap( lower = lower, sigma = sigma, 
   algorithm = "GHK", random.seed = 123, nGHK = 1000 )
print( plghk )
all.equal( plg, plghk, tol = 1e-3, check.attributes = FALSE )
all.equal( plm, plghk, tol = 1e-3, check.attributes = FALSE )


######## partly lower, partly upper ##########
lower2 <- c( -Inf, 0.3, -Inf )
upper2 <- c( -0.3, Inf, -0.5 )
# Genz + Bretz (default)
pbg <- mvProbit:::pmvnormWrap( lower = lower2, upper = upper2, sigma = sigma, 
   algorithm = GenzBretz(), random.seed = 123 )
print( pbg )

# Miwa
pbm <- mvProbit:::pmvnormWrap( lower = lower2, upper = upper2, sigma = sigma, 
   algorithm = Miwa, random.seed = 123 )
print( pbm )
all.equal( pbg, pbm, tol = 1e-3, check.attributes = FALSE )

# GHK
pbghk <- mvProbit:::pmvnormWrap( lower = lower2, upper = upper2, sigma = sigma, 
   algorithm = "GHK", random.seed = 123, nGHK = 1000 )
print( pbghk )
all.equal( pbg, pbghk, tol = 1e-3, check.attributes = FALSE )
all.equal( pbm, pbghk, tol = 1e-4, check.attributes = FALSE )

