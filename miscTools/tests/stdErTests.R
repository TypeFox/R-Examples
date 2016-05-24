## testing stdEr() methods
library( miscTools )
data(cars)

# lm()
lmRes <- lm(dist ~ speed, data=cars)
summary( lmRes )
stdEr( lmRes )

# nls()
nlsRes <- nls( dist ~ b0 + b1 * speed^b2, start = c( b0=0, b1=1, b2=1 ),
   data = cars )
summary( nlsRes )
stdEr( nlsRes )
