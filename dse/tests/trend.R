
require("dse")
Sys.info()
DSEversion()
 
fuzz <- 1e-6
digits <- 18
all.ok <- TRUE  

test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")

###################################################

# test with various TREND settings.

###################################################


# Set inputs and outputs for the ARMA model fit and test periods 
arma.fit.input <- c(105.3332, 105.3573, 105.3113, 105.1493, 105.1209, 105.2111, 104.9161,
 		     105.3654, 105.4682, 105.6789, 105.6297, 106.0155, 105.8454, 105.4322,
 		     105.6062, 106.0739, 106.1109, 105.4470, 104.9739, 105.3427, 105.4305,
 		     105.2563, 104.8501, 105.0358, 105.2827, 104.8977)

arma.fit.output <- c(106.0376, 106.0514, 106.0716, 106.0570, 106.0442, 106.0414, 106.0375,
 		      106.0169, 106.0268, 106.0670, 106.1169, 106.1544, 106.1898, 106.2252,
 		      106.2605, 106.2959, 106.3324, 106.3974, 106.3460, 106.2357, 106.1897,
 		      106.1811, 106.1556, 106.1130, 106.0805, 106.0791)

arma.pred.input <- c(104.9916, 104.8207, 104.8936, 104.8767, 104.9435, 104.8885, 104.9217,
 		      104.9029, 104.9508, 105.0065, 105.0557, 105.1982, 105.3392, 105.4007,
 		      105.6212, 105.5979, 105.2410, 105.4832, 105.8735, 105.5944, 105.1063,
 		      104.9809, 105.0821, 104.9362, 105.3037, 105.2322)

arma.pred.output <- c(106.0528, 106.0293, 106.0053, 105.9850, 105.9697, 105.9604, 105.9509,
 		       105.9430, 105.9357, 105.9314, 105.9333, 105.9420, 105.9640, 105.9994,
 		       106.0290, 106.0855, 106.1265, 106.1197, 106.1245, 106.1893, 106.2118,
 		       106.1503, 106.0883, 106.0511, 106.0194, 106.0221)

# Set TSdata object
arma.fit.TSdata <- TSdata(input = arma.fit.input, output = arma.fit.output)

# Fit the model
arma.model.without.trend <- estVARXls(arma.fit.TSdata, max.lag=1, trend=F)

arma.model.with.trend    <- estVARXls(arma.fit.TSdata, max.lag=1, trend=T)

# Apply the model for the test period
arma.pred.TSdata <- TSdata(input = arma.pred.input, output = arma.pred.output[1:2]) 

arma.pred.without.trend <- forecast(TSmodel(arma.model.without.trend), arma.pred.TSdata)

arma.pred.with.trend	 <- forecast(TSmodel(arma.model.with.trend), arma.pred.TSdata)


cat("Test without trend:\n")
z <- sum(arma.pred.without.trend$forecast[[1]][,1])
error <- max(abs((z - 2543.33359644740904)))
cat("   test value:\n")
print(z, digits=18)

if ( fuzz < error) 
     {cat("   error:\n")
      print(error, digits=18)
      all.ok <- FALSE  
     }

cat("Test with trend:\n")
z <- sum(arma.pred.with.trend$forecast[[1]][,1])
error <- max(abs((z - 2544.53598936732942)))
cat("   test value:\n")
print(z, digits=18)

if ( fuzz < error) 
     {cat("   error:\n")
      print(error, digits=18)
      all.ok <- FALSE  
     }

 tfplot(arma.pred.without.trend)
 tfplot(arma.pred.with.trend)
 tfplot(arma.pred.without.trend$forecast[[1]], arma.pred.with.trend$forecast[[1]])


if (! all.ok) stop("some tests FAILED")

