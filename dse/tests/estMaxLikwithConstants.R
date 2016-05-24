require("dse")
Sys.info()
DSEversion()
 
fuzz <- 1e-6
digits <- 18
all.ok <- TRUE  

test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")


# examples from Olivier Briet


#fix constants:
arma.model <- ARMA(A=c(1, -0.8), B=c(1,0.5), C=NULL, TREND=NULL)
simulated.data <- simulate(arma.model, sampleT=1000, rng=test.rng)
proposed.arma.model <- ARMA(A=c(1, -0.1), B=c(1,0.2), C=NULL, TREND=NULL,
       constants=list(B=array(c(TRUE,TRUE),c(2,1,1))))

estMaxLik.result <-estMaxLik(simulated.data, proposed.arma.model)

print(estMaxLik.result)
print(estMaxLik.result$estimates$like[1], digits=18)
print(TSmodel(estMaxLik.result)$B, digits=18)

good <-  c(1439.50913062493419, 1, 0.2) 
error <- max(abs(good - 
    c(estMaxLik.result$estimates$like[1],TSmodel(estMaxLik.result)$B)))
cat("max. error ", max(error), "\n")
 
if (any(is.na(error)) || any(is.nan(error)) || fuzz < error) all.ok <- FALSE  
 

# trend (p-vector form):
arma.model <- ARMA(A=c(1, -0.8), B=c(1,0.5), C=NULL, TREND=array(c(10),c(1,1)))
simulated.data <- simulate(arma.model, sampleT=1000, rng=test.rng)
estMaxLik.result <- estMaxLik(simulated.data, arma.model)

print(estMaxLik.result)
print(estMaxLik.result$estimates$like[1], digits=18)
print(TSmodel(estMaxLik.result)$TREND, digits=18)

good <-  c(1403.50786987187257 ,10.7031874021547484)
error <- max(abs(good - 
    c(estMaxLik.result$estimates$like[1], TSmodel(estMaxLik.result)$TREND)))
cat("max. error ", max(error), "\n")
 
if (any(is.na(error)) || any(is.nan(error)) || fuzz < error) all.ok <- FALSE  
 

 
# trend (matrix form):

arma.model<- ARMA(A=c(1, -0.8), B=c(1), C=NULL, TREND=matrix(c(1:10),10,1))

#N.B.: The estMaxLik computing time with matrix trend is very long for longer series
simulated.data <- simulate(arma.model, sampleT=10, rng=test.rng)

estMaxLik.result <- estMaxLik(simulated.data, arma.model)

print(estMaxLik.result)
print(estMaxLik.result$estimates$like[1], digits=18)
print(TSmodel(estMaxLik.result)$TREND, digits=18)

good <- c( 1.57484020410352699,
    1.0,  4.12020240399457460, 3.73871024870632374,  3.73115357779961387,
          4.36154038369332131, 6.04943599227131035,  8.02879757311711373,
          8.29303430238073780, 6.69123127611127888, 11.47207950229850937)

error <- max(abs(good -
  c( estMaxLik.result$estimates$like[1],TSmodel(estMaxLik.result)$TREND)))

cat("max. error ", max(error), "\n")
 
if (any(is.na(error)) || any(is.nan(error)) || fuzz < error) all.ok <- FALSE  
 
 
if (! all.ok) stop("some tests FAILED")

 
