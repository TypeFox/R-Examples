 require("stats")
 require("tfplot")

 Sys.info()


tmp <- tempfile()

z <- ts(matrix(100 + rnorm(200),100,2), start=c(1991,1), frequency=4)
tsWrite(z, file=tmp)
zz <- tsScan(tmp, nseries=2)

file.remove(tmp)

cat("max difference ", max(abs(z - zz)) )
if (max(abs(z - zz)) > 1e-10)   stop("file write and read comparison failed.")

####  tfL  ####  

if ( !all(1 == (ts(1:5) - tfL(ts(1:5)))))
       stop("default test of tfL for ts failed.")

if ( !all(1 == (as.ts(1:5) - tfL((1:5)))))
       stop("default test of tfL for non-ts vector failed.")

if ( !all(2 == (ts(1:5) - tfL(ts(1:5), p= 2))))
       stop("2 period lag test of tfL failed.")

z <- ts(1:10, start=c(1992,1), frequency=4)
if ( !all(1 == (z - tfL(z)))) stop("frequency=4 test of tfL failed.")

z <- ts(matrix(1:10,5,2), start=c(1992,1), frequency=4)
seriesNames(z) <- c("One", "Two")
if ( !all(1 == (z - tfL(z)))) stop("matrix test of tfL failed.")

####  annualizedGrowth  ####  

fuzz <- 1e-14
if ( !all(fuzz > (100/(1:4) - annualizedGrowth(ts(1:5)))))
       stop("default test of annualizedGrowth for ts failed.")

#if ( !all(fuzz >  (100/as.ts(1:4) - annualizedGrowth((1:5)))))
#       stop("default test of annualizedGrowth for non-ts vector failed.")

z <- ts(1:5, start=c(1992,1), frequency=4)
if ( !all(fuzz >  (100*((2:5 / 1:4)^4 -1) - annualizedGrowth(z)))) stop("frequency=4 test of annualizedGrowth failed.")

zz <- matrix(1:10,5,2)
z <- ts(zz, start=c(1992,1), frequency=4)
seriesNames(z) <- c("One", "Two")
if ( !all(fuzz >  (100*((zz[2:5,] / zz[1:4,])^4 -1) - annualizedGrowth(z)))) stop("matrix test of annualizedGrowth failed.")
