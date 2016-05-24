 require("dse")
 Sys.info()
 DSEversion()
 
fuzz.small <- 1e-14
fuzz.large <- 1e-10
digits <- 18
all.ok <- TRUE  


cat("dse test 7a...\n")
  d <-20
  true.roots <- c(-1/seq(d),1/seq(d),-seq(d),seq(d)) 
  A <- array(0, c(2,length(true.roots),length(true.roots)))
  A[1,,] <- diag(1,length(true.roots))
  A[2,,] <- diag(-true.roots, length(true.roots))
  # the following relies on roots using by.poly=FALSE
  # Splus needed options(expressions=1024)
   tst <-  sort(roots( ARMA(A=A, B=diag(1,length(true.roots)) ),by.poly=FALSE))
   good <- sort(true.roots)
   error <- max(Mod(good - tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error ) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }

  #  complex roots test
cat("dse test 7b...\n")

  i <- pi*(1:10)/10.1  # this is half circle, but conjs also get generated
  # div by 10.1 instead of 10 prevents a real root with = conj
  true.roots <- complex(real=cos(i), imaginary=sin(i))  # on unit circle
  # scale simplifies sorting
  true.roots <- c(true.roots*(1+.2*1:10), true.roots*(1:10)/10) 
  A <- array(0, c(3,length(true.roots),length(true.roots)))
  A[1,,] <- diag(1,length(true.roots))
  A[2,,] <- diag(-2*Re(true.roots), length(true.roots))
  A[3,,] <- diag(Re(true.roots*Conj(true.roots)), length(true.roots))
  est.roots <- roots( ARMA(A=A, B=diag(1,length(true.roots)) ))
  ec <- 0<=Im(est.roots)

  z  <- est.roots[ ec][order(Mod(est.roots[ ec]))]
  zz <-  true.roots[order(Mod(true.roots))]

  error <- max(Mod(z - zz))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error ) 
     {printTestValue(c(z), digits=18)
      printTestValue(c(zz), digits=18)
      all.ok <- FALSE  
     }

cat("dse test 7c...\n")

  z  <- est.roots[!ec][order(Mod(est.roots[!ec]))]
  zz <- Conj(true.roots)[order(Mod(true.roots))]
  error <- max(Mod(z - zz))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error ) 
     {printTestValue(c(z), digits=18)
      printTestValue(c(zz), digits=18)
      all.ok <- FALSE  
     }


  if (! all.ok) stop("some tests FAILED")

