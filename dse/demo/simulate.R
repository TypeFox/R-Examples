#######  examples from User's Guide section 5

  require("dse")

  cat("Define an ARMA TSmodel object...\n")
  arma <- ARMA(A=array(c(1,.5,.3,0,.2,.1,0,.2,.05,1,.5,.3), c(3,2,2)),
               B=array(c(1,.2,0,.1,0,0,1,.3), c(2,2,2)),
	       C=NULL) 
  arma  # or print(arma)

  cat("   generate simulated data\n")
  data.arma.sim <- simulate(arma)

  cat("   evaluate the model with the simulated data to get a TSestModel object\n")
  arma <- l(arma, data.arma.sim)
  summary(arma)
  roots(arma)
  stability(arma)

  tfplot(data.arma.sim)
  tfplot(arma)

  cat("Define a State Space TSmodel object...\n")
  ss <- SS(F=array(c(.5,.3,.2,.4),c(2,2)),
           G=NULL,
	   H=array(c(1,0,0,1),c(2,2)),
	   K=array(c(.5,.3,.2,.4),c(2,2)))
  ss   # or print(ss)

  cat("   generate simulated data\n")
  data.ss.sim <- simulate(ss)

  cat("   evaluate the model with the simulated data to get a TSestModel object\n")
  ss <- l(ss, data.ss.sim)

  summary(ss)
  roots(ss)
  stability(ss)
  
  tfplot(ss)

  cat("Convert between State Space and ARMA models\n")
  ss.from.arma <- l(toSS(arma), data.arma.sim)
  arma.from.ss <- l(toARMA(ss), data.ss.sim)
  
  summary(ss.from.arma)
  summary(arma.from.ss)
  
  stability(arma)
  stability(ss.from.arma)
