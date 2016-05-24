EMMax.mu.ratiodist <-
function(z,y, mutype=1) {
  ng  <- ncol(z)
  nr  <- nrow(z)
  yw  <- rep(y,ng)
  wgt <- as.vector(z)
  x   <- rep(0:(ng-1),each=nr)

  #start values:
  startf <- 1    # the intrinsic ratio of the two signal responses
  startc <- 0.1  # the signal background parameter
  startd <- -0.1 # the quadratic (nonlinearity of signal response) parameter;
  #startd was 0; with mutype 3 a 0 causes problems and testing shows that -0.1 is not worse than 0
  #lower values:
  lowf <- 0.04 #lower-upper f in 0.04-25 instead of 0.001-1000 gives better convergence but (very slightly) lower LL
  lowc <- 0.001 #with lowc=0, many more convergence errors
  lowd <- -0.25 #was 0; with negative lowd often infinities or convergence errors
            #with lowd=0.001 a few less convergence errors than with lowd=0, but often
            #somewhat (and in some cases much) lower likelihoods
  #upper values:
  uppf <- 25   #see lowf
  uppc <- pi/2 #gives better convergence than uppc=Inf but slightly lower LL
               #pi/4 and pi/6 go further: still better convergence but lower LL
  uppd <- 1 #was Inf, values above 1 hardly change mu values
  switch (mutype,
     { formul <- yw ~ asin(sqrt((c1+x)/(c1+x + c2+f*((ng-1)-x))))                                 # 1 : c1, c2, f
       start=list(c1=startc,c2=startc,f=startf); lower=c(lowc,lowc,lowf); upper=c(uppc,uppc,uppf) },
     { formul <- yw ~ asin(sqrt((c+x)/(c+x + c+f*((ng-1)-x))))                                    # 2 : c, f
       start=list(c=startc,f=startf); lower=c(lowc,lowf); upper=c(uppc,uppf) },
     { formul <- yw ~ asin(sqrt((c1+x+d1*x^2)/(c1+x+d1*x^2 + c2+f*((ng-1)-x)+d2*((ng-1)-x)^2 )))  # 3 : c1, c2, d1, d2, f
       start=list(c1=startc,c2=startc,f=startf,d1=startd,d2=startd); lower=c(lowc,lowc,lowf,lowd,lowd); upper=c(uppc,uppc,uppf,uppd,uppd) },
     { formul <- yw ~ asin(sqrt((c+x+d1*x^2)/(c+x+d1*x^2 + c+f*((ng-1)-x)+d2*((ng-1)-x)^2 )))     # 4 : c, d1, d2, f
       start=list(c=startc,f=startf,d1=startd,d2=startd); lower=c(lowc,lowf,lowd,lowd); upper=c(uppc,uppf,uppd,uppd) },
     { formul <- yw ~ asin(sqrt((c1+x+d*x^2)/(c1+x+d*x^2 + c2+f*((ng-1)-x)+d*((ng-1)-x)^2 )))     # 5 : c1, c2, d, f
       start=list(c1=startc,c2=startc,f=startf,d=startd); lower=c(lowc,lowc,lowf,lowd); upper=c(uppc,uppc,uppf,uppd) },
     { formul <- yw ~ asin(sqrt((c+x+d*x^2)/(c+x+d*x^2 + c+f*((ng-1)-x)+d*((ng-1)-x)^2 )))        # 6 : c, d, f
       start=list(c=startc,f=startf,d=startd); lower=c(lowc,lowf,lowd); upper=c(uppc,uppf,uppd) }
   )

  #first we try nls with the default Gauss-Newton algorithm.
  #This ususally gives a better fit than the "port" algorithm, but it fails more often
  suppressWarnings( {
    success <- tryCatch( {
      res.nls <- nls(formul, start=start, weights=wgt)
     T }, error=function(x) {F} )
    if (!success) {
      #if unsuccessful we try the "port" algorithm which allows to specify lower and upper boundaries
      res.nls <- nls(formul, start=start, weights=wgt, algorithm="port", lower=lower, upper=upper) #werkt
    }
    mu <- predict(res.nls, data.frame(x=0:(ng-1)),type="response")
  })
  mu
}
