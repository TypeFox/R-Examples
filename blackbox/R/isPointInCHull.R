isPointInCHull <- function(point, testspace001, constraints, exactArithmetic=F, returnValue="logical") { ## is point in convex hull of testspace ?
  ## FR->FR afaire: coder point -> matrice. cf vinny.pdf for syntax...
  ## input point is alsways numeric...
  ## cf "Convex Hull Revisited" in vinny.pdf: a set of constraints (an hrepr) is given directly by the vertices (a vrepr...)
  if( ! is.numeric(point)) {
    message.redef("'point' argument in 'isPointInCHull()' should be of class 'numeric'; however")
    message.redef(paste(" however, it is of class", class(point), "with value:"))
    message.redef(paste(signif(point, 6)))
    stop.redef() #safe side, but we should not have reached this point
  }
  if(is.null(names(point))) {
    stop.redef("(!) is.null(names(point)) in isPointInCHull().")
  }
  if( ! missing(constraints)) {
    if ( ! is.numeric(constraints$a[1]) &&  ! is.character(constraints$a[1])) {
      stop.redef("(!) From isPointInCHull(): 'constraints' given but not appropriately numeric or character.")
    }
    if(is.null(colnames(constraints$a))) {
      stop.redef("(!) is.null(colnames(constraints$a)) in isPointInCHull().")
    }
  } else {
    if (missing(testspace001)) {
      stop.redef("(!) From isPointInCHull(): neither 'constraints' nor 'testspace001' given.")
    } else if ( ! is.numeric(testspace001[1]) && ! is.character(testspace001[1])) {
      stop.redef("(!) From isPointInCHull(): 'testspace001' given but not appropriately numeric or character.")
      if(is.null(colnames(testspace001))) {
        stop.redef("(!) is.null(colnames(testspace001)) in isPointInCHull().")
      }
    }
  }
  # Ideally we should first check whether col numbers are identical,
  # and this would make interactive checks of individuals points painful.
  #if(any(names(point) != colnames(testspace001)[-(1:3)])) {
  #    stop.redef("(!) names(point) != colnames(testspace001)[-(1:3)] in isPointInCHull().")
  #}
  if (missing(constraints)) {
    if (returnValue=="vector") stop.redef("(!) From isPointInCHull: option returnValue=vector incompatible with missing constraints")
    if ( ! exactArithmetic) { ## then first try floating point
      if (is.character(testspace001[1])) testspace001 <- q2d(testspace001) ## stupidly inefficient case...
      if (is.character(point)) point <- q2d(point)
      hrep <- rbind(testspace001, c(0, 1, 1, -point))
      hpt <- c(-1, point)
      out <- try(lpcdd(hrep, hpt, minimize=F)) ## out$optimal should be negative if point is within
      if( (! inherits(out,"try-error")) && ( ! is.null(out$optimal.value)) ) {
        if (returnValue=="logical") rv <- !as.logical(as.numeric(out$optimal.value))
        if (returnValue=="max") rv <- as.numeric(out$optimal.value)
      }
    } ## ELSE : either exactArithmetic was required or floating point failed to return a valid result
    ## exactArithmetic is more accurate but slower
    if (is.numeric(testspace001[1])) testspace001 <- d2q(testspace001)
    if (is.numeric(point)) point <- d2q(point)
    hrep <- rbind(testspace001, c("0", "1", "1", qneg(point)))
    hpt <- c("-1", point)
    out <- try(lpcdd(hrep, hpt, minimize=F))
    if(inherits(out,"try-error")) { ## failure with exact arithmetic (a likely trivial cause: incompatible dimensions)
      if ( ! exactArithmetic) message.redef("call to lpcdd(...) (floating point), and")
      message.redef("call to lpcdd(...) (rational) failed at point:")
      message.redef(q2d(point))
      save(hrep, point, file=paste("lpcddData_", .blackbox.data$options$jobSampleNbr, sep=""))
      message.redef(paste("'hrep' and 'point' are saved in file lpcddData_", .blackbox.data$options$jobSampleNbr, sep=""))
    } else if(is.null(out$optimal.value)) {
      if ( ! exactArithmetic) message.redef("call to lpcdd(...) (floating point) failed, and")
      message.redef("call to lpcdd(, ...) (rational) returned NULL '$optimal.value' at point:")
      message.redef(q2d(point))
      save(hrep, point, file=paste("lpcddData_", .blackbox.data$options$jobSampleNbr, sep=""))
      message.redef(paste("'hrep' and 'point' are saved in file lpcddData_", .blackbox.data$options$jobSampleNbr, sep=""))
    } else {
      ## out$optimal.value>0 signifies point is out of  hull, cf vinny
        if (returnValue=="logical") rv <- !as.logical(as.numeric(out$optimal.value))
      if (returnValue=="max") rv <- as.numeric(out$optimal.value)
    }
  } else { ## using constraints
    if ( ! exactArithmetic) { ## floating point
      if (is.character(constraints$a[1])) {constraints$a <- q2d(constraints$a);constraints$b <- q2d(constraints$b)} ## stupidly inefficient case...
      if (is.character(point[1])) point <- q2d(point) ## stupidly inefficient case...
      mg <- constraints$a %*% point -constraints$b
    } else { ## rational
      if (is.numeric(constraints$a[1])) {constraints$a <- d2q(constraints$a);constraints$b <- d2q(constraints$b)}
      if (is.numeric(point[1])) point <- d2q(point)
      a.theta <- qmatmult(constraints$a, matrix(point))
      mg <- q2d(qmq(a.theta, constraints$b))
      mg <- q2d(mg)
    }
    ## in hull if all a%x-b =: mg are negative (vinny S2.1)
    if (returnValue=="logical") rv <- ! any(mg>0)
    if (returnValue=="max") rv <- max(mg) ## positive if out of hull
    if (returnValue=="vector") rv <- mg
  }
  return(rv)
} ## end def isPointInCHull
