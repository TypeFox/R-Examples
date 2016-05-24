## **** 1D profile CI bounds ****
bounds1D <- function(level=log(0.05), string1, string2, lowval, hival,
                     MLval, # must be [NA in some cases] as tested below
                     CIvar, precision="rational") {
  fitobject <- blackbox.getOption("fitobject")
  fittedNames <- blackbox.getOption("fittedNames")
  fittedparamnbr <- blackbox.getOption("fittedparamnbr")
  notinKgspace <- CIvar %w/o% fittedNames ## e.g. latt2Ns2 %w/o% canonicals is latt2Ns2
  locchull <- providefullhull(notinKgspace)[[1]] ## 'full' dimensional, vertices constraints and affectedConstraints
  profiledNames <- setdiff(colnames(locchull$vertices), CIvar) ## twoNmu, g in above example
  if (length(notinKgspace)==0){ ## no extra composite hull
    if(is.na(MLval)) message.redef("(!) From bounds1D(): MLval is NA were a value is expected")
    ## MLval <- rosglobal$par[CIvar]
    givenmax <- blackbox.getOption("rosglobal")$value
  } else {
    ## ML in composite space hull
    if( ! is.na(MLval)) message.redef("(!) From bounds1D(): MLval has a value were NA is expected")
    initpt <- fitobject$x[which.max(fitobject$fitted.values), ]
    initpt <- fromFONKtoanyspace(initpt, colnames(locchull$vertices))
    ###################
    if (precision=="rational") {
      ui <- qneg(locchull$a);ci <- qneg(locchull$b)
    } else {ui <- -locchull$a;ci <- -locchull$b}    ##  within-hullness occurs for ui %*% x -ci => 0
    candidates <- generateInitpts(bestProfiledOut=initpt, vertices=locchull$vertices, ui=ui, ci=ci, hrepr=locchull$Hrep, fixedlist=NA, precision=precision)
    #        candidates <- generateInitpts(bestProfiledOut=initpt, locsubchull=locchull$vertices, ui=ui, ci=ci, fixedlist=NA, precision=precision)
    ###################
    bestresu <- list(value=- Inf)
    for (candidate in candidates) {
      resu <- optimWrapper( ##purefn,
        initval=candidate$par, gr=NULL, ## not initpt, 14/10/2011
        chullformats=locchull,
        control=list( ## parscale is provided within optimWrapper
          fnscale=-1/blackbox.getOption("scalefactor"), trace=FALSE, maxit=10000))
      if (resu$value>bestresu$value) bestresu <- resu
    } ## end loop over candidate starting points
    maxpt <- bestresu
    maxpt$par <- fromFONKtoanyspace(maxpt$par, colnames(locchull$vertices))
    MLval <- maxpt$par[CIvar]
    givenmax <- maxpt$value
  }
  lowval <- lowval+0.002*(MLval-lowval)
  hival <- hival-0.002*(hival-MLval)
  ## def objectivefn
  shift <- - givenmax - level
  if(fittedparamnbr==1) {
    objectivefn1D <- function(x) { ## defined to return 0 at the CI threshold
      return(purefn(x, testhull=F)+shift)
    }
    objectivefn <- objectivefn1D
  } else {
    argfixedlist <- eval(parse(text=paste("list(", CIvar, "=", NA, ")"))) ## je pourrais creer une list et nommer ensuite... mais cette syntaxe pourra resservir
    objectivefnmultiD <- function(x, return.optim=FALSE) { ## defined to return 0 at the CI threshold
      argfixedlist[1] <- x
      grmf <- profileBySubHull(fixedlist=argfixedlist, extrap=FALSE) ## metric arg
      ##############################    if ( ! (is.numeric(grmf$value) && is.finite(grmf$value))) browser()
      if(return.optim) {
        return(grmf)
      } else return(grmf$value+shift) ## then returns shifted value for uniroot
    }
    objectivefn <- objectivefnmultiD
  }
  CIlo <- NA
  CIpoints <- matrix(nrow=0, ncol=fittedparamnbr)
  fupper <- -level
  if (abs(lowval-MLval)<abs(MLval*1e-08)) { ## not an error ## lowval==MLval if relative diff is <1e-16 in absolute value
    # CIlo <- NA ## already so
  } else {
    flower <- objectivefn(lowval)
    if (is.na(flower)) { ## abnormal case
      stop.redef("From bounds1D: 'flower' is NA")
    } else {
      if (flower<0) {CIlo <- try((uniroot(objectivefn, interval=c(lowval, MLval),
                                        f.lower=flower, f.upper=fupper))$root, T)
      }
    }
    if(inherits(CIlo,"try-error")) {
      CIlo <- NA
      errmsg <- paste(string1, " for ", CIvar, " could not be computed (maybe out of sampled range)", sep="")
      message.redef(errmsg)
    }
    if ( ! is.na(CIlo)) {
      CIpoint <- numeric()
      CIpoint[CIvar] <- CIlo
      if (length(profiledNames)>0) {
        thisfit <- objectivefn(CIlo, return.optim=TRUE)
        CIpoint[profiledNames] <- thisfit$par
        CIpoint <- tofullKrigingspace(CIpoint) ## will be used by predict.HLfit
      }
      CIpoints <- rbind(CIpoints, CIpoint)
    }
  }
  CIup <- NA ## default values...
  flower <- fupper
  if (abs(hival-MLval)<abs(MLval*1e-08)) { ## not an error
    ## CIup <- NA
  } else {
    fupper <- objectivefn(hival)
    if (is.na(fupper)) { ## abnormal case
      stop.redef("From bounds1D: 'fupper' is NA")
    } else {
      if (fupper<0) {CIup <- try((uniroot(objectivefn, c(MLval, hival),
                                        f.lower=flower, f.upper=fupper))$root, TRUE)
      }
    }
    if(inherits(CIup,"try-error")) {
      CIup <- NA
      errmsg <- paste(string2, " for ", CIvar, " could not be computed (maybe out of sampled range)", sep="")
      message.redef(errmsg)
    }
    if ( ! is.na(CIup)) {
      CIpoint <- numeric()
      CIpoint[CIvar] <- CIup
      if (length(profiledNames)>0) {
        thisfit <- objectivefn(CIup, return.optim=TRUE)
        CIpoint[profiledNames] <- thisfit$par
        CIpoint <- tofullKrigingspace(CIpoint) ## will be used by predict.HLfit
      }
      CIpoints <- rbind(CIpoints, CIpoint)
    }
  }
  return(list(CIlo=CIlo, CIup=CIup, CIpoints=CIpoints))
} ## end def bounds1D
