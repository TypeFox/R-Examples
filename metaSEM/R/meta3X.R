meta3X <- function(y, v, cluster, x2, x3, av2, av3, data, intercept.constraints=NULL,
                   coef.constraints=NULL, RE2.constraints=NULL, RE2.lbound=1e-10,
                   RE3.constraints=NULL, RE3.lbound=1e-10, intervals.type=c("z", "LB"), R2=TRUE,
                   model.name="Meta analysis with ML", suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
  mf <- match.call()
  if (missing(data)) {
    data <- sys.frame(sys.parent())
  } else {
    if (!is.data.frame(data)) {
      data <- data.frame(data)
    }
  }
  my.y <- mf[[match("y", names(mf))]]
  my.v <- mf[[match("v", names(mf))]]
  my.cluster <- mf[[match("cluster", names(mf))]]  
  y <- eval(my.y, data, enclos = sys.frame(sys.parent()))  
  v <- eval(my.v, data, enclos = sys.frame(sys.parent()))
  cluster <- as.character(eval(my.cluster, data, enclos = sys.frame(sys.parent())))
  ## check if there are missing data in cluster
  if (any(is.na(cluster)))
      stop("Missing values are not allowed in \"cluster\".\n")

  if (missing(x2)&missing(x3)) 
    stop("Either \"x2\" or \"x3\" should be specified.\n")
  
  if (missing(x2)) {
    no.x2 <- 0
    x2 <- x2.labels <- l2.labels <- NULL
  } else {
    my.x2 <- mf[[match("x2", names(mf))]]
    x2 <- eval(my.x2, data, enclos = sys.frame(sys.parent()))
    if (is.vector(x2)) {
      no.x2 <- 1
      x2 <- matrix(x2, ncol=1)      
    } else {
      no.x2 <- ncol(x2)
    }
    x2.labels <- paste("x2_", 1:no.x2, sep="")
    l2.labels <- paste("l2_", 1:no.x2, sep="")
    colnames(x2) <- x2.labels
  }

  ## Check whether they are true level3 by testing the variance
  check.level3 <- function(x) { var(sapply(split(x, cluster), var, na.rm=TRUE), na.rm=TRUE)==0 }
  
  if (missing(x3)) {
    no.x3 <- 0
    x3 <- x3.labels <- l3.labels <- NULL
  } else {
    my.x3 <- mf[[match("x3", names(mf))]]
    x3 <- eval(my.x3, data, enclos = sys.frame(sys.parent()))
    if (is.vector(x3)) {
      if (!check.level3(x3))
        stop("\"x3\" is not a level 3 variable!\n")
      no.x3 <- 1
      x3 <- matrix(x3, ncol=1)
    } else {
      if (!all(apply(x3, 2, check.level3)))
        stop("Some variables in \"x3\" are not level 3 variables!\n")
      no.x3 <- ncol(x3)
    }
    x3.labels <- paste("x3_", 1:no.x3, sep="")
    l3.labels <- paste("l3_", 1:no.x3, sep="")    
    colnames(x3) <- x3.labels
  }
 
  if (missing(av2)&missing(av3)) {
    no.av2 <- no.av3 <- 0
    av2.labels <- av3.labels <- NULL
    ## Original data read
    ## Avoid errors in either x2 or x3 is NULL
    x <- cbind(x2, x3) 
    input.long <- data.frame(y, v, cluster, x)
  } else {
    if (missing(av2)) {
      no.av2 <- 0
      av2 <- av2.labels <- NULL
    } else {
      my.av2 <- mf[[match("av2", names(mf))]]
      av2 <- eval(my.av2, data, enclos = sys.frame(sys.parent()))
      if (is.vector(av2)) {
        no.av2 <- 1
        av2 <- matrix(av2, ncol=1)
      } else {
        no.av2 <- ncol(av2)
      }
      av2.labels <- paste("av2_", 1:no.av2, sep="")
      colnames(av2) <- av2.labels
    }
    if (missing(av3)) {
      no.av3 <- 0
      av3 <- av3.labels <- NULL
    } else {
      my.av3 <- mf[[match("av3", names(mf))]]
      av3 <- eval(my.av3, data, enclos = sys.frame(sys.parent()))        
      if (is.vector(av3)) {
        if (!check.level3(av3))
          stop("\"av3\" is not a level 3 variable!\n")
        no.av3 <- 1
        av3 <- matrix(av3, ncol=1)
      } else {
        if (!all(apply(av3, 2, check.level3)))
          stop("Some variables in \"av3\" are not level 3 variables!\n")
        no.av3 <- ncol(av3)
      }
      av3.labels <- paste("av3_", 1:no.av3, sep="")
      colnames(av3) <- av3.labels
    }
    x <- cbind(x2, av2, x3, av3)
    ## Original data read
    input.long <- data.frame(y, v, cluster, x)
  }

  ## Remove data with missing values in v
  miss.x <- is.na(v)
  ## Data used in the analysis
  my.long <- input.long[!miss.x, ]
  ## Reorder data according to clusters
  my.long <- my.long[order(my.long$cluster), ]

  ## starting values for x and av: excluded y, v and cluster
  mean.st <- round( apply(my.long[, -(1:3), drop=FALSE], 2, mean, na.rm=TRUE), 2 )
  
  ## starting values for within var  
  if ((no.x2+no.av2)!=0) {
    varW.st <- var( do.call(rbind, lapply( split( my.long[, c(x2.labels,av2.labels), drop=FALSE], my.long$cluster),
                                           function(x) scale(x, scale=FALSE))), na.rm=TRUE )
    varW.st <- round(varW.st, 2)
  }

  ## starting values for between var
  if (length(c(x2.labels,av2.labels,x3.labels,av3.labels))==1) {
    ## varB.st <- var( sapply(split( my.long[, c(x2.labels,av2.labels,x3.labels,av3.labels)],
    ##                               my.long[,"cluster"]), mean, na.rm=TRUE), na.rm=TRUE )
    ## no t()
    varB.st <- var( (sapply(split( my.long[, c(x2.labels,av2.labels,x3.labels,av3.labels), drop=FALSE],
                                    my.long[,"cluster"]), colMeans, na.rm=TRUE)), na.rm=TRUE )    
  } else {
    varB.st <- var( t(sapply(split( my.long[, c(x2.labels,av2.labels,x3.labels,av3.labels), drop=FALSE],
                                    my.long[,"cluster"]), colMeans, na.rm=TRUE)), na.rm=TRUE )
  }
  varB.st <- round(varB.st, 2)
     
  ## Convert long format to wide format as SEM uses wide format
  ## c() is required to convert matrix to vector when the data are balanced.
  ## c() is not required when the data are unbalanced.
  my.long$time <- c(unlist(sapply(split(my.long$y, my.long$cluster), function(x) 1:length(x))))
  long.labels <- names(my.long)
  my.wide <- reshape(my.long, timevar="time", idvar=c("cluster"), direction="wide",
                     v.names=c("y", "v", x2.labels, av2.labels), sep="_")
  
  ## ## Replace "." with "_" since OpenMx does not allow "." as variable names
  ## names(my.wide) <- sub("\\.", "_", names(my.wide))

  ## maximum no. of data in level-2 unit
  k <- max(sapply( split(my.long$cluster, my.long$cluster), length))
  
  ## NA in v is due to NA in y in wide format
  ## Replace with 1e10 though it does not affect the analysis as NA in y
  missing.v <- my.wide[, paste("v", 1:k, sep="_")]
  missing.v[is.na(missing.v)] <- 1e10
  my.wide[, paste("v", 1:k, sep="_")] <- missing.v
  
  ## Prepare matrices
  # intercept: mxmatrix containing the parameter
  if (is.null(intercept.constraints)) {
    Inter <- "0*Intercept"
  } else {
    if (length(intercept.constraints)>1)
      stop("A scalar is expected in \"intercept.constraints\".\n")
    Inter <- intercept.constraints
  }
  Inter <- as.mxMatrix(Inter)

  ## Variables are arranged in A, S as
  ## y, l2,av2,x2, l3,av3,x3
  
  ## Using x2 and x3 labels for l2 and l3
  intercept_x <- paste(mean.st, "*Mean_", names(mean.st), sep="")
  ## Fixed 0 for observed x2 and x3
  if (no.x2!=0)
    intercept_x <- c(intercept_x[1:(no.x2+no.av2)], rep(0,no.x2), intercept_x[-(1:(no.x2+no.av2))])
  if (no.x3!=0)
    intercept_x <- c(intercept_x, rep(0, no.x3))
  intercept_x <- matrix(intercept_x, nrow=1)
  intercept_x <- as.mxMatrix(intercept_x)

  ## names(inter.x) <- all.labels  
  Intercept1 <- mxAlgebra( cbind(Inter, intercept_x), name="Intercept1" )
  Unit1k <- mxMatrix("Unit", nrow=1, ncol=k, name="Unit1k")
  if ((no.x3+no.av3)==0) {
    ## No between structure; within is always present with y
    Imatrix <- mxAlgebra(Unit1k %x% Intercept1, name="Imatrix")
  } else {
    ## Repeat within k times + between
    ## Imatrix <- mxAlgebra( cbind( Unit1k %x% Intercept[1, 1:(1+no.x2*2+no.av2)],
    ##                              Intercept[1, -(1:(1+no.x2*2+no.av2))]), name="Imatrix")
    
    Itext <- paste("mxAlgebra(cbind(Unit1k %x% Intercept1[1, 1:(1+", no.x2, "*2+", no.av2, ")],
                                 Intercept1[1, -(1:(1+", no.x2, "*2+", no.av2, "))]), name=\"Imatrix\")", sep="")
    Imatrix <- eval(parse(text=Itext))
    
    ## Imatrix <- eval(substitute( mxAlgebra( cbind( Unit1k %x% Intercept[1, 1:(1+No.x2*2+No.av2)],
    ##                                        Intercept[1, -(1:(1+No.x2*2+No.av2))]), name="Imatrix"), list(No.x2 = no.x2, No.av2=no.av2)))
  }
  
  ## Regression coefficients
  if (is.null(coef.constraints)) {
    coef.constraints <- matrix(c(if (no.x2==0) NULL else paste("0*SlopeX2_", 1:no.x2, sep=""),
                                 if (no.x3==0) NULL else paste("0*SlopeX3_", 1:no.x3, sep="")), nrow=1)
  } else {
    ## convert into vector
    if (!is.matrix(coef.constraints))
      coef.constraints <- matrix(coef.constraints, nrow=1)
  }
    
  Beta <- as.mxMatrix(coef.constraints, name="Beta")    
   
  if (no.x2!=0) {
    Zero1x2 <- mxMatrix("Zero", nrow=1, ncol=(no.x2+no.av2), name="Zero1x2")
    ## Beta1 <- eval(substitute(mxAlgebra( cbind(Beta[1, 1:No.x2], Zero1x2, Beta[1, -(1:No.x2)]), name="Beta1"), list(No.x2 = no.x2)))
    B1text <- paste("mxAlgebra(cbind(Beta[1, 1:", no.x2, "], Zero1x2, Beta[1, -(1:",no.x2,")]), name=\"Beta1\")", sep="")
    Beta1 <- eval(parse(text=B1text))                                  
  } else {
    ## not actually use
    Zero1x2 <- mxMatrix("Zero", nrow=1, ncol=1, name="Zero1x2")
    Beta1 <- mxAlgebra( Beta, name="Beta1" )  
  }

  if (no.x3!=0) {
    Zero1x3 <- mxMatrix("Zero", nrow=1, ncol=(no.x3+no.av3), name="Zero1x3")
    Coeffs <- mxAlgebra( cbind(Beta1, Zero1x3), name="Coeffs")     
      ## coef.constraints <- c(coef.constraints, rep(0, no.x3+no.av3))
      ## coef.constraints <- t(as.matrix(coef.constraints))
  } else {
    ## not actually use
    Zero1x3 <- mxMatrix("Zero", nrow=1, ncol=1, name="Zero1x3")
    Coeffs <- mxAlgebra( Beta1, name="Coeffs") 
  }
  
  ## Including y,l2,av2,x2,l3,av3,x3
  ## Coeffs <- as.mxMatrix(coef.constraints, name="Coeffs")
  ## ## beta1 including y: y,l2,av2,x2,l3,av3,x3 by y,l2,av2,x2,l3,av3,x3
  ## beta1 <- matrix(0, nrow=(1+no.x2*2+no.av2+no.x3*2+no.av3), ncol=(1+no.x2*2+no.av2+no.x3*2+no.av3))

  beta.labels <- c("y",l2.labels,av2.labels,x2.labels,l3.labels,av3.labels,x3.labels)
  beta2 <- matrix(0, nrow=length(beta.labels), ncol=length(beta.labels), dimnames=list(beta.labels, beta.labels))
  ## Indicators for latent variables to observed variables
  if (no.x2!=0)
    beta2[x2.labels, l2.labels] <- Diag(no.x2)
  if (no.x3!=0)
    beta2[x3.labels, l3.labels] <- Diag(no.x3)
  ## Beta without y
  Beta2 <- as.mxMatrix(beta2[-1,-1], name="Beta2")
  Zeroall1 <- mxMatrix("Zero", nrow=length(beta.labels), ncol=1, name="Zeroall1")
  ## All coefficients for one unit: within+between
  Coefficients <- mxAlgebra( cbind(Zeroall1, rbind(Coeffs,
                                                   Beta2)), name="Coefficients")
  Unitkk <- mxMatrix("Unit", nrow=k, ncol=k, name="Unitkk")
  Idenkk <- mxMatrix("Iden", nrow=k, ncol=k, name="Idenkk")
  
  ## Amatrix <- mxAlgebra( rbind( cbind(Idenkk %x% Coefficients[1:(1+no.x2*2+no.av2), 1:(1+no.x2*2+no.av2)],
  ##                                    t(Unit1k) %x% Coefficients[1:(1+no.x2*2+no.av2), -(1:(1+no.x2*2+no.av2))]),
  ##                              cbind(Unit1k %x% Coefficients[-(1:(1+no.x2*2+no.av2)), 1:(1+no.x2*2+no.av2)],
  ##                                    Coefficients[-(1:(1+no.x2*2+no.av2)), -(1:(1+no.x2*2+no.av2))])), name="Amatrix")
  
  Atext <- paste("mxAlgebra( rbind( cbind(Idenkk %x% Coefficients[1:(1+", no.x2, "*2+", no.av2, "), 1:(1+", no.x2, "*2+", no.av2, ")],
                                     t(Unit1k) %x% Coefficients[1:(1+", no.x2, "*2+", no.av2, "), -(1:(1+", no.x2, "*2+", no.av2, "))]),
                               cbind(Unit1k %x% Coefficients[-(1:(1+", no.x2, "*2+", no.av2, ")), 1:(1+", no.x2, "*2+", no.av2, ")],
                                     Coefficients[-(1:(1+", no.x2, "*2+", no.av2,
                 ")), -(1:(1+", no.x2, "*2+", no.av2, "))])), name=\"Amatrix\")", sep="")
  Amatrix <- eval(parse(text=Atext))
  ## Amatrix <- eval(substitute(mxAlgebra( rbind( cbind(Idenkk %x% Coefficients[1:(1+No.x2*2+No.av2), 1:(1+No.x2*2+No.av2)],
  ##                                                    t(Unit1k) %x% Coefficients[1:(1+No.x2*2+No.av2), -(1:(1+No.x2*2+No.av2))]),
  ##                                              cbind(Unit1k %x% Coefficients[-(1:(1+No.x2*2+No.av2)), 1:(1+No.x2*2+No.av2)],
  ##                                                    Coefficients[-(1:(1+No.x2*2+No.av2)), -(1:(1+No.x2*2+No.av2))])), name="Amatrix"),
  ##                            list(No.x2=no.x2, No.av2=no.av2)))    
  
  ## v known for all large data
  ## all variables including latent variables in the analysis: used in the dataframe
  all.labels <- c(c(outer(c("y",l2.labels,av2.labels,x2.labels), 1:k, paste, sep="_")),
                  l3.labels, av3.labels, x3.labels) 
  no.all <- k*(1+no.x2*2+no.av2)+no.x3*2+no.av3
  v.labels <- c(unlist(lapply(paste("data.v_", 1:k, sep=""), function(x) c(x, rep(NA, no.x2*2+no.av2)))), rep(NA, no.x3*2+no.av3))
  ## v <- diag(v.labels)
  ## dimnames(v) <- list(all.labels, all.labels)
  V <- mxMatrix(type="Diag", nrow=no.all, ncol=no.all, free=FALSE, values=0,
                labels=v.labels, name="V")
  
  if ( length(RE2.lbound) != 1 ) {
    warning("\"RE2.lbound\" should be a scalar.")
    RE2.lbound <- 1e-10
  }
  if ( length(RE3.lbound) != 1 ) {
    warning("\"RE3.lbound\" should be a scalar.")
    RE3.lbound <- 1e-10
  }

  ## starting values
  tau2.st <-  round(var(unlist(lapply(split(my.long[, "y"], my.long$cluster), function(x) scale(x, scale=FALSE))), na.rm=TRUE), 2)
  tau3.st <-  round(var(sapply(split(my.long[, "y"], my.long[,"cluster"]), function(x) mean(x, na.rm=TRUE)), na.rm=TRUE), 2)
  
  if ( is.null(RE2.constraints) ) {
    tau2 <- paste(tau2.st, "*Tau2_2", sep="")
  } else {
    ## Tau2 <- as.mxMatrix(RE2.constraints, name="Tau2", lbound=RE2.lbound)
    if (length(RE2.constraints)!=1) 
      stop("Length of \"RE2.constraints\" is not 1.\n") 
  }

  if ( is.null(RE3.constraints) ) {
    tau3 <- paste(tau3.st, "*Tau2_3", sep="") 
  } else {
    ## Tau3 <- as.mxMatrix(RE3.constraints, name="Tau3", lbound=RE3.lbound)
    if (length(RE3.constraints)!=1) 
      stop("Length of \"RE3.constraints\" is not 1.\n")
  }

  ## Tau for both tau2 and tau3
  Tau <- as.mxMatrix( Diag(c(tau2, tau3)), lbound=matrix(c(RE2.lbound,NA,NA,RE3.lbound), nrow=2, ncol=2), name="Tau")

  if ((no.x2+no.av2)==0) {
    ## ZeroW1 and ZeroW2 are dummies, not used
    ZeroW1 <- mxMatrix("Zero", nrow=1, ncol=1, name="ZeroW1")
    ZeroW2 <- mxMatrix("Zero", nrow=1, ncol=1, name="ZeroW2")
    TauW <- mxAlgebra(Tau[1,1], name="TauW")    
  } else {
    ZeroW1 <- mxMatrix("Zero", nrow=(no.x2*2+no.av2), ncol=1, name="ZeroW1")
    ZeroW2 <- mxMatrix("Zero", nrow=(no.x2*2+no.av2), ncol=(no.x2*2+no.av2), name="ZeroW2")
    TauW <- mxAlgebra( rbind( cbind(Tau[1,1],t(ZeroW1)),
                              cbind(ZeroW1,ZeroW2) ), name="TauW" )
  }
 
  ZeroB1 <- mxMatrix("Zero", nrow=(no.x2*2+no.av2+no.x3*2+no.av3), ncol=1, name="ZeroB1")
  ZeroB2 <- mxMatrix("Zero", nrow=(no.x2*2+no.av2+no.x3*2+no.av3), ncol=(no.x2*2+no.av2+no.x3*2+no.av3), name="ZeroB2")
  TauB <- mxAlgebra( rbind( cbind(Tau[2,2],t(ZeroB1)),
                            cbind(ZeroB1,ZeroB2) ), name="TauB" )
  
  ## covW including y: y,l2,av2,x2 by y,l2,av2,x2
  w.labels <- c("y", l2.labels,av2.labels,x2.labels)
  covW <- matrix(0, nrow=length(w.labels), ncol=length(w.labels), dimnames=list(w.labels, w.labels))
  if ((no.x2+no.av2)!=0) {
    w.labels <- c(x2.labels, av2.labels)
    varW <- outer(w.labels, w.labels, function(x,y) { ifelse(x==y, paste("Wvar", x, sep=""),
                                                                 paste("Wcov", x, y, sep="")) } )
    varW <- paste(varW.st, varW, sep="*")
    varW <- matrix(varW, nrow=length(w.labels), ncol=length(w.labels))
    varW[lower.tri(varW)] <- t(varW)[lower.tri(varW)]
    covW[c(l2.labels, av2.labels), c(l2.labels, av2.labels)] <- varW
  }
  if (no.av2!=0)
    covW["y", av2.labels] <- covW[av2.labels, "y"] <- paste("0*Wcovy", av2.labels, sep="")

  ## covB including y: y,l2,av2,x2,l3,av3,x3 by y,l2,av2,x2,l3,av3,x3
  b.labels <- c("y",l2.labels,av2.labels,x2.labels,l3.labels,av3.labels,x3.labels)
  covB <- matrix(0, nrow=length(b.labels), ncol=length(b.labels), dimnames=list(b.labels, b.labels))
  b.labels <- c(x2.labels, av2.labels, x3.labels, av3.labels)
  varB <- outer(b.labels, b.labels, function(x,y) { ifelse(x==y, paste("Bvar", x, sep=""),
                                                           paste("Bcov", x, y, sep="")) } )
  varB <- paste(varB.st, varB, sep="*")
  varB <- matrix(varB, nrow=length(b.labels), ncol=length(b.labels))
  varB[lower.tri(varB)] <- t(varB)[lower.tri(varB)]
  covB[c(l2.labels, av2.labels, l3.labels, av3.labels), c(l2.labels, av2.labels, l3.labels, av3.labels)] <- varB
  if ((no.av2+no.av3)!=0)
    covB["y", c(av2.labels, av3.labels)] <- covB[c(av2.labels, av3.labels), "y"] <- paste("0*Bcovy", c(av2.labels, av3.labels), sep="")

  covW.lbound <- matrix(NA, nrow=dim(covW)[[1]], ncol=dim(covW)[[1]])
  Diag(covW.lbound) <- 1e-10
  covB.lbound <- matrix(NA, nrow=dim(covB)[[1]], ncol=dim(covB)[[1]])
  Diag(covB.lbound) <- 1e-10  
  covW <- as.mxMatrix(covW, lbound=covW.lbound)
  covB <- as.mxMatrix(covB, lbound=covB.lbound)
  
  expCW <- mxAlgebra(TauW+covW, name="expCW")
  expCB <- mxAlgebra(TauB+covB, name="expCB")  

  ## Total within
  ## TauT1 <- mxAlgebra( Unitkk %x% expCB[1:(1+no.x2*2+no.av2), 1:(1+no.x2*2+no.av2)] + Idenkk %x% expCW, name="TauT1")
  
  TauT1text <- paste("mxAlgebra( Unitkk %x% expCB[1:(1+", no.x2, "*2+", no.av2, "), 1:(1+", no.x2, "*2+", no.av2,
                     ")] + Idenkk %x% expCW, name=\"TauT1\")", sep="")
  TauT1 <- eval(parse(text=TauT1text))
  
  ## TauT1 <- eval(substitute(mxAlgebra( Unitkk %x% expCB[1:(1+No.x2*2+No.av2), 1:(1+No.x2*2+No.av2)] + Idenkk %x% expCW, name="TauT1"),
  ##                          list(No.x2=no.x2, No.av2=no.av2)))
                     
  ## Total within + between
  ## TauT2 <- mxAlgebra( rbind( cbind(TauT1, t(Unit1k) %x% expCB[1:(1+no.x2*2+no.av2), -(1:(1+no.x2*2+no.av2))]) ,
  ##                            cbind(Unit1k %x% expCB[-(1:(1+no.x2*2+no.av2)), 1:(1+no.x2*2+no.av2)], expCB[-(1:(1+no.x2*2+no.av2)), -(1:(1+no.x2*2+no.av2))]) ), name="TauT2")
  
  TauT2text <- paste("mxAlgebra( rbind( cbind(TauT1, t(Unit1k) %x% expCB[1:(1+", no.x2, "*2+", no.av2, "), -(1:(1+", no.x2, "*2+", no.av2,
                     "))]), cbind(Unit1k %x% expCB[-(1:(1+", no.x2, "*2+", no.av2, ")), 1:(1+", no.x2, "*2+", no.av2,
                     ")], expCB[-(1:(1+", no.x2, "*2+", no.av2, ")), -(1:(1+", no.x2, "*2+", no.av2,
                     "))]) ), name=\"TauT2\")", sep="")
  TauT2 <- eval(parse(text=TauT2text))
  
  ## TauT2 <- eval(substitute(mxAlgebra( rbind( cbind(TauT1, t(Unit1k) %x% expCB[1:(1+No.x2*2+No.av2), -(1:(1+No.x2*2+No.av2))]) ,
  ##                                            cbind(Unit1k %x% expCB[-(1:(1+No.x2*2+No.av2)), 1:(1+No.x2*2+No.av2)],
  ##                                                  expCB[-(1:(1+No.x2*2+No.av2)), -(1:(1+No.x2*2+No.av2))]) ), name="TauT2"),
  ##                          list(No.x2=no.x2, No.av2=no.av2)))
    
  Smatrix <- mxAlgebra( TauT2 + V, "Smatrix")

  ## TauT2 <- mxAlgebra( cbind(TauT1, t(Unit1) %x% expCB[1:(1+no.x2*2+no.av2), -(1:(1+no.x2*2+no.av2))]) , name="TauT2")
  ## TauT2 <- mxAlgebra( t(Unit1) %x% expCB[1:(1+no.x2*2+no.av2), -(1:(1+no.x2*2+no.av2))] , name="TauT2")
  
  ## my.model <- mxModel("test", Tau2, Tau3, ZeroW1, ZeroW2, TauW, ZeroB1, ZeroB2, TauB, covW, covB,
  ##                     expCW, expCB, Iden2, TauT1, Unit1, Unit2, TauT2)
  ## my <- mxRun(my.model)
  ## my@algebras
  ## my@matrices
  
  ## y, x2,x3, l2,l3, av2,av3
  if (no.x2==0) {
    Fmatrix <- create.Fmatrix( (!all.labels %in% l3.labels), name="Fmatrix" )
  } else {
    temp.labels <- c(outer(l2.labels, 1:k, function(x,y) paste(x,y,sep="_")), l3.labels)
    Fmatrix <- create.Fmatrix( !all.labels %in% temp.labels, name="Fmatrix" )
  }

  mx.model <- mxModel(model=model.name, mxData(observed=my.wide[,-1], type="raw"),
                      Inter, intercept_x, Intercept1, Unit1k,
                      Beta, Coeffs, Beta1, Beta2, Zero1x2, Zero1x3, Zeroall1, Coefficients, Unitkk,
                      Tau, ZeroW1, ZeroW2, TauW,
                      ZeroB1, ZeroB2, TauB, covW, covB,
                      expCW, expCB, Idenkk, TauT1, TauT2, Smatrix, V,
                      Amatrix, Smatrix, Imatrix, Fmatrix,                      
                      mxExpectationRAM(A="Amatrix", S="Smatrix", F="Fmatrix", M="Imatrix", dimnames=all.labels),
                      mxFitFunctionML(),
                      mxCI(c("Inter","Beta","Tau")))

  ## V <- mxMatrix(type="Diag", nrow=no.all, ncol=no.all, free=FALSE, values=0,
  ##               labels=NA, name="V")
  ## Iden <- mxMatrix("Iden", nrow=no.all, ncol=no.all, name="Iden")
  ## model.name <- "test"
  ## mx.model <- mxModel(model=model.name, 
  ##                     intercept_y, intercept_x, Intercept, Unit1k,
  ##                     Coeffs, Beta1, Zeroall1, Coefficients, Unitkk,
  ##                     Tau2, Tau3, ZeroW1, ZeroW2, TauW,
  ##                     ZeroB1, ZeroB2, TauB, covW, covB,
  ##                     expCW, expCB, Idenkk, TauT1, TauT2, Smatrix, V, Iden,
  ##                     Amatrix, Smatrix, Mmatrix, Fmatrix)
  
  ## Calculate R2
  if (R2) {
    mx0text <- "tryCatch( meta3X(y=y, v=v, cluster=cluster,"
    if (no.x2!=0)
        mx0text <- paste(mx0text, " x2=my.long[, x2.labels],", sep="")
    if (no.x3!=0)
        mx0text <- paste(mx0text, " x3=my.long[, x3.labels],", sep="")  
    if (no.av2!=0)
        mx0text <- paste(mx0text, " av2=my.long[, av2.labels],", sep="")
    if (no.av3!=0)
        mx0text <- paste(mx0text, " av3=my.long[, av3.labels],", sep="")
    mx0text <- paste(mx0text, "data=my.long, intercept.constraints=intercept.constraints,
                               coef.constraints=matrix(0, nrow=1, ncol=(no.x2+no.x3)),
                               RE2.constraints=RE2.constraints,
                               RE3.constraints=RE3.constraints,
                               intervals.type=\"z\", R2=FALSE,
                               suppressWarnings=TRUE, silent=TRUE), error = function(e) e )", sep="")
    mx0.fit <- eval(parse(text=mx0text))
    ## mx0.fit <- tryCatch( meta3ML(y=y, v=v, cluster=cluster, x2=x2, x3=x3, av2=av2,
    ##                              av3=av3, model.name="No predictor",
    ##                              intercept.constraints=intercept.constraints,
    ##                              coef.constraints=matrix(0, nrow=1, ncol=(no.x2+no.x3)),
    ##                              RE2.constraints=RE2.constraints,
    ##                              RE3.constraints=RE3.constraints,
    ##                              intervals.type="z", R2=FALSE,
    ##                              suppressWarnings=TRUE, silent=TRUE), error = function(e) e )          
  } else {
    mx0.fit <- NA
  }

  ## Return mx model without running the analysis
  if (run==FALSE) return(mx.model)
  
  intervals.type <- match.arg(intervals.type)
  # Default is z
  switch(intervals.type,
    z = mx.fit <- tryCatch( mxRun(mx.model, intervals=FALSE, suppressWarnings=suppressWarnings,
                                  silent=silent, ...), error = function(e) e ),
    LB = mx.fit <- tryCatch( mxRun(mx.model, intervals=TRUE, suppressWarnings=suppressWarnings,
                                   silent=silent, ...), error = function(e) e ) )
 
  if (inherits(mx.fit, "error")) {
    cat("Error in running the mxModel:\n")
    warning(print(mx.fit))
  }

  ## ## my.long is complete data
  ## ## FIXME: remove miss.x (is it used by meta()?)
  out <- list(call=mf, R2=R2, data.wide=my.wide, data=my.long,
              no.y=1, no.x2=no.x2, no.x3=no.x3, no.av2=no.av2, no.av3=no.av3,
              miss.x=rep(FALSE, nrow(my.long)), mx.model=mx.model,
              mx.fit=mx.fit, mx0.fit=mx0.fit, intervals.type=intervals.type)
  class(out) <- "meta3X"
  out
}

