"get.residuals" <- function(object, units="ignore")
{ options <- c("addtoresiduals", "ignore")
  unit.opt <- options[check.arg.values(units, options)]
  if (length(object$aom) == 0)
  { res <- as.vector(object$residuals)
    attr(res, which="restype") <- "Residuals"
    stdres <- FALSE
  } else
  { res <- as.vector(object$aom$R[,2])
    attr(res, which="restype") <- "Standardized conditional residuals"
    stdres <- TRUE
  }
  if (unit.opt == "addtoresiduals")
  { #Check whether units is in random model
    termno <- findterm("units", rownames(summary(object)$varcomp))
    if (termno > 0)
    { if (!stdres)
        ucoeff <- 1
      else
      { ucoeff <- summary(object)$varcomp[termno, "component"]
        ucoeff <- sqrt(ucoeff/(ucoeff^2 -1))
        #uBLUP <- object$aom$G[,2]
      }  
      uBLUP <- object$coefficients$random
      uBLUP <- ucoeff*uBLUP[grep("units_", substr(names(uBLUP), 1, 6), fixed=TRUE)]
      names(uBLUP) <- NULL
      res <- res + uBLUP
    }
  }
  res[is.nan(res)] <- NA
  return(res)
}


"variofaces.asreml" <- function(object, V, nsim=100, seed = NULL, tolerance = 1E-10, 
                                units = "ignore", update = TRUE, trace = FALSE, 
                                graphics.device=NULL, ...)
#function to do the face variogram plots, including envelopes, described by 
#Stefanova et al (2010)
#object is an asreml object from a call to asreml in which the data argument 
#   must have been set.
#V is the fitted variance matrix i.e. having the appropriate pattern and values 
#   given the model fitted and the estimates of the parameters obtained
#nsim is the number of data sets to be simulated in obtaining the envelopes
#   Note that only the results from data sets that converge in the fitting are 
#        included in computing the envelope
#tolerance is the value such that eigenvalues less than it are consdered to be zero
#... parameters to supply to plot functions called within variofaces.asreml
{ options <- c("addtoresiduals", "ignore")
  unit.opt <- options[check.arg.values(units, options)] 
 if (class(object) != "asreml")
   stop("object should be an asreml object")
  #Check V
  if (!isSymmetric(V))
    stop("Variance matrix must be symmetric")
  n <- length(object$residuals)
  if (!all(dim(V) == c(n, n))) 
      stop("V is not a square matrix whose order equals the length of the response variable")
  #use eigenvalue decomposition to establish transformation matrix
  eigdecomp <- eigen(V, symmetric = TRUE)
  eigenval <- eigdecomp$values
  if (!all(eigenval >= -tolerance * abs(max(eigenval))))
     stop("Variance matrix is not nonnegative definite")
  R <- eigdecomp$vectors %*% diag(sqrt(pmax(eigenval, 0)), n)

 
  #check and get info in supplied call
  call <- object$call
  if (!("data" %in% names(call)))
    stop("variplot.asreml assumes that data has been set in call to asreml")
  env.dat <- eval(call$data)
 
  #set up call
  object.sim <- object
  call <- object.sim$call
  elno <- grep("data", names(call))
  languageEl(call, which=elno) <- quote(env.dat)
  languageEl(call$fixed, which=2) <- quote(y.sim)
  call$trace <- trace

 #Deal with the R.param and G.param arguments
 if (update)
 { #If update, set R.param and G.param
   languageEl(call, which = "R.param") <- object$R.param
   languageEl(call, which = "G.param") <- object$G.param
 } else
 { #If R.param and G.param already set, make them NULL
   if (!is.null(languageEl(call, which = "R.param")))
     languageEl(call, which = "R.param") <- NULL
   if (!is.null(languageEl(call, which = "G.param")))
     languageEl(call, which = "G.param") <- NULL
 }
 
 #deal with args coming via ...
 #- this will overwite previously set values, except data and models are protectd
 tempcall <- list(...)
 if (length(tempcall)) 
 { for (z in names(tempcall))
   if (z == "data")
     env.dat <- data
   else
     if (z %in% c("fixed", "random", "rcov", "sparse"))
       stop("attempt to change model to be fitted")
   else
     languageEl(call, which = z) <- tempcall[[z]]
 }
 
 #investigate rcov term to see if it is a two-factor term or a 
 #three-factor term with one factor indexing sections
 rterm <- languageEl(call, which = "rcov")
 #Check for multiple terms
 rterm.form <- as.formula(rterm)
 rterm.obj <- as.terms.object(rterm, object, 
                              specials=c("at","ar1", "ar2", "ar3", "sar","sar2",
                                         "ma1", "ma2", "arma", "exp", "gau", 
                                         "cor", "corb", "corg", "diag", "us", "ante"))
 if (length(labels(rterm.obj)) != 1)
   stop("In analysing ",object$fixed.formula[[2]],
        ", the rcov model must involve a single term")
 grid.facs <- rownames(attr(rterm.obj, which = "factors"))
 grid.facs <- unlist(lapply(grid.facs, rmFunction))
 sections <- 1
 fac.sec <- NULL
 kspecial <- attr(rterm.obj, which = "specials")$at
 if (!is.null(kspecial))
 { if (length(kspecial) != 1)
     stop("Can only have a single factor defining sections of the data (using at)")
   fac.sec <- grid.facs[kspecial[1]]  
   sections <- levels(env.dat[[fac.sec]])
   grid.facs <- grid.facs[-kspecial[1]]
 }
 if (length(grid.facs) != 2)
   stop("Must have two dimensions in the rcov model")

 #Set up for simulation
 #There will be two lists, one for each face
 #The list for each face  has number of sections components
 #Each component will be an nsim x nfaci matrix, where i is the face 
 fac1 <- grid.facs[1]
 fac2 <- grid.facs[2]
 if (any(table(env.dat[[fac1]])==0) | any(table(env.dat[[fac2]])==0))
   stop("Some levels in the the supplied factors are not observed") 
 n1 <- length(levels(env.dat[[fac1]]))
 n2 <- length(levels(env.dat[[fac2]]))
 ns1 <- nsim*n1
 ns2 <- nsim*n2
 env.var <- vector(mode = "list", length = 2)
 env.var[[1]] <- vector(mode = "list", length = length(sections))
 names(env.var[[1]]) <- sections
 env.var[[1]] <- lapply(env.var[[1]][sections], 
                        function(env.var, nrow, ncol)
                          env.var <- matrix(rep(0, nrow*ncol), nrow=nrow, ncol=ncol),
                        nrow=nsim, ncol=n1)
 env.var[[2]] <- vector(mode = "list", length = length(sections))
 names(env.var[[2]]) <- sections
 env.var[[2]] <- lapply(env.var[[2]][sections], 
                        function(env.var, nrow, ncol)
                          env.var <- matrix(rep(0, nrow*ncol), nrow=nrow, ncol=ncol),
                        nrow=nsim, ncol=n2)
 
 
 #process seed argument
 if (!is.null(seed))
   set.seed(seed)
 
 #generate nsim data sets and save variogram face data
  conv <- FALSE
  res.dat <- data.frame(env.dat[[fac1]],env.dat[[fac2]])
  names(res.dat)[1:2] <- grid.facs
  res.dat[grid.facs] <- lapply(res.dat[grid.facs], as.numfac)
   if (!is.null(fac.sec))
   { res.dat <- data.frame(env.dat[[fac.sec]], res.dat)
     names(res.dat)[1] <- fac.sec
   } 
   mu <- fitted.values(object)
  for (i in 1:nsim)
  { while (!conv)
    { env.dat$y.sim <- as.vector(mu + R %*% rnorm(n))
      sim.asreml <- eval(call)
      conv <- sim.asreml$converge
    }
    conv=FALSE
    #Get residuals
    res.dat$res <- get.residuals(sim.asreml, units = units)
    #Get variogram faces for this simulated data set
    for (k in sections)
    { if (is.null(fac.sec))
        sect.dat <- res.dat
      else
        sect.dat <- res.dat[res.dat[[fac.sec]]==k, -1]
      sect.dat <- as.matrix(sect.dat)
      sim.var <- asreml::asreml.variogram(sect.dat)
      names(sim.var)[1:2] <- grid.facs
      env.var[[1]][[k]][i,] <- sim.var[sim.var[[fac2]]==0,]$gamma
      env.var[[2]][[k]][i,] <- sim.var[sim.var[[fac1]]==0,]$gamma
    }
  }
  restype <- attr(res.dat$res, which = "restype")
 
  #Get observed residuals
  res.dat$res <- get.residuals(object, units = units)
  if (restype != attr(res.dat$res, which = "restype"))
    warning("Observed and simulated residuals are not of the same type")

  #Do face plots for fac1 and fac2 in each section
  face.1 <- data.frame(matrix(nrow = 0, ncol=5))
  colnames(face.1) <- c("factor","observed","X2.5.","X50.","X97.5.")
  face.2 <- data.frame(matrix(nrow = 0, ncol=5))
  colnames(face.2) <- c("factor","observed","X2.5.","X50.","X97.5.")
  
  #Get the observed variogram and form the elements of the variofaces plot
  for (k in sections)
  { if (is.null(fac.sec))
      sect.dat <- res.dat
    else
      sect.dat <- res.dat[res.dat[[fac.sec]]==k, -1]
    sect.dat <- as.matrix(sect.dat)
    object.var <- asreml::asreml.variogram(sect.dat)
    names(object.var)[1:2] <- grid.facs
    
    #Form data.frame for current section
    face.1 <- rbind(face.1,
                    data.frame(factor = object.var[object.var[[fac2]]==0,][[fac1]],
                               observed = object.var[object.var[[fac2]]==0,]$gamma,   
                               t(apply(env.var[[1]][[k]], 2, quantile, 
                                       probs=c(0.025, 0.50, 0.975),
                                       na.rm = TRUE)),
                               stringsAsFactors = FALSE))
    face.2 <- rbind(face.2,
                    data.frame(factor = object.var[object.var[[fac1]]==0,][[fac2]],   
                               observed = object.var[object.var[[fac1]]==0,]$gamma,   
                               t(apply(env.var[[2]][[k]], 2, quantile, 
                                       probs=c(0.025, 0.50, 0.975),
                                       na.rm = TRUE)), 
                               stringsAsFactors = FALSE))
  }
  names(face.1)[1] <- fac1
  names(face.2)[1] <- fac2
 
  if (!is.null(fac.sec))
  { face.1[[fac.sec]] <- factor(rep(1:length(sections), each=(n1)), labels=sections)
    face.2[[fac.sec]] <- factor(rep(1:length(sections), each=(n2)), labels=sections)
  }

  #Do plots
  if (!is.null(graphics.device) )
    do.call(graphics.device, list(record = FALSE))
 
  p <- ggplot(data=face.1) +
       theme_bw() +
       geom_line(aes_string(x= fac1, y = "observed"), size=1) +
       geom_point(aes_string(x= fac1, y = "observed"), size=3) +
       geom_line(aes_string(x=fac1, y = "X2.5."), colour="red") +
       geom_line(aes_string(x=fac1, y = "X97.5."), colour="red") +
       geom_line(aes_string(x=fac1, y = "X50."), colour="blue") +
       labs(title = paste("Variogram face of",restype,"for",fac1,sep=" "), 
            x = paste(fac1,"differences", sep=" "),y=NULL)
  if (!is.null(fac.sec))
    p <- p + facet_wrap(as.formula(paste("~",fac.sec)), ncol=2, scales="free")
  print(p)

  p <- ggplot(data=face.2) +
       theme_bw() +
       geom_line(aes_string(x= fac2, y = "observed"), size=1) +
       geom_point(aes_string(x= fac2, y = "observed"), size=3) +
       geom_line(aes_string(x=fac2, y = "X2.5."), colour="red") +
       geom_line(aes_string(x=fac2, y = "X97.5."), colour="red") +
       geom_line(aes_string(x=fac2, y = "X50."), colour="blue") +
       labs(title = paste("Variogram face of",restype,"for",fac2,sep=" "), 
            x = paste(fac2,"differences", sep=" "), y=NULL)
  if (!is.null(fac.sec))
    p <- p + facet_wrap(as.formula(paste("~",fac.sec)), ncol=2, scales="free")
  print(p)
 
   #return data frames containing the variogram values on which the plots are based
  invisible(list(face1 = face.1, face2 = face.2))
}

"simulate.asreml" <- function(object, nsim=100, seed = NULL, means=NULL, V, 
                              tolerance = 1E-10, update = TRUE, trace = FALSE, 
                              which = "data", units = "ignore", ...)
  #function to obtain simulated data corresponding to a fitted asreml model
  #object is an asreml object from a call to asreml in which the data argument 
  #   must have been set.
  #means is a vector of predictions for fixed terms - should include spline terms
  #   if null, fitted values from fitted model are used.
  #V is the fitted variance matrix i.e. having the appropriate pattern and values 
  #   given the model fitted and the estimates of the parameters obtained
  #nsim is the number of data sets to be simulated in obtaining the envelopes
  #   Note that only the results from data sets that converge in the fitting are 
  #        included in computing the envelope
  #tolerance is the value such that eigenvalues less than it are consdered to be zero
#... parameters to supply to plot functions called within variofaces.asreml
{ options <- c("data", "fitted", "residuals", "all")
  opt <- options[unlist(lapply(which, check.arg.values, options=options))]
  if ("all" %in% opt)
    opt <-  c("data", "fitted", "residuals")
  options <- c("addtoresiduals", "ignore")
  unit.opt <- options[check.arg.values(units, options)]
  if (class(object) != "asreml")
    stop("object should be an asreml object")
  n <- length(object$residuals)
  if (!is.null(means) & length(means) != n)
    stop("The lengths of means  and the response variable are not the same")
  #Check V
  if (!isSymmetric(V))
    stop("Variance matrix must be symmetric")
  if (!all(dim(V) == c(n, n))) 
    stop("V is not a square matrix whose order equals the length of the response variable")
  #use eigenvalue decomposition to establish transformation matrix
  eigdecomp <- eigen(V, symmetric = TRUE)
  eigenval <- eigdecomp$values
  if (!all(eigenval >= -tolerance * abs(max(eigenval))))
    stop("Variance matrix is not nonnegative definite")
  R <- eigdecomp$vectors %*% diag(sqrt(pmax(eigenval, 0)), n)
  
  
  #check and get info in supplied call
  call <- object$call
  if (!("data" %in% names(call)))
    stop("simulate.asreml assumes that data has been set in call to asreml")
  env.dat <- eval(call$data)
  
  #set up call
  object.sim <- object
  call <- object.sim$call
  elno <- grep("data", names(call))
  languageEl(call, which=elno) <- quote(env.dat)
  languageEl(call$fixed, which=2) <- quote(y.sim)
  call$trace <- trace
  
  #Deal with the R.param and G.param arguments
  if (update)
  { #If update, set R.param and G.param
    languageEl(call, which = "R.param") <- object$R.param
    languageEl(call, which = "G.param") <- object$G.param
  } else
  { #If R.param and G.param already set, make them NULL
    if (!is.null(languageEl(call, which = "R.param")))
      languageEl(call, which = "R.param") <- NULL
    if (!is.null(languageEl(call, which = "G.param")))
      languageEl(call, which = "G.param") <- NULL
  }
  
  #deal with args coming via ...
  #- this will overwite previously set values, except data and models are protectd
  tempcall <- list(...)
  if (length(tempcall)) 
  { for (z in names(tempcall))
    if (z == "data")
      env.dat <- data
    else
      if (z %in% c("fixed", "random", "rcov", "sparse"))
        stop("attempt to change model to be fitted")
    else
      languageEl(call, which = z) <- tempcall[[z]]
  }
  
  #Set up for simulation and get observed residuals and fitted values
  #There will be a data.fame for each quantity saved, as well as the supplied 
  #data.frame to which the requested quantites have been added
  if ("data"%in% opt)
  { data <- data.frame(matrix(nrow = n, ncol=nsim))
    names(data) <- paste("data",1:nsim,sep=".")
  }
  if ("fitted"%in% opt)
  { fitted <- data.frame(matrix(nrow = n, ncol=nsim))
    names(fitted) <- paste("fitted",1:nsim,sep=".")
    env.dat$fitted <- fitted(object)
  } 
  if ("residuals"%in% opt)
  { residuals <- data.frame(matrix(nrow = n, ncol=nsim))
    names(residuals) <- paste("residuals",1:nsim,sep=".")
    #Get observed residuals
    env.dat$residuals <- get.residuals(object, units = units)
    restype <- attr(env.dat$residuals, which = "restype")
  }
  
  #Set up expectation
  if (is.null(means))
    mu <- fitted.values(object)
  else
    mu <- means
  
  #Check type of residuals in simulation
  env.dat$y.sim <- as.vector(mu + R %*% rnorm(n))
  sim.asreml <- eval(call)
  if ("residuals" %in% opt)
    if (((length(sim.asreml$aom) == 0) & (restype == "Standardized conditional residuals")) | 
          ((length(sim.asreml$aom) != 0) & (restype == "Residuals")))
      warning("Observed and simulated residuals are not of the same type\n",
              "- check setting of aom")
  
  #process seed argument
  if (!is.null(seed))
    set.seed(seed)
  
  #generate nsim data sets and save variogram face data
  conv <- FALSE
  for (i in 1:nsim)
  { while (!conv)
  { env.dat$y.sim <- as.vector(mu + R %*% rnorm(n))
    sim.asreml <- eval(call)
    conv <- sim.asreml$converge
  }
  conv=FALSE
  if ("data"%in% opt)
    data[i] <- env.dat$y.sim
  #Get residuals
  if ("residuals"%in% opt)
    residuals[i] <- get.residuals(sim.asreml, units=units)
  if ("fitted"%in% opt)
    fitted[i] <- fitted.values(sim.asreml)
  }
  
  out <- vector("list", length = 0)
  if ("residuals" %in% opt | "fitted" %in% opt)
  { env.dat <- env.dat[,-match("y.sim",names(env.dat))]
    out[["observed"]] <- env.dat
  }
  if ("data" %in% opt)
    out[["data"]] <- data
  if ("fitted" %in% opt)
    out[["fitted"]] <- fitted
  if ("residuals" %in% opt)
    out[["residuals"]] <- residuals
  #return data frames containing the variogram values on which the plots are based
  invisible(out)
}


"plotvariofaces.asreml" <- function(data, residuals, restype="Residuals", ...)
  #function to do the face variogram plots, including envelopes, described by 
  #Stefanova et al (2010)
  #data is a data.frame with 3 or 4 columns. If there are 4, the first is taken to index
  #     sections of the data, the second and the third to index the x and y dimensions 
  #     underlying the observations and the fourth the observed residuals. If there are 3, 
  #     then the same as for 4 columns, except that the factor index the sections has 
  #     been omitted.
  #residuals is a data frame with the same initial columns as data, the observed residuals 
  #     having been ommitted. Following the initial columns are nsim sts of residuals 
  #     derived from simulated data. 
  #restype is a character string describing the type of residuals supplied.
{ 
  if (!is.data.frame(data) | !is.data.frame(residuals))
    stop("Both data and residuals should be data frames")
  n <- nrow(data)
  nvars <- ncol(data)
  if (nvars < 3 & nvars > 4)
    stop("The number of columns in data is not 3 or 4")
  nsim <- ncol(residuals) - nvars
  if (nsim < 1)
    stop("Could not find columns with simulated residuals")
  #Check compatibility
  if (n != nrow(residuals))
    stop("The number of rows in data and residuals are not equal")
  
  #Check that initial columns in data and residuals match
  if (!setequal(names(data)[1:(nvars-1)], names(residuals)[1:(nvars-1)]))
    stop("The names of the initial ",nvars-1," columns do not match")
  if (!all(unlist(lapply(data[1:(nvars-1)], is.factor))))
    stop("Some of the initial ",nvars-1," columns of data are not factors")
  if (!all(unlist(lapply(residuals[1:(nvars-1)], is.factor))))
    stop("Some of the initial ",nvars-1," columns of residuals are not factors")
  for (i in 1:(nvars-1))
    if (!all(data[i] == residuals[i]))
      stop("The values in column ",i," of data and residuals do not match")
  
  #Identify columns in data
  if (nvars == 3)
  { fac.sec <- NULL
    sections <- 1
    grid.facs <- colnames(data[1:2])
    residuals <- colnames(data[3])
  } 
  else
  { fac.sec <- colnames(data[1])
    sections <- levels(data[[fac.sec]])
    grid.facs <- colnames(data[2:3])
    residuals <- colnames(data[4])
  }
  
  #Form variogram faces for each set of residuals
  #There will be two lists, one for each face
  #The list for each face  has number of sections components
  #Each component will be an nsim x nfaci matrix, where i is the face 
  fac1 <- grid.facs[1]
  fac2 <- grid.facs[2]
  if (any(table(data[[fac1]])==0) | any(table(data[[fac2]])==0))
    stop("Some levels in the the supplied factors are not observed") 
  n1 <- length(levels(data[[fac1]]))
  n2 <- length(levels(data[[fac2]]))
  ns1 <- nsim*n1
  ns2 <- nsim*n2
  variog <- vector(mode = "list", length = 2)
  variog[[1]] <- vector(mode = "list", length = length(sections))
  names(variog[[1]]) <- sections
  variog[[1]] <- lapply(variog[[1]][sections], 
                        function(variog, nrow, ncol)
                          variog <- matrix(rep(0, nrow*ncol), nrow=nrow, ncol=ncol),
                        nrow=nsim, ncol=n1)
  variog[[2]] <- vector(mode = "list", length = length(sections))
  names(variog[[2]]) <- sections
  variog[[2]] <- lapply(variog[[2]][sections], 
                        function(variog, nrow, ncol)
                          variog <- matrix(rep(0, nrow*ncol), nrow=nrow, ncol=ncol),
                        nrow=nsim, ncol=n2)
  data[grid.facs] <- lapply(data[grid.facs], as.numfac) 
  residuals[grid.facs] <- lapply(residuals[grid.facs], as.numfac) 
  
  for (k in sections)
  { sect.dat <- data.frame(data[fac1], data[fac2])
    if (!is.null(fac.sec))
      sect.dat <- sect.dat[data[[fac.sec]] == k]
    names(sect.dat)[1:2] <- grid.facs
    for (i in 1:nsim)
    { if (is.null(fac.sec))
      sect.dat$res <- as.vector(resid$residuals[[(nvars+i-1)]])
      else
        sect.dat$res <- resid$residuals[data[[fac.sec]] == k, (nvars+i-1)]
      sim.var <- asreml::asreml.variogram(sect.dat)
      names(sim.var)[1:2] <- grid.facs
      variog[[1]][[k]][i,] <- sim.var[sim.var[[fac2]]==0,]$gamma
      variog[[2]][[k]][i,] <- sim.var[sim.var[[fac1]]==0,]$gamma
    }
  }
  
  #Do face plots for fac1 and fac2 in each section
  face.1 <- data.frame(matrix(nrow = 0, ncol=5))
  colnames(face.1) <- c("factor","observed","X2.5.","X50.","X97.5.")
  face.2 <- data.frame(matrix(nrow = 0, ncol=5))
  colnames(face.2) <- c("factor","observed","X2.5.","X50.","X97.5.")
  
  #Get the observed variogram and form the elements of the variofaces plot
  for (k in sections)
  { if (is.null(fac.sec))
    sect.dat <- data
    else
      sect.dat <- data[data[[fac.sec]]==k, -1]
    object.var <- asreml::asreml.variogram(sect.dat)
    names(object.var)[1:2] <- grid.facs
    
    #Form data.frame for current section
    face.1 <- rbind(face.1,
                    data.frame(factor = object.var[object.var[[fac2]]==0,][[fac1]],
                               observed = object.var[object.var[[fac2]]==0,]$gamma,   
                               t(apply(variog[[1]][[k]], 2, quantile, 
                                       probs=c(0.025, 0.50, 0.975),
                                       na.rm = TRUE)),
                               stringsAsFactors = FALSE))
    face.2 <- rbind(face.2,
                    data.frame(factor = object.var[object.var[[fac1]]==0,][[fac2]],   
                               observed = object.var[object.var[[fac1]]==0,]$gamma,   
                               t(apply(variog[[2]][[k]], 2, quantile, 
                                       probs=c(0.025, 0.50, 0.975),
                                       na.rm = TRUE)), 
                               stringsAsFactors = FALSE))
  }
  names(face.1)[1] <- fac1
  names(face.2)[1] <- fac2
  
  if (!is.null(fac.sec))
  { face.1[[fac.sec]] <- factor(rep(1:length(sections), each=(n1)), labels=sections)
    face.2[[fac.sec]] <- factor(rep(1:length(sections), each=(n2)), labels=sections)
  }
  
  #Do plots
  p <- ggplot(data=face.1) +
    theme_bw() +
    geom_line(aes_string(x= fac1, y = "observed"), size=1) +
    geom_point(aes_string(x= fac1, y = "observed"), size=3) +
    geom_line(aes_string(x=fac1, y = "X2.5."), colour="red") +
    geom_line(aes_string(x=fac1, y = "X97.5."), colour="red") +
    geom_line(aes_string(x=fac1, y = "X50."), colour="blue") +
    labs(title = paste("Variogram face of",restype,"for",fac1,sep=" "), 
         x = paste(fac1,"differences", sep=" "),y=NULL)
  if (!is.null(fac.sec))
    p <- p + facet_wrap(as.formula(paste("~",fac.sec)), ncol=2, scales="free")
  print(p)
  
  p <- ggplot(data=face.2) +
    theme_bw() +
    geom_line(aes_string(x= fac2, y = "observed"), size=1) +
    geom_point(aes_string(x= fac2, y = "observed"), size=3) +
    geom_line(aes_string(x=fac2, y = "X2.5."), colour="red") +
    geom_line(aes_string(x=fac2, y = "X97.5."), colour="red") +
    geom_line(aes_string(x=fac2, y = "X50."), colour="blue") +
    labs(title = paste("Variogram face of",restype,"for",fac2,sep=" "), 
         x = paste(fac2,"differences", sep=" "), y=NULL)
  if (!is.null(fac.sec))
    p <- p + facet_wrap(as.formula(paste("~",fac.sec)), ncol=2, scales="free")
  print(p)
  
  #return data frames containing the variogram values on which the plots are based
  invisible(list(face1 = face.1, face2 = face.2))
}
