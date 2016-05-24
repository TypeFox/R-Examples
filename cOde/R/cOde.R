#' Generate C code for a function and compile it
#' 
#' @param f Named character vector containing the right-hand sides of the ODE. You may use the key word
#' \code{time} in your equations for non-autonomous ODEs. 
#' @param forcings Character vector with the names of the forcings
#' @param outputs Named character vector for additional output variables, 
#' see arguments \code{nout} and \code{outnames} of \link[deSolve]{lsode}
#' @param jacobian Character, either "none" (no jacobian is computed), "full" (full jacobian 
#' is computed and written as a function into the C file) or "inz.lsodes" (only the non-zero elements
#' of the jacobian are determined, see \link[deSolve]{lsodes})
#' @param rootfunc Named character vector. The root function (see \link[deSolve]{lsoda}). Besides the
#' variable names (\code{names(f)}) also other symbols are allowed that are treated like new
#' parameters.
#' @param boundary data.frame with columns name, yini, yend specifying the boundary condition set-up. NULL if not a boundary value problem
#' @param compile Logical. If FALSE, only the C file is written
#' @param fcontrol Character, either \code{"nospline"} (default, forcings are handled by deSolve) or \code{"einspline"}
#' (forcings are handled as splines within the C code based on the einspline library).
#' @param nGridpoints Integer, defining the number of grid points between tmin and tmax where the ODE
#' is computed in any case. Indicates also the number of spline nodes if \code{fcontrol = "einspline"}.
#' @param precision Numeric. Only used when \code{fcontrol = "einspline"}.
#' @param modelname Character. The C file is generated in the working directory and is named <modelname>.c.
#' If \code{NULL}, a random name starting with ".f" is chosen, i.e. the file is hidden on a UNIX system.
#' @param verbose Print compiler output to R command line.
#' @details The function replaces variables by arrays \code{y[i]}, etc. and replaces "^" by pow() 
#' in order to have the correct C syntax. The file name of the C-File is derived from \code{f}. 
#' I.e. \code{funC(abc, ...} will generate a file abc.c in the current directory. 
#' Currently, only explicit ODE specification is supported, i.e. you need to have the right-hand sides of the ODE.
#' 
#' @return the name of the generated shared object file together with a number of attributes
#' @examples 
#' \dontrun{
#' # Exponential decay plus constant supply
#' f <- c(x = "-k*x + supply")
#' func <- funC(f, forcings = "supply")
#' 
#' # Example 2: root function
#' f <- c(A = "-k1*A + k2*B", B = "k1*A - k2*B")
#' rootfunc <- c(steadyState = "-k1*A + k2*B - tol")
#' 
#' func <- funC(f, rootfunc = rootfunc, modelname = "test")
#' 
#' yini <- c(A = 1, B = 2)
#' parms <- c(k1 = 1, k2 = 5, tol = 0.1)
#' times <- seq(0, 10, len = 100)
#' 
#' odeC(yini, times, func, parms)
#' }
#' @export
#' @importFrom utils packageVersion
funC <- function(f, forcings=NULL, outputs=NULL, 
                 jacobian=c("none", "full", "inz.lsodes", "jacvec.lsodes"), 
                 rootfunc = NULL, boundary= NULL, 
                 compile = TRUE, fcontrol = c("nospline", "einspline"),
                 nGridpoints = 500, precision=1e-5, modelname = NULL,
                 verbose = FALSE) {
  
  f <- unclass(f)
  constraints <- NULL # Might be an interesting option in the future
  myattr <- attributes(f)
  equations <- f
  
  if("names"%in%names(myattr)) myattr <- myattr[-which(names(myattr)=="names")]
  
  if(is.null(modelname)) modelname <- paste(c(".f", sample(c(letters, 0:9), 8, TRUE)), collapse="")
  dllname <- modelname
  filename <- paste0(dllname, ".c")
  
  ## If boundary conditions are given, sort for leftbc first
  if(!is.null(boundary)) {
    leftbc <- boundary$name[!is.na(boundary$yini)]
    f <- c(f[names(f)%in%leftbc], f[!names(f)%in%leftbc])
  }
  
  ## Check which kind of forcings are used
  fcontrol <- match.arg(fcontrol)
  forcings.t <- paste(forcings, "t", sep=".")
  forc.replace <- forcings
  forc.t.replace <- forcings.t
  if(fcontrol == "einspline" & !is.null(forcings)) {
    forc.replace <- paste0("x[", 1:length(forcings)-1, "]")
    forc.t.replace <- paste0("xdot[", 1:length(forcings.t)-1, "]")
  }
  
  ## Analyze f by parser
  
  variables <- names(f)
  symbols <- getSymbols(c(f, rootfunc, constraints))
  parameters <- symbols[!symbols%in%c(variables, forcings, names(constraints), names(rootfunc), "time")]
  jac <- NULL
  inz <- NULL
  
  jacobian <- match.arg(jacobian)
  if(jacobian != "none") jac  <- jacobianSymb(f)
  if(jacobian %in% c("inz.lsodes", "jacvec.lsodes")) {
    jac.matrix <- matrix(jac, length(f), length(f))
    inz <- apply(jac.matrix, 2, function(v) which(v != "0"))
    inz <- do.call(rbind, lapply(1:length(inz), function(j) if(length(inz[[j]]) > 0) cbind(i = inz[[j]], j = j)))
  }
  not.zero.jac <- which(jac != "0")
  
  dv <- length(variables)
  dp <- length(parameters)
  if(is.null(forcings)) di <- 0 else di <- length(forcings)
  if(is.null(constraints)) dc <- 0 else dc <- length(constraints)
  if(is.null(outputs)) do <- 0 else do <- length(outputs)
  if(is.null(rootfunc)) dr <- 0 else dr <- length(rootfunc)
  
  
  ## Replace powers and symbols to get correct C syntax
  CSyntax <- function(f) {
    f <- replaceOperation("^", "pow", f)
    f <- replaceSymbols(variables, paste0("y[", 1:length(variables)-1, "]"), f)
    f <- replaceSymbols(names(constraints), paste0("cons[", 1:dc-1, "]"), f)
    f <- replaceSymbols(forcings, forc.replace, f)
    f <- replaceSymbols(forcings.t, forc.t.replace, f)
    return(f)
  }
  
  f <- CSyntax(f)
  
  if(jacobian %in% c("full", "jacvec.lsodes")) 
    jac <- CSyntax(jac)
  
  if(!is.null(constraints)) 
    constraints <- CSyntax(constraints)
  
  if(!is.null(outputs)) 
    outputs <- CSyntax(outputs)
  
  if(!is.null(rootfunc))
    rootfunc <- CSyntax(rootfunc)

  
  
  
  ## ------------ write C code -------------
  
  mypath <- system.file(package="cOde")
  splinefile <- paste0("cat ", mypath,"/code/splineCreateEvaluate.c")
  includings <- c("#include <R.h>",
                  "#include <math.h>")
  
  if(fcontrol == "einspline") 
    includings <- c(includings, "#include <einspline/nubspline.h>")
  
  definitions <- paste0("#define ", c(parameters, paste0("y",0:(dv-1),"_0")), " parms[", 0:(dv+dp-1),"]")
  if(!is.null(forcings)) 
    definitions <- c(definitions, paste0("#define ", forcings, " forc[", 0:(di-1),"]"))
  
  
  sink(filename)
  cat("/** Code auto-generated by cOde", as.character(packageVersion("cOde")), "**/\n")
  cat(paste(includings, "\n"))
  cat("\n")
  cat(paste("static double parms[", dv+dp,"];\n", sep=""))
  cat(paste("static double forc[", di,"];\n", sep=""))
  cat(paste("static double cons[", dc,"];\n", sep=""))
  cat("static double range[2];\n")
  cat("\n")
  cat(paste("#define nGridpoints",nGridpoints,"\n"))
  cat(paste("#define nSplines", di, "\n"))
  cat(paste("#define precision", precision, "\n"))
  cat("\n")
  cat(paste(definitions, "\n"))
  cat("#define tmin range[0]\n")
  cat("#define tmax range[1]\n")
  cat("\n")
  if(fcontrol == "einspline") cat(paste(system(splinefile, intern=TRUE), "\n"))
  cat("\n")
  cat("void initmod(void (* odeparms)(int *, double *)) {\n")
  cat(paste("\t int N=", dv+dp,";\n",sep=""))
  cat("\t odeparms(&N, parms);\n")
  cat("}\n")
  cat("\n")
  cat("void initforc(void (* odeforcs)(int *, double *)) {\n")
  cat(paste("\t int N=", di,";\n",sep=""))
  cat("\t odeforcs(&N, forc);\n")
  cat("}\n")
  cat("\n")
  
  
  ## Derivative function
  
  cat("/** Derivatives (ODE system) **/\n")
  cat("void derivs (int *n, double *t, double *y, double *ydot, double *RPAR, int *IPAR) {\n")
  cat("\n")
  if(fcontrol == "einspline") {
    cat("\t double x[nSplines];\n")
    cat("\t double xdot[nSplines];\n")
    cat("\t evaluateSplines(t, x, xdot);\n")
  }
  cat("\t double time = *t;\n")
  if(!is.null(constraints)) 
    cat(paste("\t cons[", 1:dc-1, "] = ", constraints, ";\n", sep=""))
  
  cat("\n")
  #if(length(reductions)>0) cat(paste("\t double ", reductions, ";\n", sep=""))
  cat(paste("\t ydot[", 0:(dv-1),"] = ", f,";\n", sep=""))
  cat("\n")
  
  # Return forcings and other outputs (only for IVP)
  if(is.null(boundary)) {
    if(di > 0){
      cat(paste0("\t RPAR[", 0:(di-1),"] = ", forc.replace,";\n"))
    }
    if(do > 0){
      cat("\t for(int i= ",di,"; i < ",do+di,"; ++i) RPAR[i] = 0;\n")
      non.zero.outputs <- which(outputs != "0")
      for(i in non.zero.outputs) 
        cat(paste0("\t RPAR[", di + i - 1, "] = ", outputs[i], ";\n")) 
    }
  }
    
  cat("}\n")
  cat("\n")
  
  ## Jacobian of deriv
  if(jacobian == "full") {
    cat("/** Jacobian of the ODE system **/\n")
    cat("void jacobian (int *n, double *t, double *y, double * df, double *RPAR, int *IPAR) {\n")
    cat("\n")
    cat("\n")
    if(fcontrol == "einspline") {
      cat("\t double x[nSplines];\n")
      cat("\t double xdot[nSplines];\n")
      cat("\t evaluateSplines(t, x, xdot);\n")
    }
    cat("\n")
    cat("double time = *t;\n")
    cat("\t int i;\n")
    cat("for(i=0; i<(*n) * (*n); i++) df[i] = 0.;\n")
    cat(paste("\t df[", not.zero.jac-1,"] = ", jac[not.zero.jac],";\n", sep=""))
    cat("\n")
    cat("}\n")
    cat("\n")
  }
  
  ## Jacvec of deriv
  if(jacobian == "jacvec.lsodes") {
    vecs <- lapply(1:dv, function(i) matrix(jac, ncol=dv, nrow=dv)[,i])
    not.zero.vec <- lapply(vecs, function(v) which(v != "0"))
    not.zero.columns <- which(sapply(not.zero.vec, function(v) length(v) > 0))
    cat("/** Jacobian vector of the ODE system **/\n")
    cat("void jacvec (int *neq, double *t, double *y, int *j, int *ian, int *jan, double *pdj, double *yout, int *iout) {\n")
    cat("\n")
    cat("\n")
    if(fcontrol == "einspline") {
      cat("\t double x[nSplines];\n")
      cat("\t double xdot[nSplines];\n")
      cat("\t evaluateSplines(t, x, xdot);\n")
    }
    cat("\t double time = *t;\n")
    cat("\t int i;\n")
    cat("\t for(i=0; i<*neq; i++) pdj[i] = 0.;\n")
    
    j <- not.zero.columns[1]
    cat(paste("\t if(*j ==", j, ") {\n"))
    cat(paste("\t pdj[", not.zero.vec[[j]]-1,"] = ", vecs[[j]][not.zero.vec[[j]]],";\n", sep=""))
    cat("\t }\n")
    for(j in not.zero.columns[-1]) {
      cat(paste("\t else if(*j ==", j, ") {\n"))
      cat(paste("\t pdj[", not.zero.vec[[j]]-1,"] = ", vecs[[j]][not.zero.vec[[j]]],";\n", sep=""))
      cat("\t }\n")
    }
    
    cat("\n")
    cat("}\n")
    cat("\n")
    
  }
  
  if(!is.null(rootfunc)) {
    
    cat("/** Root function **/\n")
    cat("void myroot(int *neq, double *t, double *y, int *ng, double *gout, double *out, int *ip ) {\n")
    cat("\n")
    cat("\n")
    if(fcontrol == "einspline") {
      cat("\t double x[nSplines];\n")
      cat("\t double xdot[nSplines];\n")
      cat("\t evaluateSplines(t, x, xdot);\n")
    }
    cat("\t double time = *t;\n")
    cat(paste("\t gout[", 0:(dr-1),"] = ", rootfunc,";\n", sep=""))
    cat("\n")
    cat("}\n")
    
  }
  
  if(!is.null(boundary)) {
    
    ## Check length of boundary conditions
    nbc <- length(which(!is.na(c(boundary$yini, boundary$yend))))
    if(nbc != dv) {
      sink()
      warning("Number of boundary conditions not correct\n")
      return()
    }
    
    boundary <- boundary[match(variables, boundary$name),]
    
    leftbc <- which(!is.na(boundary$yini))
    rightbc <- which(!is.na(boundary$yend))
    myorder <- c(leftbc, rightbc)
    
    ## Boundary Condition (for compatibility with bvpSolve)
    
    cat("/** Boundary Conditions **/\n")
    cat("void gsub(int *i, int *n, double *z, double *g, double *RPAR, int *IPAR) {\n")
    cat("\n")
    cat(paste("\t if (*i==", 1,") *g=z[", myorder[1]-1, "]-y", 0, "_0;\n", sep=""))
    if(dv>1) cat(paste("\t else if (*i==", 2:dv,") *g=z[", myorder[-1]-1, "]-y", 2:dv-1, "_0;\n", sep=""))
    cat("\n")
    cat("}\n")
    cat("\n")
    
    ## Jacobian of Boundary Condition (for compatibility with bvpSolve)
    
    cat("/** Jacobian of the Boundary Conditions **/\n")
    cat("void dgsub(int *i, int *n, double *z, double *dg, double *RPAR, int *IPAR) {\n")
    cat("\n")
    cat("\t int j;\n")
    cat("\t for (j = 0; j< *n; j++) dg[j] = 0;\n")
    
    cat(paste("\t if (*i==", 1,") dg[", myorder[1]-1, "] = 1.;\n", sep=""))
    if(dv>1) cat(paste("\t else if (*i==", 2:dv,") dg[", myorder[-1]-1, "]=1.;\n", sep=""))
    cat("\n")
    cat("}\n")
  }
  
  
  
  sink()
  
  ## ----------- compile C code and load shared object file---------
  
  if(compile & fcontrol == "nospline") 
    shlibOut <- system(paste0(R.home(component="bin"), "/R CMD SHLIB ", filename), intern = TRUE)
  if(compile & fcontrol == "einspline") 
    shlibOut <- system(paste0(R.home(component="bin"), "/R CMD SHLIB ", filename, " -leinspline"), intern = TRUE)
  if (verbose) {
    cat(shlibOut)
  }
  
  .so <- .Platform$dynlib.ext
  soExists <- file.exists(paste0(dllname, .so))
  if(soExists) {
    try(dyn.unload(paste0(dllname, .so)), silent = TRUE)
    dyn.load(paste0(dllname, .so))
  }
  
  ## ----------- function return -----------
  
  
  f <- dllname
  attributes(f) <- c(attributes(f), myattr)
  
  attr(f, "equations") <- equations
  attr(f, "variables") <- variables
  attr(f, "parameters") <- parameters
  attr(f, "forcings") <- forcings
  attr(f, "outputs") <- outputs
  attr(f, "jacobian") <- jacobian
  attr(f, "inz") <- inz
  attr(f, "boundary") <- boundary
  attr(f, "rootfunc") <- rootfunc
  attr(f, "nGridpoints") <- nGridpoints
  attr(f, "fcontrol") <- fcontrol
  
  return(f)
  
}



#' Dynamically load DLL with automatic unloading of conflicting DLLs
#' 
#' @param func result from funC(), contains the information about the DLL name to be loaded
#' @param cfunction character, denoting the C function name.
#' @details If the C function name is already part of another loaded DLL, the corresponding DLL is
#' unloaded and the desired func DLL is loaded instead.
#' @export
loadDLL <- function(func, cfunction="derivs") {
 
  .so <- .Platform$dynlib.ext
  checkDLL <- try(getNativeSymbolInfo(cfunction), silent=TRUE)
  if(inherits(checkDLL, "try-error")) {
    dyn.load(paste0(func, .so))
    cat("Shared object is loaded and ready to use\n")
  } else if((checkDLL$package)[[1]] != func) {
    test <- try(dyn.unload(paste0((checkDLL$package)[[1]], .so)), silent = TRUE)
    if(!inherits(test, "try-error")) 
      warning("Conflicting shared object was unloaded and new one is loaded")
    dyn.load(paste0(func, .so))
  }
  
}

#' Generate interpolation spline for the forcings and write into list of matrices
#' 
#' @param func result from funC()
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric)
#' @return list of matrices with time points and values assigned to the forcings interface of deSolve
#' @details Splines are generated for each name in forcings and both, function value and first
#' derivative are evaluated at the time points of the data frame.
#' @examples
#' f <- c(x = "-k*x + a - b")
#' func <- funC(f, forcings = c("a", "b"))
#' forcData <- rbind(
#'   data.frame(name = "a", time = c(0, 1, 10), value = c(0, 5, 2)),
#'   data.frame(name = "b", time = c(0, 5, 10), value = c(1, 3, 6)))
#' forc <- setForcings(func, forcData) 
#' @export
setForcings <- function(func, forcings) {
  
  #loadDLL(func)
  
  inputs <- attr(func, "forcings")
  fcontrol <- attr(func, "fcontrol")
  nGridpoints <- attr(func, "nGridpoints")
  trange <- range(forcings$time)
  tspan <- seq(trange[1], trange[2], len=nGridpoints)
  
  times <- NULL
  values <- NULL
  
  
  out <- do.call(c, lapply(inputs, function(i) {
    
    t <- forcings[forcings$name == i, "time"]
    x <- forcings[forcings$name == i, "value"]
    
    if(fcontrol == "nospline") {
      mat <- list(cbind(t, x))
      names(mat) <- i
    }
    if(fcontrol == "einspline") {
      myfun <- splinefun(t, x)
      out <- myfun(tspan)
      mat <- list(cbind(tspan, out))
      names(mat) <- i
    }
 
    return(mat)
    
  }))
  
  
  if(fcontrol == "einspline") {
   
    times <- do.call(c, lapply(out, function(o) o[,1]))
    values <- do.call(c, lapply(out, function(o) o[,2]))
    cfunc <- getNativeSymbolInfo("createSplines")
    .C(cfunc, as.double(times), as.double(values))
    
  }
  
  return(out)
  
}


#' Interface to ode()
#' 
#' @param y named vector of type numeric. Initial values for the integration
#' @param times vector of type numeric. Integration times
#' @param func return value from funC()
#' @param parms named vector of type numeric. 
#' @param ... further arguments going to \code{ode()}
#' @details See deSolve-package for a full description of possible arguments
#' @return matrix with times and states
#' @example inst/examples/example1.R
#' @export
odeC <- function(y, times, func, parms, ...) {
  
  nGridpoints <- attr(func, "nGridpoints")
  times.inner <- seq(min(c(times, 0)), max(times), len=nGridpoints)
  times.inner <- sort(unique(c(times, times.inner)))
  which.times <- match(times, times.inner)
  yout <- c(attr(func, "forcings"), names(attr(func, "outputs")))
  
  y <- y[attr(func, "variables")]
  parms <- parms[attr(func, "parameters")]
  parms <- c(parms, rep(0, length(y)))
  
  arglist <- list(y = y, times = times.inner, func = "derivs", parms = parms, dllname = func, initfunc = "initmod")
  
  
  if (attr(func, "jacobian") == "full")
    arglist <- c(arglist, list(jacfunc = "jacobian"))
    
  if (attr(func, "jacobian") == "inz.lsodes") {
    inz <- attr(func, "inz")
    lrw <- 20 + 3*dim(inz)[1] + 20*length(y)
    arglist <- c(arglist, list(sparsetype = "sparseusr", inz = inz, lrw = lrw))
  }
  
  if (attr(func, "jacobian") == "jacvec.lsodes") {
    inz <- attr(func, "inz")
    arglist <- c(arglist, list(sparsetype = "sparseusr", jacvec = "jacvec", inz = inz))
  }
    
  if (!is.null(attr(func, "forcings")) & attr(func, "fcontrol") == "nospline") 
    arglist <- c(arglist, list(initforc = "initforc"))
 
  if(!is.null(attr(func, "rootfunc")))
    arglist <- c(arglist, list(rootfunc = "myroot", nroot = length(attr(func, "rootfunc"))))
   
  if (!is.null(yout)) {
    arglist <- c(arglist, list(nout = length(yout), outnames = yout))
  }
  
  # Replace arguments and add new ones
  moreargs <- list(...)
  if(any(names(moreargs)=="forcings") & attr(func, "fcontrol") == "einspline") 
    moreargs <- moreargs[-which(names(moreargs)=="forcings")]
  
  i <- match(names(moreargs), names(arglist))
  is.overlap <- which(!is.na(i))
  is.new <- which(is.na(i))
  arglist[i[is.overlap]] <- moreargs[is.overlap]
  arglist <- c(arglist, moreargs[is.new])
  
  #loadDLL(func)
  
  out <- do.call(deSolve::ode, arglist)
  out.index <- unique(c(which.times[which.times <= nrow(out)], nrow(out)))
  out <- matrix(out[out.index, ], nrow = length(out.index), dimnames = list(NULL, colnames(out)))
  
  return(out)
  
  
}


#' Interface to bvptwp()
#' 
#' 
#' @param yini named vector of type numeric. Initial values to be overwritten.
#' @param x vector of type numeric. Integration times
#' @param func return value from funC() with a boundary argument. 
#' @param yend named vector of type numeric. End values to be overwritten.
#' @param parms named vector of type numeric. The dynamic parameters.
#' @param xguess vector of type numeric, the x values
#' @param yguess matrix with as many rows as variables and columns as x values
#' @param ... further arguments going to \code{bvptwp()}
#' @details See bvpSolve-package for a full description of possible arguments
#' @return matrix with times and states
#' @example inst/examples/example4.R
#' @export
bvptwpC <- function(yini=NULL, x, func, yend=NULL, parms, xguess=NULL, yguess=NULL,  ...) {
  
  #loadDLL(func)
  
  dynpar <- parms[attr(func, "parameters")]
  boundary <- attr(func, "boundary")
  leftbc <- boundary$name[!is.na(boundary$yini)]
  rightbc <- boundary$name[!is.na(boundary$yend)]
  
  ## Fill yini/yend with values from func. If yini/yend are given,
  ## set their values.
  bini <- boundary$yini
  names(bini) <- boundary$name
  bini <- bini[!is.na(bini)]
  
  bend <- boundary$yend
  names(bend) <- boundary$name
  bend <- bend[!is.na(bend)]
  
  if(!is.null(yini)) bini[names(yini)] <- yini
  if(!is.null(yend)) bend[names(yend)] <- yend
  if(!is.null(attr(func, "forcings")) & attr(func, "fcontrol") == "nospline") 
    initforc <- "initforc"
  else
    initforc <- NULL
  
  
  
  posbound <- c(rep(min(x), length(bini)), rep(max(x), length(bend)))
  
  
  statepars <- c(bini, bend)
  newparms <- c(dynpar, statepars)
  
  moreargs <- list(...)
  if(any(names(moreargs)=="forcings") & attr(func, "fcontrol") == "einspline") 
    moreargs <- moreargs[-which(names(moreargs)=="forcings")]

  print(initforc)
   
  out <- do.call(bvpSolve::bvptwp, c(moreargs, list(
    x = x, parms = newparms, xguess = xguess, yguess = yguess, posbound=posbound,
    func = "derivs", jacfunc = "jacobian", bound = "gsub", jacbound = "dgsub", 
    initfunc = "initmod", initforc = initforc,
    dllname = func,
    ncomp = length(statepars)
  )))
  

  colnames(out) <- c("x", attr(func, "variables"))
  
  return(out)
  
  
}

