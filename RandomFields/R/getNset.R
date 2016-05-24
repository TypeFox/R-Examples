
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



### get&set functions using the C interface
### and RFgetNset functions

summary.RFopt <- function(object, ...) {
  object <- lapply(object, function(z) z[order(names(z))])
  object <- object[c(1, 1 + order(names(object[-1])))]
  class(object) <- "summary.RFopt"
  object
}


print.summary.RFopt <- function(x, ...) {
  str(x, give.attr=FALSE, ...) #
  invisible(x)
}

print.RFopt <- function(x, ...) {
  print.summary.RFopt(summary.RFopt(x, ...)) #
  invisible(x)
}

summary.RFoptElmnt <- function(object, ...) {
  object <- object[order(names(object))]
  class(object) <- "summary.RFoptElmt"
  object
}

print.summary.RFoptElmnt <- function(x, ...) {
  str(x, give.attr=FALSE, ...) #
  invisible(x)
}

print.RFoptElmnt <- function(x, ...) {
  print.summary.RFoptElmnt(summary.RFoptElmnt(x, ...)) #
  invisible(x)
}


RFoptions <- function(..., no.readonly=TRUE) {
##  on.exit(.C("RelaxUnknownRFoption", FALSE))
##  .C("RelaxUnknownRFoption", TRUE)
  opt <- lapply(.External("RFoptions", ...),
                function(x) {
                  class(x) <- "RFoptElmnt"
                  x
                })
  if (length(opt)!=0) {      
    class(opt) <-  "RFopt"
    if (!no.readonly) {
      opt$readonly <- list(covmaxchar=MAXCHAR,
                           covnr=GetCurrentNrOfModels(FALSE),
                           maxdim=c(cov=MAXCOVDIM, mle=MAXMLEDIM,
                               simu=MAXSIMUDIM, ce=MAXCEDIM, tbm=MAXTBMSPDIM,
                               mpp=MAXMPPDIM, hyper=MAXHYPERDIM,
                               nug=MAXNUGGETDIM, vario=MAXVARIODIM),
                           maxmodels=MAXFIELDS,
                           methodmaxchar=METHODMAXCHAR,
                           methodnr=Forbidden
                           )
    }
  }
  if (length(opt)==0) invisible(opt) else opt
}

GetCurrentNrOfModels <- function(init=TRUE) {
  res<- .C("GetCurrentNrOfModels", as.integer(init), nr=as.integer(1))$nr
}

internal.rfoptions <- function(..., REGISTER=FALSE, COVREGISTER=as.integer(NA),
                               RELAX=FALSE){
 # Print(list(...))
  RFopt <- list()
  RFopt[[1]] <- .External("RFoptions")
  if (is.logical(REGISTER)) {
    REGISTER <- if (REGISTER) RFopt[[1]]$registers$register else as.integer(NA)
  }
  RFopt[[1]]$general$storing <-
    c(RFopt[[1]]$general$storing, REGISTER, COVREGISTER)
  l <- list(asList=FALSE, ...)
#  l <- list(...)
  if (length(l) > 1) { ## wegen asList; sonst  > 0 !!
    storing <- (substr(names(l), 1, 3) == "sto" |
                substr(names(l), 1, 9) == "general.sto")
    if (any(storing)) last <- rev(which(storing))[1]
    if (any(storing) && !l[[last]]) {
      for (p in which(storing)) l[[p]] <- c(FALSE, REGISTER, COVREGISTER)
    }
    on.exit(.C("RelaxUnknownRFoption", FALSE))
    .C("RelaxUnknownRFoption", RELAX)
    .External("RFoptions", LIST=l)
    RFopt[[2]] <- .External("RFoptions")    
  } else {
    RFopt[[2]] <- RFopt[[1]]
  }
  return(RFopt)
}

xylabs <- function(x, y, T=NULL, units=NULL) {
  if (is.null(units)) units <- RFoptions()$coords$coordunits
  xlab <- if (is.null(x)) NULL
          else if (units[1]=="") x else paste(x, " [", units[1], "]", sep="")
  ylab <- if (length(y)==0) NULL
          else if (units[2]=="") y else paste(y, " [", units[2], "]", sep="")
  Tlab <- if (length(T)==0) NULL
          else if (units[3]=="") T else paste(T, " [", units[3], "]", sep="")
  return (list(xlab=xlab, ylab=ylab, Tlab=Tlab))
}

add.units <- function(x,  units=NULL) {
    if (is.null(x)) return(NULL)
  if (is.null(units)) units <- RFoptions()$coords$varunits
  return(ifelse(units=="", x, paste(x, " [", units, "]", sep="")))
}



InitModel <- function(reg, model, dim, NAOK=FALSE){ # ok
  for (y in list(double(0), matrix(nrow=dim, ncol=3, as.double(1:3)))) {
    vdim <- try(.Call("Init",
                      MODEL_USER,
                      model,
                      list(x=matrix(nrow=dim, ncol=3, as.double(1:3)), #0 nur dummies
                           y=y, #1 y 
                           as.double(0), #2 T
                           FALSE, #3 grid
                           as.integer(dim), #4 spatdim
                           FALSE, #5 Zeit
                           FALSE #6 distances
                           ),
                      NAOK=NAOK, # ok
                      PACKAGE="RandomFields"), silent=TRUE)
    if (is.numeric(vdim)) return(vdim)
    msg <- strsplit(vdim[[1]], "\n")[[1]][2]
    if (RFoptions()$general$printlevel >= PL_ERRORS) cat(msg, "\n")
  }
  stop(msg)
  ##  stop("model could not be initialized")
}


resolve.register <- function(register){
  if (missing(register) || length(register) == 0) {
    register <- .C("GetCurrentRegister", reg=integer(1))$reg
    if (register < 0) stop("model that has been used right now cannot be determined or no model has been used up to now")
  }
  if (!is.numeric(register)) {
 #   register <- deparse(substitute(register))   
    register <-
      switch(register,
             "RFcov" = MODEL_USER,
             "RFcovmatrix" = MODEL_USER,
             "RFfctn" = MODEL_USER,
             "RFpseudovariogram" =  MODEL_USER,
             "RFvariogram" = MODEL_USER,
             
             "RFdistr" = MODEL_USER,
             "RFddistr" = MODEL_USER,
             "RFpdistr" = MODEL_USER,
             "RFqdistr" = MODEL_USER,
             "RFrdistr" = MODEL_USER,
             
             "RFcrossvalidate" =  MODEL_MLE,
             "RFfit" =  MODEL_MLE,
             "RFgui" =  MODEL_GUI,
             "RFinterpolate" =  MODEL_KRIGE,
             "RFratiotest" =  MODEL_MLE,
             "RFsimulate" =  RFoptions()$registers$register,
             stop("register unknown")
             )
  }
  stopifnot(is.numeric(register))
  if (register < 0) stop("'register' does not have a non-negative value.")
  return(register)
}



print_RFgetModelInfo <- function(x, max.level=99,
                                 give.attr=FALSE, ...) {
  if (is.null(attr(x, "level"))) {
    y <- x
    y$minmax <- NULL
    str(y, give.attr=FALSE) #
    types <- sort(unique(x$minmax$type))
    if (length(types) > 0) {
      cat(" $ minmax: \n")
      print(x$minmax, justify="right") #
      cat(" pmin/pmax : bound usually met in practice\n",
          "type\n",
          paste("   type =", formatC(types, width=2), ":",
                TYPEOF_PARAM_NAMES[types + 1], "\n"),
          "NAN/bayes : internal\n",
          "min/max   : mathematically valid interval for the parameter\n",
          "omin/omax : whether the interval is open to the left/right\n",
          "col/row   : size of parameter\n")
    }
  } else {
    str(object = x,  max.level=max.level, give.attr=give.attr, ...) #
  }
  x
}

print.RFgetModelInfo <- function(x, ...) {
  print_RFgetModelInfo(x, ...) 
}


RFgetModelInfo <- function(...) {
  x <- list(...)
  if (length(x) > 0 &&
      (is(x[[1]], "RMmodel") || is(x[[1]], "list"))) RFgetModelInfo_model(...)
  else RFgetModelInfo_register(...)
}

RFgetModelInfo_model <- function(model, dim = 1, Time = FALSE,
                                 kernel = FALSE, exclude_trend = TRUE, ...) {
  Reg <- MODEL_AUX
  relax <- isFormulaModel(model)
  RFoptOld <- internal.rfoptions(..., RELAX=relax)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  if (!exclude_trend) {
    stop("'exclude_trend = FALSE' not programmed yet.")
  }

#  if (FALSE)    {
#
#    library(RandomFields, lib="~/TMP"); RFoptions(print = 1)
#    for (f in dir(pattern="*.R")) {print(f); source(f)}
#
 #   model <- RMexp(scale=4, var=2) + RMnugget(var=3) + RMtrend(mean=1)
#    z <- RFsimulate(model, 1:4, storing=TRUE)
 #   RFgetModel(show.call=!FALSE)
#
#  #  (z <- RFgetModelInfo(RMwhittle(scale=NA, nu=NA) + NA))
#   #  RFgetModelInfo()
#    
#  }

  dim <- as.integer(dim)
  intern <- try(.Call("SetAndGetModelInfo", Reg,
                          list("Dummy", PrepareModel2(model)),
                          dim, FALSE, as.logical(kernel), as.logical(Time), dim,
                          as.integer(10), ## ehemals RFoptions(short=10)
                          TRUE, as.logical(exclude_trend),
                          PACKAGE="RandomFields"))
  if (class(intern) == "try-error") return(0)
  intern$effect <- intern$xdimOZ <- intern$matrix.indep.of.x <- NULL

  intern$minmax <- as.data.frame(intern$minmax)
  storage.mode(intern$minmax$omin) <- "logical"
  storage.mode(intern$minmax$omax) <- "logical"
  storage.mode(intern$minmax$bayes) <- "logical"
  storage.mode(intern$minmax$NAN) <- "integer"
  storage.mode(intern$minmax$col) <- "integer"
  storage.mode(intern$minmax$row) <- "integer"
 
  
  class(intern) <- "RFgetModelInfo"
  intern
}

RFgetModelInfo_register<-
  function(register, level=1, 
           spConform=RFoptions()$general$spConform,
           which.submodels = c("user", "internal",
               "call+user", "call+internal",
               "user.but.once", "internal.but.once",
               "user.but.once+jump", "internal.but.once+jump", "all"),
           modelname=NULL) {  
  ## positive values refer the covariance models in the registers
  ##define MODEL_USER : 0  /* for user call of Covariance etc. */
  ##define MODEL_SIMU : 1  /* for GaussRF etc */ 
  ##define MODEL_INTERN  : 2 /* for kriging, etc; internal call of covariance */
  ##define MODEL_MLE  : 3
  ##define MODEL_BOUNDS  4 : - /* MLE, lower, upper */
  ## level + 10: auch die call fctn !
  ## [ in RF.h we have the modulus minus 1 ]
    
  register <- resolve.register(if (missing(register)) NULL else
                               if (is.numeric(register)) register else
                               deparse(substitute(register)))
#  Print(which.submodels, as.list(RFgetModelInfo)$which.submodels[-1])  
  w <- (if (!hasArg("which.submodels")) 0
        else pmatch(match.arg(which.submodels),
                    as.list(RFgetModelInfo_register)$which.submodels[-1]) - 1)
   if (is.na(w) || length(w) ==0)
     stop("the value for 'which.submodels' does not match. Note that the definitoin has been changed in version 3.0.70")

  cov <- .Call("GetExtModelInfo", as.integer(register), as.integer(level),
               as.integer(spConform),
               as.integer(w),
               PACKAGE="RandomFields")

  if (!is.null(modelname)) {
    cov <- search.model.name(cov, modelname, 0)
  }
  class(cov) <- "RFgetModelInfo"
  attr(cov, "level") <- level
 
  return(cov)
}


RFgetModel <- function(register, explicite.natscale, show.call=FALSE)  {
  register <- resolve.register(if (missing(register)) NULL else
                               if (is.numeric(register)) register else
                               deparse(substitute(register)))
  modus <- (if (missing(explicite.natscale)) GETMODEL_AS_SAVED else
            if (explicite.natscale)  GETMODEL_DEL_NATSC else
            GETMODEL_SOLVE_NATSC)

  which <- if (is.logical(show.call)) {
    if (show.call) "call+user" else "user"
  } else {
    show.call
  }
  
  m <- GetModel(register=register, modus=modus, which.submodels = which)
  class(m) <-  "RM_model"
  m
}
           
           
GetModel <- function(register, modus=GETMODEL_DEL_NATSC,
                     spConform=RFoptions()$general$spConform,      
                     which.submodels = c("user", "internal",
                         "call+user", "call+internal",
                         "user.but.once", "internal.but.once",
                         "user.but.once+jump", "internal.but.once+jump", "all"),
                     do.notreturnparam=FALSE, solve_random = FALSE){
  
  ## modus:
  ##  AS_SAVED : Modell wie gespeichert
  ##  DEL_NATSC : Modell unter Annahme PracticalRange>0 (natsc werden geloescht)
  ##  SOLVE_NATSC : natscale soweit wie moeglich zusammengezogen (natsc werden
  ##               drauf multipliziert; Rest wie gespeichert)
  ##  DEL_MLE : nur natscale_MLE werden geloescht
  ##  SOLVE_MLE : nur natscale_MLE  zusammengezogen (natsc werden
  ##               drauf multipliziert; Rest wie gespeichert)
  ## 

  ## do.notreturnparam : if true, then also parameters with the flag
  ##                      DONOTRETURN are returned

  ## spConform : only the names of the models
  ##              if > 1 than also "+" is forced to be "RMplus"

  
  register <- resolve.register(if (missing(register)) NULL else
                               if (is.numeric(register))  register else
                               deparse(substitute(register)))
   w<-pmatch(match.arg(which.submodels),
            as.list(GetModel)$which.submodels[-1]) - 1
  if (missing(register)) register <- 0
  
  model <- .Call("GetModel", as.integer(register), as.integer(modus),
                 as.integer(spConform),
                 as.integer(w),
                 as.logical(solve_random),
                 as.integer(do.notreturnparam),
                 PACKAGE="RandomFields")
  return(model)
}

GetModelRegister <- function(name) { ## obsolete
  stopifnot(is.character(name))
  return(as.integer(.C("GetModelRegister", name, integer(1),
                       PACKAGE="RandomFields")[[2]]))
}


RFgetModelNames <- function(type = RC_TYPENAMES, domain = RC_DOMAIN_NAMES,
                            isotropy = RC_ISONAMES, operator = c(TRUE, FALSE),
                            monotone = RC_MONOTONE_NAMES,
                            implied_monotonicities = length(monotone) == 1,
                            finiterange = c(TRUE, FALSE, NA),
                            valid.in.dim = c(1, Inf), 
                            vdim = c(1, 5),
                            group.by,
                            simpleArguments = FALSE,
                            internal,
                            newnames
                            ){ #, .internal=FALSE) {

  group.names <- c("type", "domain", "isotropy", "operator",
                   "monotone", "finiterange", "valid.in.dim", "vdim")

  if (hasArg(internal)) {
    return(PrintModelList(operators=operator, internal = internal,
                          newstyle=missing(newnames) || newnames))
  }
  if (!missing(newnames) && !newnames) {
    if (hasArg(type) || hasArg(domain) || hasArg(isotropy) || hasArg(operator)
        || hasArg(monotone) || hasArg(finiterange) || hasArg(valid.in.dim)
        || hasArg(vdim) || hasArg(group.by))
      stop("use 'newnames=FALSE' without further parameters or in combination with 'internal'")
    return (.Call("GetAllModelNames", PACKAGE="RandomFields"))
  }
  

  if (!(length(valid.in.dim) %in% 1:2)) stop("'valid.in.dim' has wrong size.")
  if (length(valid.in.dim) == 1) valid.in.dim <- c(valid.in.dim, Inf)
  
  if (!(length(vdim) %in% 1:2)) stop("'vdim' has wrong size.")
  if (length(vdim) == 1) vdim <- rep(vdim, 2)

  debug <- !TRUE
  if (hasArg(type)) type <- TYPENAMES[pmatch(type, TYPENAMES)]
  if (hasArg(domain)) domain <- DOMAIN_NAMES[pmatch(domain, DOMAIN_NAMES)]
  if (hasArg(isotropy)) isotropy <- ISONAMES[pmatch(isotropy, ISONAMES)]
  if (hasArg(monotone)) monotone <- MONOTONE_NAMES[pmatch(monotone,
                                                           MONOTONE_NAMES)]
  if (any(is.na(pmatch(type, TYPENAMES))))
    stop(paste("'", type, "'", " is not a valid category", sep=""))
  if (any(is.na(pmatch(domain, DOMAIN_NAMES))))
    stop(paste("'", domain, "'", " is not a valid category", sep=""))
  if (any(is.na(pmatch(isotropy, ISONAMES))))
    stop(paste("'", isotropy, "'", " is not a valid category", sep=""))
  if (any(is.na(pmatch(monotone, MONOTONE_NAMES))))
    stop(paste("'", monotone, "'", " is not a valid category", sep=""))


  if (missing(group.by)) {
    if (length(type) == 1) {
      if (type == TYPENAMES[TcfType + 1]) type <- c(type, "undefined") # to do
      else if (type == TYPENAMES[PosDefType + 1])
        type <- c(TYPENAMES[TcfType + 1], type, "undefined") # to do
      else if (type == TYPENAMES[VariogramType + 1])
        type <- c(TYPENAMES[c(TcfType, PosDefType) + 1], type,"undefined")#to do
     }
     if (!hasArg("group.by")) group.by <- if (length(type) == 1) NULL else 'type'
  }
 
  if (length(group.by) > 0) {
    group.idx <- pmatch(group.by, group.names)
    if (any(is.na(group.idx)))
      stop("'group.by' can be equal to '",
           paste(group.names, collapse="', '"), "'")
    group.by <- group.names[group.idx]
    group.idx <- group.idx[1]
  }

  
  if (group <- !is.null(group.by)) {  
    FUN <- function(string){
      args <- list(type=type, domain=domain, isotropy=isotropy,
                   operator=operator, monotone=monotone,
                   implied_monotonicities = implied_monotonicities &&
                   group.by[1] != "monotone",
                   finiterange=finiterange, valid.in.dim=valid.in.dim,
                   vdim=vdim,
                   group.by =
                     if (group && length(group.by) > 1) group.by[-1] else NULL
                   )
      args[[group.idx]] <- string
      list(do.call("RFgetModelNames", args))
    }
    li <- sapply(get(group.by[1]), FUN=FUN)
    if (is.null(names(li)))
      names(li) <- paste(group.by[1], get(group.by[1]), sep="=")
    length <- unlist(lapply(li, FUN=length))
    li <- li[length>0]
    return(li)
  } # matches  if (hasArg(group.by)) {
  
  if (implied_monotonicities) {
    mon <- MONOTONE_NAMES[-1 : (MISMATCH - 1)]
    for (i in MONOTONE:NORMAL_MIXTURE)
      if (mon[i] %in% monotone) monotone <- c(monotone, mon[i + 1])
    if (mon[BERNSTEIN] %in% monotone) monotone <- c(monotone, mon[MONOTONE])
    monotone <- unique(monotone)
  }

  internals <-  c("RMcovariate", "RMfixcov")
  envir <- as.environment("package:RandomFields")
  ls.RFmg <- ls(envir=envir)
  ls.RFmg <- c(ls.RFmg[substr(ls.RFmg, 1, 1) == "R"],
               paste("i", internals, sep=""))
  idx <- logical(len <- length(ls.RFmg))
  for (i in 1:len) {
    fun <-  do.call(":::", list("RandomFields", ls.RFmg[i]))  ##get(ls.RFmg[i], envir=envir)
    idx[i] <- is.function(fun) && is(fun, class2="RMmodelgenerator")
    
    if (!idx[i]) next
       idx[i] <- (!all(is.na(pmatch(fun["type"], type, 
                                        duplicates.ok=TRUE) &
                           pmatch(fun["isotropy"], isotropy, 
                                  duplicates.ok=TRUE))) &&
          #      !any(pt[!is.na(pt)] %in% pi[!is.na(pi)]) &&
                !all(is.na(pmatch(domain, fun["domain"]))) &&
                fun["operator"] %in% operator &&
                !all(is.na(pmatch(monotone, fun["monotone"]))) &&
                (!simpleArguments || fun["simpleArguments"]) &&
                fun["finiterange"] %in% finiterange &&
                (fun["maxdim"] < 0 ||
                 (fun["maxdim"] >= valid.in.dim[1] &&
                  fun["maxdim"] <= valid.in.dim[2]))  &&
                (fun["vdim"] < 0 ||
                 (fun["vdim"] >= vdim[1] && fun["vdim"] <= vdim[2]))  &&
                ls.RFmg[i] != ZF_INTERNALMIXED                
                )
  }
  ls.RFmg[(len-length(internals) + 1):len] <- internals
  
  return(unique(sort(ls.RFmg[idx])))
}


RFgetMethodNames <- function() {
  RFgetModelNames(type=TYPENAMES[c(GaussMethodType, BrMethodType) + 1])
}


RFformula <- function(f)
  return(parseModel(f))


GetProcessType <- function(model) {
  stopifnot(is.list(model))
  return(.Call("GetProcessType", MODEL_INTERN, model))
}


parameter.range <- function(model, param, dim=1){
  cat("sorry not programmed yet\n")
  return(NULL)
 

#parampositions < - function(model, param, trend=NULL, dim, print=1) {
#  stopifnot(!missing(dim))
#  model <- PrepareModel(model, param, trend=trend, nugget.remove=FALSE)
#  .Call("Get NA Positions", reg, model, as.integer(dim), as.integer(dim),
#        FALSE, FALSE, as.integer(print), PACKAGE="RandomFields")
#}
#  pm <- PrepareModel(model=model, param=param, nugget.remove=FALSE)        
#  storage.mode(dim) <- "integer"
#  ResGet <- .Call("SetAnd  ?? GetModelInfo",
#                  reg,
#                  pm, dim, Zeit, dim, FALSE, MaxNameCharacter=254,
#                  TRUE, TRUE,
#                  PACKAGE="RandomFields")
#  minmax <- ResGet$minmax[, 1:2]
#  dimnames(minmax) <-
#    list(attr(ResGet$minmax, "dimnames")[[1]], c("min", "max"))
#  return(minmax)
}



mergeWithGlobal <- function(dots) {
  if (exists(par.storage, envir=par.storage.env )) {
    current <- dots
    dots <- get(par.storage, envir=par.storage.env )
    names.current <- names(current)
    if (length(current) > 0) {
      for (i in 1:length(current))
        dots[[names.current[i]]] <- current[[i]]
    }
  }
  dots
}


RFpar <- function(...) {
  #l <- eval(substitute(list(...)))
  
  l <- list(...)
  if (length(l) == 1 && is.null(l[[1]]) && length(names(l)) ==0) {
    assign(par.storage, list(), envir=par.storage.env )
    return(NULL)
  }

  if (!exists(par.storage, envir=par.storage.env ))
    assign(par.storage, list(), envir=par.storage.env)

  par <- get(par.storage, envir=par.storage.env )
  if (length(l) == 0) return(par)

  n <- names(l)
  for (i in 1:length(l)) {
    par[[n[i]]] <- l[[i]]
  }
  assign(par.storage, par, envir=par.storage.env ) 
}
