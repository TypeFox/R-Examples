## These functions handle SOM learning
initSOM <- function(dimension=c(5,5), topo=c("square"),
                    radius.type=c("gaussian", "letremy"),
                    dist.type=switch(match.arg(radius.type), 
                                     "letremy"="letremy", 
                                     "gaussian"="euclidean"),
                    type=c("numeric", "relational", "korresp"), 
                    mode=c("online"), affectation=c("standard", "heskes"),
                    maxit=500, nb.save=0, verbose=FALSE, proto0=NULL, 
                    init.proto=switch(type,
                                      "numeric"="random",
                                      "relational"="obs",
                                      "korresp"="random"),
                    scaling=switch(type,
                                   "numeric"="unitvar",
                                   "relational"="none",
                                   "korresp"="chi2"), eps0=1) {
  type <- match.arg(type)
  radius.type <- match.arg(radius.type)
  affectation <- match.arg(affectation)
  scaling <- match.arg(scaling,
                       c("unitvar", "none", "center", "chi2", "cosine"))
  dist.type <- match.arg(dist.type, c("letremy", "maximum", "euclidean",
                                      "manhattan", "canberra", "binary",
                                      "minkowski"))
  
  if (dist.type=="letremy" && radius.type=="gaussian") {
    dist.type <- "euclidean"
    warning("dist.type value replaced to 'euclidean' for Gaussian radius\n
            ('letremy' is not allowed)\n", 
            call.=TRUE, immediate.=TRUE)
  }
  
  init.proto <- match.arg(init.proto, c("random", "obs", "pca"))
  
  # check scaling compatibility
  if (type=="korresp" && scaling!="chi2") {
    scaling <- "chi2"
    warning("scaling value replaced: must be 'chi2' for 'korresp' type\n", 
            call.=TRUE, immediate.=TRUE)
  }
  if (type=="relational" && ! scaling %in% c("none", "cosine")) {
    scaling <- "none"
    warning("Wrong scaling for 'relational' SOM ; set to 'none'\n", call.=TRUE, 
            immediate.=TRUE)
  }
  if (type=="numeric" && scaling %in% c("chi2", "cosine"))
    stop(paste0("scaling='", scaling,
                "' is only implemented for 'korresp' type\n"), call.=TRUE)
  
  # check init.proto compatibility
  if (type=="korresp" && init.proto=="pca")
    stop("'init.proto' cannot be 'pca' for 'korresp' type\n", call.= TRUE)
  
  # check proto0
  if (!is.null(proto0)) {
    if (type=="relational") {
      if (sum(proto0<0)>0)
        stop("initial prototypes given by user do not match chosen type.
             Prototype values must be greater than 0.\n", call.=TRUE)
      if (sum(rowSums(proto0)!=1)>0)
        stop("initial prototypes given by user do not match chosen type.
             Prototype row sums for 'relational' must be equal to 1\n", 
             call.=TRUE)
    }
    if (type=="korresp") {
      if (min(proto0)<0 || max(proto0)>1)
        stop("initial prototypes given by user do not match chosen type.
               Prototypes for 'korresp' must have values between 0 and 1\n", 
             call.=TRUE)
    }
  }
  
  params <- list("the.grid"=initGrid(dimension,match.arg(topo),
                                     dist.type),
                 type=type, mode=match.arg(mode), affectation=affectation,
                 maxit=maxit, nb.save=nb.save, proto0=proto0,
                 init.proto=init.proto, scaling=scaling, 
                 radius.type=radius.type, verbose=verbose, eps0=eps0)
  
  class(params) <- "paramSOM"
  
  return(params)
}

print.paramSOM <- function(x,...){
  cat("\n  Parameters of the SOM\n\n")
  cat("    SOM mode                       : ", x$mode, "\n")
  cat("    SOM type                       : ", x$type, "\n")
  cat("    Affectation type               : ", x$affectation, "\n")
  cat("    Grid                           : ")
  print(x$the.grid)
  cat("    Number of iterations           : ", x$maxit, "\n")
  cat("    Number of intermediate backups : ", x$nb.save, "\n")
  if(identical(x$proto0,NULL)){
    cat("    Initializing prototypes method : ", x$init.proto, "\n")
  } else {
    cat("    Initial prototypes given by user\n")
  }
  cat("    Data pre-processing type       : ", x$scaling, "\n")
  cat("    Neighbourhood type             : ", x$radius.type, "\n")
  cat("\n")
}

summary.paramSOM <- function(object,...){
  cat("\nSummary\n\n")
  cat("  Class                            : ", class(object),"\n")
  print(object)
}