## __________________________________________________________
##
## SIMoNe
##
## (c) Statistique & Genome
## 12/02/2010
##______________________________________________________

simone <- function (X,
                    type = "steady-state",
                    clustering = FALSE,
                    tasks     = factor(rep(1,nrow(X))),
                    control   = setOptions() ) {
  
  ## ___________________________________________________________
  ##
  ## PRETREATMENTS: getting/checking the model parameters
  
  normalize      <- control$normalize
  verbose        <- control$verbose
  penalties      <- control$penalties
  penalty.min    <- control$penalty.min
  penalty.max    <- control$penalty.max
  n.penalties    <- control$n.penalties
  edges.max      <- control$edges.max
  edges.steady   <- control$edges.steady
  edges.coupling <- control$edges.coupling
  edges.sym.rule <- control$edges.sym.rule
  clusters.crit  <- control$clusters.crit
  clusters.qmin  <- control$clusters.qmin
  clusters.qmax  <- control$clusters.qmax
  clusters.meth  <- control$clusters.meth
    
  ## Centering (always) and normalizing if required (TRUE by default)
  if (nlevels(tasks) == 1) {
    X <- scale(X,TRUE,normalize)
  } else {
    X.s <- by(X,tasks,scale.default,TRUE,normalize)
    X <- NULL
    for (t in 1:nlevels(tasks)) {
      X <- rbind(X,X.s[[t]])
    }
  }
  
  ## Time-course data (time-series)
  if (type=="time-course") {
    if (nlevels(tasks) > 1) {
      method <- "multi.var1"      
    } else {
      method <- "var1.inference"
    }
  }
  ## steady-state data (replicats)
  else {
    if (nlevels(tasks) > 1) {
      method <- "multi.gaussian"
    } else {
      ## check for graphical.lasso or neighborhood.selection
      method <- edges.steady
    }
  }
  
  ## INITIALIZING THE NETWORK INFERENCE PARAMETERS
  ctrl.inf <- OptInference(method   = method,
                           tasks    = tasks,
                           coupling = edges.coupling,
                           sym.rule = edges.sym.rule,
                           solver   = ifelse(min(table(tasks)) <= ncol(X),
                             "active.set","shooting"))
  edges.sym.rule <- ctrl.inf$sym.rule
  
  ## GETTING THE RANGE OF THE PENALTY PARAMETER
  if (is.null(penalty.max)) {
    penalty.max <- get.rho.max(X,type,tasks)
  }
  if (is.null(penalty.min)) {
    if (nlevels(tasks) > 1) {
      penalty.min <- 1e-2
    } else {
      penalty.min <- 1e-5
    }
  }
  ## check that the user's choice is on the interval of possible values
  penalty.max <- min(penalty.max,get.rho.max(X,type,tasks))
  penalty.min <- max(penalty.min,.Machine$double.eps)
  if (is.null(penalties)) {
    penalties <- seq(penalty.max,penalty.min,length=n.penalties)
    keep.all <- FALSE
  } else {
    keep.all <- TRUE
  }
  penalties[penalties < penalty.min] <- penalty.min
  penalties[penalties > penalty.max] <- penalty.max
  
  ## save the updated parameters in the control object
  control <- setOptions(normalize      = normalize,
                        verbose        = verbose,
                        penalties      = penalties,
                        penalty.min    = penalty.min,
                        penalty.max    = penalty.max,                        
                        n.penalties    = n.penalties,
                        edges.max      = edges.max,
                        edges.sym.rule = edges.sym.rule,
                        edges.steady   = edges.steady,
                        edges.coupling = edges.coupling,
                        clusters.crit  = clusters.crit,
                        clusters.qmin  = clusters.qmin,
                        clusters.qmax  = clusters.qmax,
                        clusters.meth  = clusters.meth)
  
  ## ___________________________________________________________
  ##
  ## INFER NETWORK'S CLUSTERING
  if (clustering) {
    inferred.clustering <- InferStructure(X, type, control)
    clusters <- inferred.clustering$clusters
    weights  <- inferred.clustering$weights
  } else {
    clusters <- factor(rep("N",ncol(X)))
    weights  <- matrix(1,ncol(X),ncol(X))
  }

  penalties <- penalties / min(weights)
  
  ## ___________________________________________________________
  ##
  ## RUN THE MAIN LOOP
  if (nlevels(tasks) > 1) {
    m <- paste(method,"with",edges.coupling,"coupling and",edges.sym.rule,
               "symmetrization rule applied")
  } else {
    m <- paste(method,"with",edges.sym.rule,"symmetrization rule applied")
  }
  if (verbose) {
    cat("\nNetwork Inference:",m,"\n\n")
  }
  networks <- list()
  BIC      <- c()
  AIC      <- c()
  rhos     <- c()
  n.edges  <- c()
  last.edges <- rep(-1,nlevels(tasks))
  last.crit  <- -Inf
  if (verbose) {
    cat(format(c("|  penalty","|    edges","| criteria"), width=10,
               justify="right"),"\n\n")
  }
  for (rho in penalties){
    out <- InferEdges(X, rho*weights, ctrl.inf)
    ## if Rho is too small or if the previous guess has more non zero
    ## entries than possibly activated by the Lasso, then stop here.
    if (sum(is.na(out$n.edges))){
      break
    }
    ## if the number of edges decrease, it means the inference algo
    ## fails in correctly optimising the objective function: roughly,
    ## the penalty is too small and we better stop here
    if (sum(out$n.edges  - last.edges) < 0) {
      break
    }    
    ctrl.inf$initial.guess <- out$Beta
    ## if the current penalty does not change the network's clustering,
    ## go to the following one and keep the smallest corresponding penalty
    if (!keep.all) { ## unless the user provided its own penalties
      if (all(out$n.edges == last.edges)) {
        rhos[length(rhos)] <- rho
        next
      }
    }
    BIC  <- c(BIC,out$BIC)
    AIC  <- c(AIC,out$AIC)
    rhos <- c(rhos,rho)
    last.edges <- out$n.edges
    n.edges    <- rbind(n.edges,out$n.edges)
    last.crit  <- out$loglik.l1
    networks[[length(networks)+1]] <- out$Theta
    if (verbose) {
      cat(format(list(rho,paste(last.edges,collapse=","),last.crit),
               width=10, digits=4, justify="right"),"\n")
    }
    if (max(last.edges) > edges.max) {
      break
    }
  }

  ## ___________________________________________________________
  ##
  ## RETURN THE RESULT
  return(structure(list(networks  = networks,
                        penalties = as.vector(rhos),
                        n.edges   = as.matrix(n.edges),
                        BIC       = as.vector(BIC),
                        AIC       = as.vector(AIC),
                        clusters  = as.factor(clusters),
                        weights   = as.matrix(weights),
                        control   = control), class="simone"))
}

get.rho.max <- function(X,type,tasks) {
  
  if (type =="steady-state") {
    r.max <- max(sapply(sapply(by(X,tasks,var),abs),max))
  } else {
    r.max <- 0
    for (t in levels(tasks)) {
      x <- X[tasks == t,]
      r.max <- max(r.max,max(abs(t(x[-nrow(x),]) %*% x[-1,]/(nrow(x)-1))))
    }
  } 
  return(r.max)
}

setOptions <- function(normalize      = TRUE,
                       verbose        = TRUE,
                       penalties      = NULL,
                       penalty.min    = NULL,
                       penalty.max    = NULL,
                       n.penalties    = 100, 
                       edges.max      = Inf,
                       edges.sym.rule = NULL,
                       edges.steady   = "neighborhood.selection",
                       edges.coupling = "coopLasso",
                       clusters.crit  = "BIC",
                       clusters.meth  = "bayesian",
                       clusters.qmin  = 2,
                       clusters.qmax  = 4) {
  
  return(list(normalize      = normalize,
              verbose        = verbose,
              penalties      = penalties,
              penalty.min    = penalty.min, 
              penalty.max    = penalty.max, 
              n.penalties    = n.penalties,
              edges.max      = edges.max,
              edges.sym.rule = edges.sym.rule,
              edges.steady   = edges.steady,
              edges.coupling = edges.coupling,
              clusters.crit  = clusters.crit,
              clusters.meth  = clusters.meth,
              clusters.qmax  = clusters.qmin,
              clusters.qmin  = clusters.qmax))
}
