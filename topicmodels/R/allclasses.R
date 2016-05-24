##**********************************************************
## control parameters

setClass("OPTcontrol",
         representation(iter.max = "integer",
                        tol      = "numeric"),
         prototype(iter.max = -1L,
                   tol      = sqrt(.Machine$double.eps)))

setClass("TopicModelcontrol",
         representation(seed          = "integer",
                        verbose       = "integer",
                        prefix        = "character",
                        save          = "integer",
                        nstart        = "integer",
                        best          = "logical",
                        keep          = "integer",
                        estimate.beta = "logical",
                        "VIRTUAL"),
         prototype(verbose = 0L,
                   save = 0L,
                   best = TRUE,
                   keep = 0L,
                   estimate.beta = TRUE))

init_TopicModelcontrol <- function(.Object, prefix, seed, nstart, ...) {
  if (missing(prefix)) prefix <- tempfile()
  if (missing(seed)) {
    if (missing(nstart)) {
      nstart <- 1L
      seed <- as.integer(Sys.time())
    } else {
      seed <- sample(seq_len(10^6), nstart)
    }
  } else if (missing(nstart)) nstart <- length(seed)
        
  list(.Object = .Object, prefix = prefix, seed = seed, nstart = nstart, ... = ...)
}

setMethod("initialize", "TopicModelcontrol", function(.Object, prefix, seed, nstart, ...) {
  args <- init_TopicModelcontrol(.Object, prefix, seed, nstart, ...)
  .Object <- do.call("callNextMethod", args)
  invisible(.Object)
})

setClass("VEMcontrol",
         representation(var    = "OPTcontrol",
                        em     = "OPTcontrol",
                        initialize = "character",
                        "VIRTUAL"),
         prototype(var        = new("OPTcontrol", iter.max = 500L, tol = 10^-6),
                   em         = new("OPTcontrol", iter.max = 1000L, tol = 10^-4),
                   initialize = "random"))

setMethod("initialize", "VEMcontrol", function(.Object, initialize = "random", ...) {
  initialize <- match.arg(initialize, c("random", "seeded", "model"))
  args <- init_TopicModelcontrol(.Object, ...)
  .Object <- do.call("callNextMethod",
                     c(args, initialize = initialize))
  invisible(.Object)
})

setClass("LDAcontrol",
         representation(alpha = "numeric",
                        "VIRTUAL"),
         contains = "TopicModelcontrol")

setClass("LDA_VEMcontrol",
         representation(estimate.alpha   = "logical"),
         contains = c("LDAcontrol", "VEMcontrol"), 
         prototype(estimate.alpha = TRUE))

setMethod("initialize", "LDA_VEMcontrol", function(.Object, prefix, initialize = "random", ...) {
  if (missing(prefix)) prefix <- tempfile()
  .Object <- callNextMethod(.Object = .Object, prefix = prefix, initialize = initialize, ...)
  invisible(.Object)
})

setClass("LDA_Gibbscontrol",
    representation(delta  = "numeric",
                   iter   = "integer",
                   thin   = "integer",
                   burnin = "integer",
                   initialize = "character"),
         contains = "LDAcontrol", 
         prototype(delta   = 0.1,
                   verbose = 0L,
                   iter    = 2000L,
                   burnin  = 0L,
                   nstart  = 1L,
                   best    = TRUE,
                   initialize = "random"))

setMethod("initialize", "LDA_Gibbscontrol", function(.Object, initialize = "random", seed = as.integer(NA), ...) {
  initialize <- match.arg(initialize, c("random", "beta", "z"))
  .Object <- callNextMethod(.Object = .Object, initialize = initialize, seed = seed, ...)
  if (length(.Object@thin) == 0) .Object@thin <- .Object@iter
  invisible(.Object)
})

setClass("CTM_VEMcontrol",
         representation(cg = "OPTcontrol"),
         contains = c("TopicModelcontrol", "VEMcontrol"),
         prototype(cg                   = new("OPTcontrol", iter.max = 500L, tol = 10^-5),
                   shrinkage.covariance = FALSE))

setMethod("initialize", "CTM_VEMcontrol", function(.Object, prefix, initialize = "random", ...) {
  if (missing(prefix)) prefix <- tempfile()
  .Object <- callNextMethod(.Object = .Object, prefix = prefix, initialize = initialize, ...)
  invisible(.Object)
})
##**********************************************************
## Topic Models Objects

setClass("TopicModel",
   representation(call            = "call",
                  Dim             = "integer",                  
                  control         = "TopicModelcontrol",
                  k               = "integer",
                  terms           = "ANY",
                  documents       = "ANY",
                  beta            = "matrix",
                  gamma           = "matrix",
                  wordassignments = "ANY",
                  loglikelihood   = "numeric",
                  iter            = "integer",
                  logLiks         = "numeric",
                  n               = "integer",
                  "VIRTUAL"))

setClass("VEM",
         contains = "TopicModel",
         representation("VIRTUAL"))

setClass("LDA",
         representation(alpha = "numeric",
                        "VIRTUAL"),
         contains = "TopicModel")

setClass("LDA_VEM",
         representation(),
         contains = c("LDA", "VEM"),
         prototype(control = new("LDA_VEMcontrol")))

setClass("Gibbs",
         contains = "TopicModel",
         representation("VIRTUAL"))

setClass("LDA_Gibbs",
         representation(seedwords = "ANY",
                        z = "integer"),
         contains = c("LDA", "Gibbs"),
         prototype(control = new("LDA_Gibbscontrol")))

setClass("Gibbs_list",
         representation(fitted = "list"))

setClass("CTM",
         representation(mu    = "numeric",
                        Sigma = "matrix",
                        "VIRTUAL"),
         contains = "TopicModel")

setClass("CTM_VEM",
         representation(nusquared = "matrix"),
         contains = c("CTM", "VEM"),
         prototype(control = new("CTM_VEMcontrol")))
         
setMethod("show", signature(object = "TopicModel"), function(object) {
  cat("A", class(object), "topic model with", object@k, "topics.\n")
})

