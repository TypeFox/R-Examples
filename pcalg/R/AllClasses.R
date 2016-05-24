##################################################
### Part 1 : S4 classes used by pc and r/fci
##################################################

## $Id: AllClasses.R 342 2015-07-22 14:26:29Z mmaechler $

setClass("gAlgo",
         representation(call = "call",
                        n = "integer",
                        max.ord = "integer",
                        n.edgetests= "numeric",
                        sepset= "list",
                        pMax= "matrix"), "VIRTUAL")


setClass("fciAlgo",
         representation(amat = "matrix", allPdsep = "list",
                        n.edgetestsPDSEP = "numeric", max.ordPDSEP = "integer"),
         contains = "gAlgo")

setClass("pcAlgo",
         representation(graph = "graph", zMin = "matrix"), ## zMin for compatibility
         contains = "gAlgo")

## Methods

##' Extract graph part of an R object:
setGeneric("getGraph", function(x) as(x,"graph"))
setMethod("getGraph", "matrix", function(x) as(x, "graphAM"))
if(FALSE) {## if we would importFrom("Matrix", ....) in NAMESPACE
    setMethod("getGraph", "sparseMatrix", function(x) as(x, "graphNEL"))
    setMethod("getGraph", "Matrix", function(x) as(x, "graphAM"))
}
setMethod("getGraph", "pcAlgo", function(x) x@graph)
setMethod("getGraph", "fciAlgo", function(x) as(x@amat, "graphAM"))

setMethod("summary", "pcAlgo",
          function(object) {
 	    cat("\nObject of class 'pcAlgo', from Call: \n",
                deparse(object@call),
 		"\n\nNmb. edgetests during skeleton estimation:\n")
            cat("===========================================\n")
            cat("Max. order of algorithm: ",object@max.ord,
                "\nNumber of edgetests from m = 0 up to m =",object@max.ord,
                ": ",object@n.edgetests)
            nbrs <- vapply(object@graph@edgeL, function(x) length(x$edges), 1L)
            cat("\n\nGraphical properties of skeleton:\n")
            cat("=================================\n")
            cat("Max. number of neighbours: ", max(nbrs),
                "at node(s)", which(nbrs==max(nbrs)),
                "\nAvg. number of neighbours: ",mean(nbrs),"\n")
          })

setMethod("summary", "fciAlgo",
          function(object) {
 	    cat("Object of class 'fciAlgo'\n\n")
            cat("Call: \n=====\n", deparse(object@call))
            cat("\n\nNmb. edgetests during skeleton estimation:\n==========================================")
            cat("\nMax. order of algorithm: ",object@max.ord,
                "\nNumber of edgetests from m = 0 up to m =",object@max.ord,
                ": ",object@n.edgetests)
            cat("\n\nAdd. nmb. edgetests when using PDSEP:\n=====================================")
            cat("\nMax. order of algorithm: ",object@max.ordPDSEP,
                "\nNumber of edgetests from m = 0 up to m =",object@max.ordPDSEP,
                ": ",object@n.edgetestsPDSEP)

            myLength <- function(x) if(is.null(x)) NA_integer_ else length(x)
            cat("\n\nSize distribution of SEPSET:")
            myTab <- table(sapply(object@sepset,
                                  function(x) vapply(x, myLength, 1L)),
                           useNA = "always")
            print(myTab)

            cat("\nSize distribution of PDSEP:")
            print(table(vapply(object@allPdsep, length, 1L)))
          })


setMethod("show", "pcAlgo",
	  function(object) {
	    cat("Object of class 'pcAlgo', from Call: \n", deparse(object@call),"\n")
            amat <- as(object@graph, "matrix")
            amat2 <- amat + 2*t(amat)
            ude <- sum(amat2 == 3)/2
            de <- sum(amat2 == 1)
            cat("Number of undirected edges: ", ude, "\n")
            cat("Number of directed edges:   ", de, "\n")
            cat("Total number of edges:      ", de + ude, "\n")
	    invisible(object)
	  })


print.fciAlgo <- function(x, zero.print = ".", ...) {
    cat("Object of class 'fciAlgo', from Call:", deparse(x@call),
        "\nAdjacency Matrix G:",
        "G[i,j] = 1/2/3 if edge mark of edge i-j at j is circle/head/tail.",
        "", sep="\n")
    print.table(x@amat, zero.print=zero.print, ...)
    invisible(x)
}

setMethod("show", "fciAlgo", function(object) print.fciAlgo(object))

## -> ../man/pcAlgo-class.Rd
setMethod("plot", signature(x = "pcAlgo"),
          function(x, y, main = NULL, zvalue.lwd = FALSE,
                   lwd.max = 7, labels = NULL, ...)
	{
          check.Rgraphviz()

          if(is.null(main))
              main <- deparse(x@call)
          attrs <- nodeAttrs <- list()
          p <- numNodes(G <- x@graph)
          if (!is.null(labels)) {
              attrs$node <- list(shape = "ellipse", fixedsize = FALSE)
              names(labels) <- nodes(G)
              nodeAttrs$label <- labels
          }

          if (zvalue.lwd && numEdges(G) != 0) {
	      lwd.mat <-
		  if(is.matrix(Z <- x@zMin) && all(dim(Z) == p)) Z
		  else ## from newer pc(): 'zMin' is deprecated there, but pMax corresponds:
		      qnorm(x@pMax/2, lower.tail=FALSE)
	      lwd.mat <- lwd.max * lwd.mat/max(lwd.mat)
	      z <- Rgraphviz::agopen(G, name = "lwdGraph",
				     nodeAttrs = nodeAttrs, attrs = attrs)
	      for (i in seq_along(z@AgEdge)) {
		  z@AgEdge[[i]]@lwd <- lwd.mat[as.integer(z@AgEdge[[i]]@head),
					       as.integer(z@AgEdge[[i]]@tail)]
	      }
              Rgraphviz::plot(z, main = main, ...)
          } else {
              Rgraphviz::plot(G, nodeAttrs = nodeAttrs, main = main,
                              attrs = attrs, ...)
          }
      })

setMethod("plot", signature(x = "fciAlgo"),
          function(x, y, main = NULL, ...)
      {
          check.Rgraphviz()

          if(is.null(main))
	      main <- deparse(x@call)
	  else ## see also below
	      warning("main title cannot *not* be set yet [Rgraphviz::plot() deficiency]")
          amat <- x@amat
          g <- as(amat,"graphNEL")
          nn <- nodes(g)
          p <- numNodes(g)
          n.edges <- numEdges(g)
          ah.list <- at.list <- rep("none",n.edges)
          counter <- 0
          list.names <- NULL
          amat[amat==1] <- "odot"
          amat[amat==2] <- "normal"
          amat[amat==3] <- "none"
          for (i in seq_len(p-1)) {
              for (j in (i+1):p) {
                  x <- nn[i]
                  y <- nn[j]
                  if (amat[x,y]!=0) {
                      counter <- counter + 1
                      ah.list[[counter]] <- amat[x,y]
                      at.list[[counter]] <- amat[y,x]
                      list.names <- c(list.names,paste(x,"~",y,sep=""))
                  }
              }
          }
          names(ah.list) <- names(at.list) <- list.names
	  edgeRenderInfo(g) <- list(arrowhead= ah.list,
				    arrowtail= at.list)
          ## XXX Sep/Oct 2010  --- still current -- FIXME ??
          ## XXX undid change by MM, since edge marks didn't work anymore
          ## XXX "known bug in Rgraphviz, but not something they may fix soon"
	  ## Rgraphviz::plot(g, main = main, ...)
          Rgraphviz::renderGraph(Rgraphviz::layoutGraph(g))
      })

#######################################################
### Part 2 : Reference classes and Methods used by GIES
#######################################################

##' Virtual base class for all parametric causal models.
##' The meaning of the "params" depends on the model used.
setRefClass("ParDAG",
    fields = list(
        .nodes = "vector",
        .in.edges = "list",
        .params = "list"),

    validity = function(object) {
      if (anyDuplicated(object$.nodes))
        return("The node names must be unique")
      if (any(names(object$.in.edges) != object$.nodes))
        return("The elements of 'in.edges' must be named after the nodes.")
      if (!all(sapply(object$.in.edges, is.numeric)))
        return("The vectors in 'in.edges' must contain numbers.")

      edgeRange <- range(unlist(object$.in.edges))
      if (object$edge.count() > 0 &&
          (edgeRange[1] < 1 || edgeRange[2] > object$node.count()))
        return("Invalid range of edge sources.")

      return(TRUE)
    },

    methods = list(
        #' Constructor
        initialize = function(nodes, in.edges = NULL, params = list()) {
          .nodes <<- nodes

          if (is.null(in.edges))
            .in.edges <<- replicate(length(nodes), integer(0), simplify = FALSE)
          else
            .in.edges <<- lapply(1:length(in.edges), function(i) as.integer(in.edges[[i]]))
          names(.in.edges) <<- nodes

          .params <<- params
        },

        #' Yields the number of nodes
        node.count = function() {
          length(.nodes)
        },

        #' Yields the total number of edges in the graph
        edge.count = function() {
          sum(sapply(.in.edges, length))
        },

        #' Simulates (draws a sample of) interventional (or observational) data
        simulate = function(n, target = integer(0), int.level = numeric(0)) {
          stop("simulate() is not implemented in this class.")
        },

        #' Fits parameters by MLE using a scoring object
        mle.fit = function(score) {
          .params <<- score$global.mle(.self)
        }
        ),

    "VIRTUAL")

#' Coercion to a graphNEL instance
setAs("ParDAG", "graphNEL",
      function(from) {
          reverseEdgeDirections(
              new("graphNEL",
                  nodes = from$.nodes,
                  edgeL = from$.in.edges,
                  edgemode = "directed"))
      })

#' Coercion to a (logical) matrix
setAs("ParDAG", "matrix",
      function(from) {
          i.p <- seq_len(p <- from$node.count())
          in.edge <- from$.in.edges
	  vapply(i.p, function(i) i.p %in% in.edge[[i]], logical(p))
      })

#' Plot method (needs Rgraphviz to work!!)
setMethod("plot", "ParDAG",
    function(x, y, ...) {
      if (!validObject(x))
        stop("The parametric DAG model to be plotted is not valid")

      if (missing(y))
        y <- "dot"
      invisible(plot(as(x, "graphNEL"), y, ...))
    })


#' Virtual base class for all scoring classes
setRefClass("Score",
    fields = list(
        decomp = "logical",
        c.fcn = "character",
        pp.dat = "list",
        .pardag.class = "character"),

    validity = function(object) {
      ## Check if targets are valid (i.e., unique)
      targets.tmp <- object$pp.dat$targets
      for (i in seq(along = targets.tmp)) {
        targets.tmp[[i]] <- sort(unique(targets.tmp[[i]]))
        if (length(targets.tmp[[i]]) != length(object$pp.dat$targets[[i]]))
          return("Target variables must not be listed multiple times.")
      }
      if (length(unique(targets.tmp)) != length(targets.tmp))
        return("Targets must not be listed multiple times.")

      ## Check whether data is available from all intervention targets
      if (unique(object$pp.dat$target.index) != 1:length(object$pp.dat$targets))
        return("Data from all intervention targets must be available")

      ## Check if dimensions of target.index and data conincide
      if (length(object$pp.dat$target.index) != nrow(object$pp.dat$data))
        return("Length of target index vector does not coincide with sample size.")

      return(TRUE)
    },

    methods = list(
        #' Constructor
        #'
        #' Note: all arguments must have a default value for inheritance,
        #' see ?setRefClass; apart from that, the default values are meaningless
        initialize = function(data = matrix(1, 1, 1),
            targets = list(integer(0)),
            target.index = rep(as.integer(1), nrow(data)),
            ...) {
          ## Order by ascending target indices (necessary for certain scoring objects)
          if (is.unsorted(target.index))
            perm <- order(target.index)
          else
            perm <- 1:length(target.index)

          pp.dat$targets <<- lapply(targets, sort)
          pp.dat$target.index <<- target.index[perm]
          pp.dat$data <<- data[perm, ]
          pp.dat$vertex.count <<- ncol(data)
          pp.dat$total.data.count <<- as.integer(nrow(data))

          ## Declare scores as not decomposable "by default"
          decomp <<- FALSE

          ## No C++ scoring object by default
          c.fcn <<- "none"

          ## R function objects
          pp.dat$local.score <<- function(vertex, parents) local.score(vertex, parents)
          pp.dat$global.score <<- function(edges) global.score(vertex, parents)
          pp.dat$local.mle <<- function(vertex, parents) local.mle(vertex, parents)
          pp.dat$global.mle <<- function(edges) global.mle(vertex, parents)

          callSuper(...)
        },

        #' Checks whether a vertex is valid
        validate.vertex = function(vertex) {
          stopifnot(is.whole(vertex))
          stopifnot(abs(vertex - round(vertex)) < sqrt(.Machine$double.eps))
          stopifnot(1 <= vertex && vertex <= pp.dat$vertex.count)
        },

        #' Checks whether a vector is a valid list of parents
        validate.parents = function(parents) {
          stopifnot(all(parents %in% 1:pp.dat$vertex.count))
          stopifnot(anyDuplicated(parents) == 0)
        },

        #' Getter and setter function for the targets
        getTargets = function() {
          pp.dat$targets
        },

        setTargets = function(targets) {
          pp.dat$targets <<- lapply(targets, sort)
        },

        #' Creates a list of options for the C++ functions for the internal
        #' calculation of scores and MLEs
        c.fcn.options = function(DEBUG.LEVEL = 0) {
          list(DEBUG.LEVEL = DEBUG.LEVEL)
        },

        #' Calculates the local score of a vertex and its parents
        local.score = function(vertex, parents, ...) {
          stop("local.score is not implemented in this class.")
        },

        #' Calculates the global score of a DAG which is only specified
        #' by its list of in-edges
        global.score.int = function(edges, ...) {
          ## Calculate score in R
          if (c.fcn == "none")
            sum(sapply(1:pp.dat$vertex.count,
                    function(i) local.score(i, edges[[i]], ...)))
          ## Calculate score with the C++ library
          else
            .Call("globalScore", c.fcn, pp.dat, edges, c.fcn.options(...), PACKAGE = "pcalg")
        },

        #' Calculates the global score of a DAG
        global.score = function(dag, ...) {
          global.score.int(dag$.in.edges, ...)
        },

        #' Calculates the local MLE for a vertex and its parents
        local.mle = function(vertex, parents, ...) {
          stop("local.mle is not implemented in this class.")
        },

        #' Calculates the global MLE
        global.mle = function(dag, ...) {
          ## Calculate score in R
          if (c.fcn == "none") {
              in.edge <- dag$.in.edges
              lapply(1:pp.dat$vertex.count,
                     function(i) local.mle(i, in.edge[[i]], ...))
          }
          ## Calculate score with the C++ library
          else
            .Call("globalMLE", c.fcn, pp.dat, dag$.in.edges, c.fcn.options(...),
                  PACKAGE = "pcalg")
        }
        ),

    "VIRTUAL")

#' l0-penalized log-likelihood for Gaussian models, with freely
#' choosable penalty lambda.
#' Special case: BIC where \lambda = 1/2 \log n (default value for lambda)
setRefClass("GaussL0penIntScore",
    contains = "Score",

    fields = list(
        .format = "character"),

    validity = function(object) {
      p <- ncol(object$pp.dat$data)
      if (!is.null(object$pp.dat$scatter)) {
          ## Data storage with precalculated scatter matrices
          if (unique(object$pp.dat$scatter.index) != 1:length(object$pp.dat$scatter))
            return("The index list of distinct scatter matrices has an invalid range.")
          if (any(sapply(object$pp.dat$scatter, function(mat) dim(mat) != c(p, p))))
            return("The scatter matrices have invalid dimensions.")
        }

      return(TRUE)
    },

    methods = list(
        #' Constructor
        initialize = function(data = matrix(1, 1, 1),
            targets = list(integer(0)),
            target.index = rep(as.integer(1), nrow(data)),
            lambda = 0.5*log(nrow(data)),
            intercept = FALSE,
            format = c("raw", "scatter"),
            use.cpp = TRUE,
            ...) {
          ## Store supplied data in sorted form
          callSuper(data = data, targets = targets, target.index = target.index, ...)

          ## Number of variables
          p <- ncol(data)

          ## l0-penalty is decomposable
          decomp <<- TRUE

          ## Underlying causal model class: Gaussian
          .pardag.class <<- "GaussParDAG"

          ## Store different settings
          pp.dat$lambda <<- lambda

          ## Store data format. Currently supporting scatter matrices
          ## and raw data only (recommended for high-dimensional data)
          .format <<- match.arg(format)
          ## If format not specified by user, choose it based on dimensions
          ## TODO: check if this choice is reasonable...
          if (length(format) > 1)
            .format <<- ifelse(p >= nrow(data) || p >= 500, "raw", "scatter")
          ## TODO change following line as soon as "raw" format is implemented and tested
          .format <<- "scatter"

          ## Use C++ functions if requested
          if (use.cpp)
            c.fcn <<- ifelse(.format == "scatter", "gauss.l0pen.scatter", "gauss.l0pen.raw")

          ## Add column of ones to data matrix to calculate scatter matrices;
          ## this allows the computation of an intercept if requested
          pp.dat$intercept <<- intercept
          data <- cbind(pp.dat$data, 1)# take matrix that is already pre-processed,
          # having reordered rows!

          ## Create scatter matrices for different targets
          ti.lb <- c(sapply(1:length(pp.dat$targets), function(i) match(i, pp.dat$target.index)),
              length(pp.dat$target.index) + 1)
          scatter.mat <- lapply(1:length(pp.dat$targets),
              function(i) crossprod(data[ti.lb[i]:(ti.lb[i + 1] - 1), , drop = FALSE]))

          ## Find all interventions in which the different variables
          ## are _not_ intervened
          non.ivent <- matrix(FALSE, ncol = p, nrow = length(pp.dat$targets))
          pp.dat$scatter.index <<- integer(p)
          pp.dat$data.count <<- integer(p)
          max.si <- 0
          for (i in 1:p) {
            ## Generate indices of (distinct) scatter matrices
            non.ivent[ , i] <- sapply(seq_along(pp.dat$targets),
                                      function(j) i %nin% pp.dat$targets[[j]])
            pp.dat$scatter.index[i] <<- max.si + 1
            j <- 1
            while (j < i) {
              if (all(non.ivent[, i] == non.ivent[, j])) {
                pp.dat$scatter.index[i] <<- pp.dat$scatter.index[j]
                j <- i
              }
              j <- j + 1
            }
            if (pp.dat$scatter.index[i] == max.si + 1)
              max.si <- max.si + 1

            ## Count data samples from "non-interventions" at i
            pp.dat$data.count[i] <<- sum(ti.lb[which(non.ivent[, i]) + 1] - ti.lb[which(non.ivent[, i])])
          }

          ## Calculate the distinct scatter matrices for the
          ## "non-interventions"
          pp.dat$scatter <<- lapply(1:max.si,
             function(i) Reduce("+", scatter.mat[non.ivent[, match(i, pp.dat$scatter.index)]]))
        },

        #' Calculates the local score of a vertex and its parents
        local.score = function(vertex, parents, ...) {
          ## Check validity of arguments
          validate.vertex(vertex)
          validate.parents(parents)

          ## Calculate score in R
          if (c.fcn == "none") {
            ## If an intercept is allowed, add a fake parent node
            parents <- sort(parents)
            if (pp.dat$intercept)
              parents <- c(pp.dat$vertex.count + 1, parents)

            sigma2 <- pp.dat$scatter[[pp.dat$scatter.index[vertex]]][vertex, vertex]
            if (length(parents) != 0) {
              b <- pp.dat$scatter[[pp.dat$scatter.index[vertex]]][vertex, parents]
              sigma2 <- sigma2 - as.numeric(b %*% solve(pp.dat$scatter[[pp.dat$scatter.index[vertex]]][parents, parents], b))
            }

            return(-0.5*pp.dat$data.count[vertex]*(1 + log(sigma2/pp.dat$data.count[vertex])) - pp.dat$lambda*(1 + length(parents)))
          }
          ## Calculate score with the C++ library
          else
            return(.Call("localScore", c.fcn, pp.dat, vertex, parents, c.fcn.options(...), PACKAGE = "pcalg"))
        },

        #' Calculates the local MLE for a vertex and its parents
        local.mle = function(vertex, parents, ...) {
          ## Check validity of arguments
          validate.vertex(vertex)
          validate.parents(parents)

          ## Calculate score in R
          if (c.fcn == "none") {
            ## If an intercept is allowed, add a fake parent node
            parents <- sort(parents)
            if (pp.dat$intercept)
              parents <- c(pp.dat$vertex.count + 1, parents)

            sigma2 <- pp.dat$scatter[[pp.dat$scatter.index[vertex]]][vertex, vertex]
            if (length(parents) != 0) {
              beta <- solve(pp.dat$scatter[[pp.dat$scatter.index[vertex]]][parents, parents],
                  pp.dat$scatter[[pp.dat$scatter.index[vertex]]][vertex, parents])
              sigma2 <- sigma2 - pp.dat$scatter[[pp.dat$scatter.index[vertex]]][vertex, parents] %*% beta
            }
            else
              beta <- numeric(0)

            if (pp.dat$intercept)
              return(c(sigma2/pp.dat$data.count[vertex], beta))
            else
              return(c(sigma2/pp.dat$data.count[vertex], 0, beta))
          }
          ## Calculate score with the C++ library
          else
            return(.Call("localMLE", c.fcn, pp.dat, vertex, parents, c.fcn.options(...), PACKAGE = "pcalg"))
        }
        )
    )

##' Observational score as special case
setRefClass("GaussL0penObsScore",
    contains = "GaussL0penIntScore",

    methods = list(
        #' Constructor
        initialize = function(data = matrix(1, 1, 1),
            lambda = 0.5*log(nrow(data)),
            intercept = FALSE,
            use.cpp = TRUE,
            ...) {
          callSuper(data = data,
              targets = list(integer(0)),
              target.index = rep(as.integer(1), nrow(data)),
              lambda = lambda,
              intercept = intercept,
              use.cpp = use.cpp,
              ...)
          }
        )
    )

#' Interventional essential graph
setRefClass("EssGraph",
    fields = list(
        .nodes = "vector",
        .in.edges = "list",
        .targets = "list",
        .score = "Score"
    ),

    validity = function(object) {
      ## Check nodes
      if (any(names(object$.in.edges) != object$.nodes)) {
        return("The elements of 'in.edges' must be named after the nodes.")
      }

      ## Check in-edges
      if (!all(sapply(object$.in.edges, is.numeric))) {
        return("The vectors in 'in.edges' must contain numbers.")
      }
      if (!all(unique(unlist(object$.in.edges)) %in% 1:object$node.count())) {
        return(sprintf("Invalid edge source(s): edge sources must be in the range 1:%d.",
          object$node.count()))
      }

      ## Check targets
      if (anyDuplicated(object$.targets)) {
        return("Targets are not unique.")
      }
      if (!all(unique(unlist(object$.targets)) %in% 1:object$node.count())) {
        return(sprintf("Invalid target(s): targets must be in the range 1:%d.",
          object$node.count()))
      }

      ## Check score
      if (!is.null(score <- object$getScore())) {
        targets <- object$getTargets()
        if (length(score$getTargets()) != length(targets) ||
            !all.equal(duplicated(c(targets, score$getTargets())),
                    rep(c(FALSE, TRUE), each = length(targets)))) {
          return("Targets do not coincide with that of the scoring object.")
        }
      }

      return(TRUE)
    },

    methods = list(
        #' Constructor
        initialize = function(nodes,
            in.edges = replicate(length(nodes), integer(0)),
            targets = list(integer(0)),
            score = NULL) {
          ## Store nodes names
          if (missing(nodes)) {
            stop("Argument 'nodes' must be specified.")
          }
          .nodes <<- as.character(nodes)

          ## Store in-edges
          # TODO: improve error checking; possibly put it into separate function
          stopifnot(is.list(in.edges) && length(in.edges) == length(nodes))
          .in.edges <<- in.edges
          names(.in.edges) <<- .nodes

          ## Store targets
          setTargets(targets)

          ## Store score
          setScore(score)
        },

        #' Yields the number of nodes
        node.count = function() {
          length(.nodes)
        },

        #' Yields the total number of edges in the graph
        edge.count = function() {
          sum(vapply(.in.edges, length, 1L))
        },

        #' Getter and setter functions for score object
        getScore = function() {
          .score
        },

        setScore = function(score) {
          if (!is.null(score)) {
            .score <<- score
          }
        },

        #' Getter and setter functions for targets list
        getTargets = function() {
          .targets
        },

        setTargets = function(targets) {
          .targets <<- lapply(targets, sort)
        },

        #' Creates a list of options for the C++ function "causalInference";
        #' internal function
        causal.inf.options = function(caching = TRUE,
            turning = TRUE,
            maxDegree = integer(0),
            maxSteps = 0,
            childrenOnly = integer(0),
            fixedGaps = NULL,
            verbose = 0) {
          list(caching = caching,
              turning = turning,
              maxDegree = maxDegree,
              maxSteps = maxSteps,
              childrenOnly = childrenOnly,
              fixedGaps = fixedGaps,
              DEBUG.LEVEL = as.integer(verbose))
        },

        #' Performs one greedy step
        greedy.step = function(direction = c("forward", "backward", "turning"), verbose = FALSE, ...) {
          stopifnot(!is.null(score <- getScore()))

          ## Cast direction
          direction <- match.arg(direction)
          alg.name <- switch(direction,
              forward = "GIES-F",
              backward = "GIES-B",
              turning = "GIES-T")

          new.graph <- .Call("causalInference",
              .in.edges,
              score$pp.dat,
              alg.name,
              score$c.fcn,
              causal.inf.options(caching = FALSE, maxSteps = 1, verbose = verbose, ...),
              PACKAGE = "pcalg")
          if (identical(new.graph, "interrupt"))
            return(FALSE)

          if (new.graph$steps > 0) {
            .in.edges <<- new.graph$in.edges
            names(.in.edges) <<- .nodes
          }

          return(new.graph$steps == 1)
        },

        greedy.search = function(direction = c("forward", "backward", "turning")) {
          stopifnot(!is.null(score <- getScore()))

          ## Cast direction
          direction <- match.arg(direction)
          alg.name <- switch(direction,
              forward = "GIES-F",
              backward = "GIES-B",
              turning = "GIES-T")

          new.graph <- .Call("causalInference",
              .in.edges,
              score$pp.dat,
              alg.name,
              score$c.fcn,
              causal.inf.options(caching = FALSE),
              PACKAGE = "pcalg")
          if (identical(new.graph, "interrupt"))
            return(FALSE)

          if (new.graph$steps > 0) {
            .in.edges <<- new.graph$in.edges
            names(.in.edges) <<- .nodes
          }

          return(new.graph$steps)
        },

        #' Performs a causal inference from an arbitrary start DAG
        #' with a specified algorithm
        caus.inf = function(algorithm, ...) {
          stopifnot(!is.null(score <- getScore()))
          stopifnot(algorithm %in% c("GIES", "GIES-F", "GIES-B", "GIES-T", "GIES-STEP", "GDS", "SiMy"))

          new.graph <- .Call("causalInference",
              .in.edges,
              score$pp.dat,
              algorithm,
              score$c.fcn,
              causal.inf.options(...),
              PACKAGE = "pcalg")

          if (identical(new.graph, "interrupt"))
            return(FALSE)
          else {
            .in.edges <<- new.graph$in.edges
            names(.in.edges) <<- .nodes
            return(TRUE)
          }
        },

        #' Performs GIES from an arbitrary start DAG
        gies = function(...) caus.inf("GIES", ...),

        #' Performs GDS from an arbitrary start DAG
        gds = function(...) caus.inf("GDS", ...),

        #' DP search of Silander and Myllymäki (ignores the start DAG!)
        silander = function(...) caus.inf("DP", ...),

        #' Calculates the parameters of a DAG via MLE (wrapper function only)
        mle.fit = function(dag) {
          stopifnot(!is.null(score <- getScore()))
          dag$.params <- score$global.mle(dag)
          return(dag)
        },

        #' Yields a representative (estimating parameters via MLE)
        repr = function() {
          stopifnot(!is.null(score <- getScore()))
          in.edges <- .Call("representative", .in.edges, PACKAGE = "pcalg")
          result <- new(score$.pardag.class, nodes = .nodes, in.edges = in.edges)
          result$.params <- score$global.mle(result)

          return(result)
        },

        #' Calculates an optimal intervention target
        #'
        #' @param   max.size    maximum target size; allowed values: 1, p (= # nodes)
        ## TODO document that function... or better: provide a documented wrapper function
        opt.target = function(max.size) {
          .Call("optimalTarget", .in.edges, max.size, PACKAGE = "pcalg")
        }
        ))

##' Coercion to a graphNEL instance
.ess2graph <- function(from)
    reverseEdgeDirections(new("graphNEL",
                              nodes = from$.nodes,
                              edgeL = from$.in.edges,
                              edgemode = "directed"))

setAs("EssGraph", "graphNEL", .ess2graph)
setAs("EssGraph", "graph", .ess2graph)

## NOTE: Coercion to SparseMatrix is more efficient via
## ----  via "graphNEL" for larger p :
##' Coercion to a (logical) matrix
setAs("EssGraph", "matrix",
      function(from) {
          ip <- seq_len(p <- from$node.count())
          in.edge <- from$.in.edges
          vapply(ip, function(i) ip %in% in.edge[[i]], logical(p))
      })

#' Plot method (needs Rgraphviz to work!!)
## TODO maybe adapt method to make sure that undirected edges are not plotted as
## bidirected
setMethod("plot", "EssGraph",
    function(x, y, ...) {
      if (!validObject(x))
        stop("Invalid parametric DAG model (\"EssGraph\")")
      if (missing(y))
        y <- "dot"
      invisible(plot(.ess2graph(x), y, ...))
    })

#' Gaussian causal model
setRefClass("GaussParDAG",
    contains = "ParDAG",

    validity = function(object) {
      if (any(names(object$.params) != object$.nodes))
        return("The elements of 'params' must be named after the nodes.")
      if (!all(sapply(1:object$node.count(),
          function(i) length(object$.params[[i]]) == length(object$.in.edges[[i]]) + 2)))
        return("The number of parameters does not match the number of in-edges.")

      return(TRUE)
    },

    methods = list(
        #' Yields the intercept
        intercept = function() {
          sapply(.params, function(par.vec) par.vec[2])
        },

        #' Sets the intercept
        set.intercept = function(value) {
          for (i in 1:node.count())
            .params[[i]][2] <<- value[i]
        },

        #' Yields the error variances
        err.var = function() {
          sapply(.params, function(par.vec) par.vec[1])
        },

        #' Sets the error variances
        set.err.var = function(value) {
          for (i in 1:node.count())
            .params[[i]][1] <<- value[i]
        },

        #' Yields the weight matrix w.r.t. an intervention target
        #'
        #' TODO add a method for sparse matrices...
        weight.mat = function(target = integer(0)) {
          ## Fill in weights
          p <- node.count()
          target <- as.integer(sort(target))
          result <- matrix(0, p, p)
          for (i in 1:p)
            if (as.integer(i) %nin% target)
              result[.in.edges[[i]], i] <- .params[[i]][-(1:2)]

          ## Set row and column names
          rownames(result) <- .nodes
          colnames(result) <- .nodes

          return(result)
        },

        #' Yields an observational or interventional covariance matrix
        #'
        #' @param   target    intervention target
        #' @param   ivent.var variances of the intervention variables
        #' @return  (observational or interventional) covariance matrix
        cov.mat = function(target = integer(0), ivent.var = numeric(0)) {
          A <- -weight.mat()
          A[, target] <- 0
          diag(A) <- 1
          A <- solve(A)

          all.var <- err.var()
          all.var[target] <- ivent.var

          return(t(A) %*% diag(all.var) %*% A)
        },

        #' Simulates (draws a sample of) interventional (or observational)
        #' data
        #'
        #' @param   n
        #' @param   target
        #' @param   int.level   intervention level: values of the intervened
        #'                      variables. Either a vector of the same length
        #'                      as "target", or a matrix with dimensions
        #'                      n x length(target)
        #' @return  a vector with the simulated values if n = 1, or a matrix
        #'          with rows corresponding to different samples if n > 1
        simulate = function(n, target = integer(0), int.level = numeric(0)) {
          ## Error terms, intercepts, and intervention levels
          if (n == 1) {
            Y <- rnorm(node.count(), mean = intercept(), sd = sqrt(err.var()))
          } else {
            Y <- matrix(rnorm(n*node.count(), mean = intercept(), sd = sqrt(err.var())), ncol = n)
          }
          if (length(target) > 0) {
            if (length(int.level) %nin% c(length(target), n*length(target))) {
              stop("int.level must either be a vector of the same length as target, or a matrix of dimension n x length(target)")
            }
            if (is.matrix(int.level)) {
              int.level <- t(int.level)
            }
            if (n == 1) {
              Y[target] <- int.level
            } else {
              Y[target, ] <- int.level
            }
          }

          ## Modified weight matrix (w.r.t. intervention target)
          D <- - t(weight.mat(target))
          diag(D) <- 1.

          ## Calculate results: simulation samples
          result <- solve(D, Y)
          if (n == 1) {
            result
          } else {
            t(result)
          }
        }
        )
    )

#' Coercion from a weight matrix
setAs("matrix", "GaussParDAG",
    def = function(from) {
      p <- nrow(from)
      if (!isAcyclic(from))
        stop("Input matrix does not correspond to an acyclic DAG.")
      edgeL <- lapply(1:p, function(i) which(from[, i] != 0))
      new("GaussParDAG",
          nodes = as.character(1:p),
          in.edges = edgeL,
          param = lapply(1:p, function(i) c(0, 0, from[edgeL[[i]], i])))
    })

#' Coercion from a "graphNEL" object
setAs("graphNEL", "GaussParDAG",
    def = function(from) {
      ## Perform coercion via weight matrix
      A <- as(from, "matrix")
      as(A, "GaussParDAG")
    })

#' Predict interventional or observational data points.  Intervention values
#' must be provided, the value of all non-intervened variables is calculated
#'
#' @param   object    an instance of GaussParDAG
#' @param   newdata   list with two entries:
#'                    target:     list of intervention targets (or single
#'                                intervention target)
#'                    int.level:  list of intervention levels (or single
#'                                vector of intervention levels)
#' @return  a matrix with rows containing the predicted values, or a vector,
#'          if a single prediction is requested
setMethod("predict", "GaussParDAG",
    function(object, newdata) {
      ## Check validity of parameters
      if (!validObject(object))
        stop("The parametric DAG model to be plotted is not valid")
      if (!is.list(newdata$target)) {
        if (is.list(newdata$int.level))
          stop("The two entries of newdata must both be vectors or both be lists.")
        newdata$target <- list(newdata$target)
        newdata$int.level <- list(newdata$int.level)
      }
      stopifnot(is.list(newdata$target),
          is.list(newdata$int.level),
          length(newdata$target) == length(newdata$int.level),
          all(sapply(newdata$target, length) == sapply(newdata$int.level, length)))

      if (length(newdata$target > 1))
        fit <- matrix(0, nrow = length(newdata$target), ncol = object$node.count())
      for (i in 1:length(newdata$target)) {
        ## Calculate predition for i-th target
        y <- object$intercept()
        y[newdata$target[[i]]] <- newdata$int.level[[i]]
        D <- -object$weight.mat(newdata$target[[i]])
        diag(D) <- 1.
        if (length(newdata$target > 1))
          fit[i, ] <- solve(D, y)
        else
          fit <- solve(D, y)
      }

      fit
    })
