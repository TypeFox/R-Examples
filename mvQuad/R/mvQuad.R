
# !diagnostics style=[true]
# !diagnostics level=[all]



# GetPredefinedTypes <- function(print.types = FALSE){
#    rule.names <- names(QuadRules)
#    rules <- vector(mode = "list", length = length(rule.names))
#
#    for (i in seq_along(rule.names)){
#       tmp <- NULL
#       for (j in seq_len(length(QuadRules[[rule.names[i]]]))){
#          len <- length(QuadRules[[rule.names[i]]][[j]][[2]])
#          if (len==0) len <- NA
#          tmp <- c(tmp, len)
#       }
#       rules[[i]] <- list(levels = tmp)
#       if (print.types){
#          cat("-----", rule.names[i], "---------------------- \n",
#              "levels: \n ",
#              tmp, "\n\n" )
#       }
#    }
#    names(rules) <- rule.names
#    invisible(rules)
# }

.hardCodedTypes <- c("cNC1", "cNC2", "cNC3", "cNC4", "cNC5", "cNC6",
                    "oNC0", "oNC1", "oNC2", "oNC3")

.extGaussQuad <- c("GLe", "GLa", "GHe", "GHN")

.preDefinedTypes <- structure(list(nLe = structure(list(levels = c(1L, 3L, 3L, 7L, 7L, 7L, 15L,
                                                                   15L, 15L, 15L, 15L, 15L, 31L, 31L, 31L, 31L, 31L, 31L, 31L,
                                                                   31L, 31L, 31L, 31L, 31L, 63L)), .Names = "levels"),
                                   GKr = structure(list(levels = c(NA, NA, 3L, NA, 5L, NA, 7L, NA, 9L, NA, 11L,
                                                                   NA, 13L, NA, 15L, NA, 17L, NA, 19L, NA, 21L, NA, 23L,
                                                                   NA, 25L, NA, 27L, NA, 29L)), .Names = "levels"),
                                   nHe = structure(list(levels = c(1L, 3L, 3L, 7L, 9L, 9L, 9L, 9L, 17L, 19L,
                                                        19L, 19L, 19L, 19L, 19L, 31L, 33L, 35L, 35L, 35L, 35L,
                                                        35L, 35L, 35L, 35L)), .Names = "levels"),
                                   nHN = structure(list(levels = c(1L, 3L, 3L, 7L, 9L, 9L, 9L, 9L, 17L, 19L,
                                                                   19L, 19L, 19L, 19L, 19L, 31L, 33L, 35L, 35L, 35L, 35L,
                                                                   35L, 35L, 35L, 35L)), .Names = "levels"),
                                   Leja = structure(list(levels = 1:141), .Names = "levels")),
                              .Names = c("nLe", "GKr", "nHe", "nHN", "Leja"))


#' reads a quadrature-rule from a text file
#'
#' \code{readRule} reads a quadrature-rule from a text file
#'
#' @param file file name of the text file containing the quadrature rule
#'
#' @details The text file containing the quadrature rule has to be formatted in the following way:
#'
#' The first line have to declare the domain \code{initial.domain a b}, where a and b denotes the lower and upper-bound for the integration domain.
#' This can be either a number or '-Inf'/'Inf' (for example \code{initial.domain 0 1} or \code{initial.domain 0 Inf})
#'
#' Every following line contains one single node and weight belonging to one level of the rule (format: 'level' 'node' 'weight').
#' This example shows the use for the "midpoint-rule" (levels: 1 - 3).
#'
#' > \code{initial.domain 0 1}
#'
#' > \code{1 0.5 1}
#'
#' > \code{2 0.25 0.5}
#'
#' > \code{2 0.75 0.5}
#'
#' > \code{3 0.166666666666667 0.333333333333333}
#'
#' > \code{3 0.5 0.333333333333333}
#'
#' > \code{3 0.833333333333333 0.333333333333333}
#'
#'
#' @return Returns an object of class 'customRule', which can be used for creating a 'NIGrid' (\code{\link{createNIGrid}})
#'
#' @seealso
#' \code{\link{createNIGrid}}
#'
#' @examples
#' \dontrun{myRule <- readRule(file="midpoint_rule.txt")}
#' \dontrun{nw <- createNIGrid(d=1, type = myRule.txt, level = 2)}

readRule <- function(file=NULL){
  if (is.null(file)) stop("filename necessary")

  tmp <- readLines(file)

  tmp.index <- which(grepl("initial.domain", tmp))
  if (length(tmp.index)!=1) stop("no or more than one 'initial.domain' defined in your file")
  initial.domain <- matrix(as.numeric(strsplit(tmp[tmp.index], split = " ", fixed = T)[[1]][-1]), ncol=2)
  tmp <- strsplit(tmp, split = " ", fixed = T)
  tmp <- do.call(rbind, tmp)
  colnames(tmp) <- c("l", "n", "w")

  suppressWarnings(levels <- unique(na.omit(as.numeric(tmp[,"l"]))))

  rule <- vector(mode = "list", length = max(levels))

  for (i in levels){
    tmp.l <- subset(tmp, tmp[,1]==i)
    rule[[i]] <- list(n = as.numeric(tmp.l[,2]),
                      w = as.numeric(tmp.l[,3]),
                      features = list(initial.domain=initial.domain))
  }

  class(rule) <- "costumRule"
  return(rule)
}


.grid1D <- function(type, level=1) {

   if (class(type)=="costumRule"){
     if (!(level %in% c(1:length(type)))){
       stop(paste("degree (",level,") for user defined rule not supported \n") )
     } else {
       return(type[[level]])
     }
   }

  if (class(type)!="character") stop("type of quadrature rule not appropriate defined")

   if (!((type %in% c(.hardCodedTypes, names(.preDefinedTypes), .extGaussQuad)))) {
      if (!existsFunction(type)) {
         stop("type of quadrature rule is not supported")
      } else {
         tmp.fun <- match.fun(type)
         return(tmp.fun(level))
      }
   }

   if (type %in% .hardCodedTypes) {
      quad.type <- substr(type,1,3)

      # closed-Newton-Cotes forumla
      if (quad.type=="cNC") {

         degree <- as.numeric(substr(type,4,4))
         if (!(degree %in% c(1:6)) ) stop("degree of closed-Newton-Cotes formula must between 1 and 6")

         if (degree == 1) weights <- c(1/2, 1/2)
         if (degree == 2) weights <- c(1/6, 4/6, 1/6)
         if (degree == 3) weights <- c(1/8, 3/8, 3/8, 1/8)
         if (degree == 4) weights <- c(7/90, 32/90, 12/90, 32/90, 7/90)
         if (degree == 5) weights <- c(19/288, 75/288, 50/288, 50/288, 75/288, 19/288)
         if (degree == 6) weights <- c(41/840, 216/840, 27/840, 272/840, 27/840, 216/840, 41/840)

         interv_length <- 1/level
         weights <- weights * interv_length
         steps <- interv_length/degree
         n <- seq(0,1,steps)
         w <- numeric(length(n))

         for (i in 1:level) {
            w[((i-1)*degree+1):(i*degree+1)] <- w[((i-1)*degree+1):(i*degree+1)] + weights
         }

         return(list(n = n,w = w, features=list(initial.domain=matrix(c(0,1), ncol = 2))))
      }

      # open-Newton-Cotes Forumla
      if (quad.type=="oNC") {
         degree <- as.numeric(substr(type,4,4))
         if (!(degree %in% c(0:3)) ) stop("degree of open-Newton-Cotes formula must between 0 and 3")

         if (degree == 0) {
            weights <- c(1)
            nodes <- c(1/2)
         }
         if (degree == 1) {
            weights <- c(1/2, 1/2)
            nodes <- c(1/4, 3/4)
         }
         if (degree == 2) {
            weights <- c(3/8, 2/8, 3/8)
            nodes <- c(1/6, 1/2, 5/6)
         }
         if (degree == 3) {
            weights <- c(13/48, 11/48, 11/48, 13/48)
            nodes <- c(1/8, 3/8, 5/8, 7/8)
         }
         interv_length <- 1/level
         weights <- weights * interv_length
         nodes <- nodes * interv_length

         n = numeric((degree+1)*level)
         w = rep(weights,level)

         for (i in 1:level){
            n[((i-1)*(degree+1)+1):(i*(degree+1))] <- (i-1)*interv_length + nodes
         }

         return(list(n = n, w = w, features=list(initial.domain=matrix(c(0,1), ncol = 2))))
      }
   }

  if (type %in% names(.preDefinedTypes)) {
    if (!(level %in% c(1:length(QuadRules[[type]]))) ){
      stop(paste("degree (",level,") for rule '", type, "' not supported \n") )
    } else {
      return(QuadRules[[type]][[level]])
    }
  }

  if (type %in% .extGaussQuad) {
    if (type=="GLe") {
      tmp.nw <- gauss.quad(level, kind="legendre")
      tmp.nw$nodes <- tmp.nw$nodes / 2 + 0.5
      tmp.nw$weights <- tmp.nw$weights / 2
      tmp.dom <- list(initial.domain=matrix(c(0,1), ncol = 2))
    }

    if (type=="GLa") {
      tmp.nw <- gauss.quad(level, kind="laguerre")
      tmp.dom <- list(initial.domain=matrix(c(0,Inf), ncol = 2))
    }

    if (type=="GHe") {
      tmp.nw <- gauss.quad(level, kind="hermite")
      tmp.nw$nodes <- tmp.nw$nodes * sqrt(2)
      tmp.nw$weights <- tmp.nw$weights * sqrt(2) / sqrt(exp(-(tmp.nw$nodes^2)))
      tmp.dom <- list(initial.domain=matrix(c(-Inf, Inf), ncol = 2))
    }
    if (type=="GHN") {
      tmp.nw <- gauss.quad(level, kind="hermite")
      tmp.nw$nodes <- tmp.nw$nodes * sqrt(2)
      tmp.nw$weights <- tmp.nw$weights * sqrt(2) / sqrt(exp(-(tmp.nw$nodes^2)))*dnorm(tmp.nw$nodes)
      tmp.dom <- list(initial.domain=matrix(c(-Inf, Inf), ncol = 2))
    }

    return(list(n = as.numeric(tmp.nw$nodes), w = as.numeric(tmp.nw$weights), features=tmp.dom))
  }

}

.fgridnD <- function(dim, type, level) {
  le <- numeric(dim)
  domain <- matrix(NA, nrow = dim, ncol = 2)
  storage <- vector(mode = "list", length = dim)
  for (i in 1:dim) {
    tmp <- .grid1D(type[[i]], level[i])
    le[i] <- length(tmp$w)
    domain[i, ] <- tmp$features$initial.domain
    storage[[i]] <- tmp
  }

  No.Points <- prod(le)
  w <- rep(1, No.Points)
  n <- matrix(NA, nrow = No.Points, ncol = dim)

  for (i in 1 :dim) {
    storage[[i]]
    n[ ,i] = rep(storage[[i]]$n, each = prod(le[1:i - 1]), times = (No.Points/prod(le[1:i])))
    w = w * rep(storage[[i]]$w, each = prod(le[1:i - 1]), times = (No.Points/prod(le[1:i])))
  }
  return(list(n = n, w = w, features = list(initial.domain = domain)))
}

.SpGrSeq <- function(dimension, norm, level.trans){
  seq.vec       <- rep(0, dimension)
  a             <- norm - dimension
  seq.vec[ 1 ]  <- a
  fs            <- matrix(NA, nrow=choose(norm - 1, dimension - 1), ncol=dimension )
  fs[1, ]       <- seq.vec
  cnt           <- 1
  cnt.seq       <- 1

  while ( seq.vec[ dimension ] < a ) {
    if ( cnt == dimension ) {
      for (i in seq(cnt-1, 1, -1) ) {
        cnt <- i
        if ( seq.vec[ i ] != 0 ) { break }
      }
    }
    seq.vec[ cnt ]  <- seq.vec[ cnt ] - 1
    cnt             <- cnt + 1
    seq.vec[ cnt ]  <- a - sum( seq.vec[1:(cnt-1)] )
    if ( cnt < dimension ) {
      seq.vec[ (cnt+1):dimension ] <- rep(0, dimension - cnt)
    }
    cnt.seq <- cnt.seq + 1
    fs[cnt.seq, ] = seq.vec
  }

  fs <- fs + 1
  fs <- level.trans(fs)

  return(fs)
}


.sgridnD <- function(dim, type, level, level.trans){
  minq <- max(0,level-dim)
  maxq <- level-1

  storage <- vector(mode="list", length= sum(choose(c(minq:maxq) + dim - 1, dim - 1)))

  cnt <- 0
  for (q in (minq:maxq)) {
    bq  <- (-1)^(maxq-q)*choose(dim-1,dim+q-level)

    GridSeq <- .SpGrSeq(dim,dim+q,level.trans)
    for (i in 1:dim(GridSeq)[1]){
      res <- .fgridnD(dim,type,GridSeq[i,])

      cnt <- cnt + 1
      storage[[cnt]] <- data.table::data.table(res[[1]], bq*res[[2]])
    }
  }
  tmp <- data.table::rbindlist(storage)
  tmp.names <- c(paste("n", c(1:dim), sep="_"), "w.c")
  data.table::setnames(tmp, tmp.names)
  data.table::setkeyv(x = tmp, cols = tmp.names[-(dim+1)])

  w.c <- NULL

  tmp[, w:=sum(w.c), by = eval(tmp.names[-(dim+1)])]
  tmp <- unique(tmp, by = eval(tmp.names[-(dim+1)]))

  tmp <- tmp[w!=0]

  n <- tmp[,tmp.names[-(dim+1)], with=FALSE]
  w <- tmp[, w]

  row.names(n) <- NULL
  colnames(n) <- NULL

  return(list(n = as.matrix(n), w = as.matrix(w), features=res[[3]]))
}

#' creates a grid for numerical integration.
#'
#' \code{createNIGrid} Creates a grid for multivariate numerical integration.
#' The Grid can be based on different quadrature- and construction-rules.
#'
#' @param dim number of dimensions
#' @param type quadrature rule (see Details)
#' @param level accuracy level (typically number of grid points for the underlying 1D quadrature rule)
#' @param ndConstruction character vector which denotes the construction rule
#' for multidimensional grids.
#' \describe{
#' \item{\code{product}}{for product rule, returns a "full grid" (default)}
#' \item{\code{sparse}}{for combination technique, leads to a regular "sparse grid".}}
#'
#' @param level.trans logical variable denotes either to take the levels as number
#' of grid points (FALSE = default) or to transform in that manner that number of
#' grid points = 2^(levels-1) (TRUE). Alternatively \code{level.trans} can be a function, which takes (n x d)-matrix and returns
#' a matrix with the same dimensions (see the example; this feature is particularly useful for the 'sparse' construction rule,
#' to account for different importance of the dimensions).
#'
#' @details The following quadrature rules are supported (build-in).
#' \describe{
#'  \item{\code{cNC1, cNC2, ..., cNC6}}{closed Newton-Cotes Formula of degree 1-6 (1=trapezoidal-rule; 2=Simpson's-rule; ...),
#'  initial interval of integration: [0, 1]}
#'  \item{\code{oNC0, oNC1, ..., oNC3}}{open Newton-Cote Formula of degree 0-3 (0=midpoint-rule; ...),
#'  initial interval of integration: [0, 1]}
#'  \item{\code{GLe, GKr}}{Gauss-Legendre and Gauss-Kronrod rule for an initial interval of integration: [0, 1]}
#'  \item{\code{nLe}}{nested Gauss-Legendre rule for an initial interval of integration: [0, 1] (Knut Petras (2003). Smolyak cubature of given polynomial degree with few nodes for increasing dimension. Numerische Mathematik 93, 729-753)}
#'  \item{\code{GLa}}{Gauss-Laguerre rule for an initial interval of integration: [0, INF)}
#'  \item{\code{GHe}}{Gauss-Hermite rule for an initial interval of integration: (-INF, INF)}
#'  \item{\code{nHe}}{nested Gauss-Hermite rule for an initial interval of integration: (-INF, INF) (A. Genz and B. D. Keister (1996). Fully symmetric interpolatory rules for multiple integrals over infinite regions with Gaussian weight." Journal of Computational and Applied Mathematics 71, 299-309)}
#'  \item{\code{GHN}, \code{nHN}}{(nested) Gauss-Hermite rule as before but weights are multiplied by the standard normal density (\eqn{\hat(w)_i = w_i * \phi(x_i)}).}
#'  \item{\code{Leja}}{Leja-Points for an initial interval of integration: [0, 1]}}
#' The argument \code{type} and \code{level} can also be vector-value, different for each dimension (the later only for "product rule"; see examples)
#'
#' @return Returns an object of class 'NIGrid'. This object is basically an environment
#' containing nodes and weights and a list of features for this special grid. This
#' grid can be used for numerical integration (via \code{\link{quadrature}})
#' @references
#' Philip J. Davis, Philip Rabinowitz (1984): Methods of Numerical Integration
#'
#' F. Heiss, V. Winschel (2008): Likelihood approximation by numerical integration on sparse grids, Journal of Econometrics
#'
#' H.-J. Bungartz, M. Griebel (2004): Sparse grids, Acta Numerica
#' @seealso
#' \code{\link[=rescale.NIGrid]{rescale}}, \code{\link{quadrature}}, \code{\link[=print.NIGrid]{print}}, \code{\link[=plot.NIGrid]{plot}} and \code{\link[=size.NIGrid]{size}}
#' @examples
#' ## 1D-Grid --> closed Newton-Cotes Formula of degree 1 (trapeziodal-rule)
#' myGrid <- createNIGrid(dim=1, type="cNC1", level=10)
#' print(myGrid)
#' ## 2D-Grid --> nested Gauss-Legendre rule
#' myGrid <- createNIGrid(dim=2, type=c("GLe","nLe"), level=c(4, 7))
#' rescale(myGrid, domain = rbind(c(-1,1),c(-1,1)))
#' plot(myGrid)
#' print(myGrid)
#' myFun <- function(x){
#'    1-x[,1]^2*x[,2]^2
#' }
#' quadrature(f = myFun, grid = myGrid)
#' ## level transformation
#' levelTrans <- function(x){
#'   tmp <- as.matrix(x)
#'   tmp[, 2] <- 2*tmp[ ,2]
#'   return(tmp)
#' }
#' nw <- createNIGrid(dim=2, type="cNC1", level = 3,
#'    level.trans = levelTrans, ndConstruction = "sparse")
#' plot(nw)
createNIGrid <- function(dim=NULL, type=NULL, level=NULL,
                         ndConstruction="product", level.trans=NULL){
## to-do:
#     - adaptive grid (CW, 2015-01-05)

   if (is.null(dim) | (dim%%1)!=0 | dim < 1) {
      stop("'dim' is not appropriate defined")
   }

   if (class(type)!="list"){
     if (class(type)=="character") type <- as.list(type)
     if (class(type)=="costumRule") type <- list(type)
   }

   if (is.null(type) | (length(type) > 1 & length(type) < dim)) {
      stop("'type' is not appropriate defined")
   }

   if (length(type) < dim){
      type <- replicate(type[[1]], n=dim, FALSE)
   }

   if (is.null(level) | any(level%%1!=0)){
      stop("'level' is not appropriate defined")
   }

   if (length(level) < dim){
      level <- rep(level, dim)
   }

   level <- matrix(level, ncol = dim)

   if (is.logical(level.trans)) {
      if (level.trans == TRUE){
         level.trans <- function(x){2^(x-1)}
      }
   }

   if (!is.function(level.trans)){
      level.trans <- function(x){x}
   }

   if (dim > 1){
      if (ndConstruction=="product"){
         nw <- .fgridnD(dim = dim, type = type, level = level.trans(level))
      }
      if (ndConstruction=="sparse"){
         nw <- .sgridnD(dim = dim, type = type, level = level[1], level.trans = level.trans)
      }
   } else {
      nw <- .grid1D(type = type[[1]], level = level.trans(level))
   }

   object=new.env(parent=globalenv())

   object$dim <- dim

   type.text <- lapply(type,
                       function(tmp){
                         if (class(tmp)=="costumRule") return("costum")
                         if (class(tmp)=="character")  return(tmp)
                       })
   object$type <- unlist(type.text)
   object$level <- level
   object$level.trans <- level.trans

   object$ndConstruction <- ndConstruction

   if (exists("features", where=nw)){
     object$features <- c(type="static", move=0L, nw[["features"]])
   } else {
     object$features <- list(type="static", move=0L, initial.domain = matrix(c(NA, NA), ncol=2))
   }

   if (!exists("n", where=nw) | !exists("n", where=nw)) stop("this quadrature rule doesn't provide nodes (or/and) weights")

   object$nodes <- as.matrix(nw[["n"]])
   object$weights <- as.matrix(nw[["w"]])
   class(object) <- "NIGrid"
   return(object)
}

#' copies an NIGrid-object
#'
#' \code{copyNIGrid} copies an NIGrid-object
#'
#' @param object1 original NIGrid-object
#' @param object2 destination; if NULL \code{copyNIGrid} returns a NIGrid-object
#' otherwise the \code{object2} will be overwritten.
#'
#' @return Returns a NIGrid-object or NULL
#' @examples
#' myGrid <- createNIGrid(dim=2, type="GHe", level=5)
#' myGrid.copy <- copyNIGrid(myGrid)
copyNIGrid <- function(object1, object2=NULL){
   if (is.null(object2)) {
      object2 = new.env(parent = globalenv())
      class(object2) = class(object1)
      nullFlag = TRUE
   }
   elements = ls(envir = object1)
   for (index in 1:length(elements)) {
      assign(elements[[index]], get(elements[[index]], envir = object1,
                                  inherits = FALSE), envir = object2)
   }
   if (nullFlag) {
      return(object2)
   } else {
      return(NULL)
   }
}

#' @name getNodes and getWeights
#' @rdname getNodesWeights
#'
#' @title get nodes and weights from an NIGrid-object
#' @description \code{getNodes} and \code{getWeights} extract the (potentially rescaled) nodes and weights
#'  out of an NIGrid-Object
#'
#' @param grid object of class \code{NIGrid}
#' @return Returns the nodes or weights of the given grid
#' @seealso
#' \code{\link{createNIGrid}}
#' @examples
#' myGrid <- createNIGrid(dim=2, type="cNC1", level=3)
#' getNodes(myGrid)
getNodes <- function(grid){
   if (!is(grid, "NIGrid")) { stop(" 'object' argument must be of class 'NIGrid' .") }
   if (grid$features$move==0) return(grid$nodes)

   if (grid$features$move==1) {
      n <- grid$nodes
      for (d in 1:grid$dim){
         n[,d] <- (grid$nodes[,d]-grid$features$initial.domain[d,1])*
            diff(grid$features$domain[d,])/diff(grid$features$initial.domain[d,]) +
            grid$features$domain[d,1]
      }
   }

   if (grid$features$move==2) {
      C <- as.matrix(grid$features$C)
      m <- grid$features$m


      if (grid$dim > 1) {
        if (grid$features$dec.type==0) {
          A <- diag(sqrt(diag(C)))
        }
        if (grid$features$dec.type==1) {
          res.dec <- eigen(C)
          A <- res.dec[[2]] %*% diag(sqrt(res.dec[[1]]))
        }
        if (grid$features$dec.type==2) {
          A <- t(chol(C))
        }
      } else {
        A <- as.matrix(sqrt(C))
      }

      n <- t(A %*% t(grid$nodes))
      for (d in 1:grid$dim){
         n[, d] <- n[, d] + m[d]
      }



   }
   return(n)
}


#' @rdname getNodesWeights
#' @examples
#' getWeights(myGrid)
getWeights <- function(grid){
   if (!is(grid, "NIGrid")) { stop(" 'object' argument must be of class 'NIGrid' .") }
   if (grid$features$move==0) return(grid$weights)

   if (grid$features$move==1) {
      w <- grid$weights
      for (d in 1:grid$dim){
         w <- w * diff(grid$features$domain[d,])/diff(grid$features$initial.domain[d,])
      }
   }

   if (grid$features$move==2) {
      C <- as.matrix(grid$features$C)

      if (grid$features$dec.type==0) {
         w <- grid$weights * prod(sqrt(diag(C)))
      }
      if (grid$features$dec.type==1) {
         res.dec <- eigen(C)
         w <- grid$weights * prod(sqrt(res.dec[[1]]))
      }
      if (grid$features$dec.type==2) {
         A <- t(chol(C))
         w <- grid$weights * prod(diag(A))
      }
   }

   return(w)
}

#' @name size (size.NIGrid)
#' @rdname size.NIGrid
#'
#' @title returns the size of an NIGrid-object
#' @description Returns the size of an NIGrid-object
#' @param ... other arguments passed to the specific method
#'
#' @return Returns the grid size in terms of dimensions, number of grid points and used memory
#' @examples
#' myGrid <- createNIGrid(dim=2, type="GHe", level=5)
#' size(myGrid)
#' @export
size <- function(object, ...) UseMethod("size")

#' @rdname size.NIGrid
#' @param object a grid of type \code{NIGrid}
#' @export
size.NIGrid <- function(object, ...){
   size.dim <- object$dim
   size.no  <- length(object$weights)
   size.mem <- object.size(object$nodes)+object.size(object$weights)
   size.object <- list(dim=size.dim, gridpoints=size.no, memory=size.mem)
   class(size.object) <- "NIGrid.size"
   return(size.object)
}

.print.NIGrid.size <- function(x, ...){
   cat("grid size \n",
       " # dimensions: ", x[[1]],
       "\t # gridpoints: ", x[[2]],
       "\t used memory:", format(x[[3]], units="auto"), "\n", sep="")
}

#' @rdname size.NIGrid
#' @param x object of type \code{NIGrid}
#' @examples
#' dim(myGrid)
#' @export
dim.NIGrid <- function(x){
   return(x$dim)
}


#' @rdname print.NIGrid
#' @name print (print.NIGrid)
#' @title prints characteristic information for an NIGrid-object
#' @description Prints characteristic information for an NIGrid-object
#' @param x a grid of type \code{NIGrid}
#' @param ... further arguments passed to or from other methods
#'
#' @return Prints the information for an NIGrid-object (i.a. grid size (dimensions, grid points, memory usage),
#' type and support)
#' @examples
#' myGrid <- createNIGrid(dim=2, type="GHe", level=5)
#' print(myGrid)
#' @export
print.NIGrid <- function(x, ...){
   tmp <- size(x)
   cat("Grid for Numerical Integration  \n",
       " # dimensions: ", tmp[[1]],
       "\t # gridpoints: ", tmp[[2]],
       "\t used memory:", format(tmp[[3]], units="auto"), "\n\n", sep="")

   if (x$dim > 1) cat("nD-Construction principle:  '", x$ndConstruction, "'\n", sep="")
   tmp.level <- as.vector(x$level)
   if (length(unique(x$type))==1 & length(unique(tmp.level))==1) {
      cat(" same quadrature rule for each dimension \n")
      cat(" \t type: ", unique(x$type), "\n")
      cat(" \t level: ", unique(tmp.level), "\n")
      cat(" \t initial domain: ", x$features$initial.domain[1,])
   } else {
      cat(" individual quadrature rule for each dimension \n")
      for (d in 1:x$dim){
         cat(" dim =", d, "--> \t type:", x$type[d], "\t level: ",
             tmp.level[d], "\t initial domain:", x$features$initial.domain[d,], "\n")
      }
   }

   if (x$features$move!=0){
      if (x$features$move == 1){
         cat("\n Grid rescaled. New Domain: \n")
         tmp <- as.matrix(x$features$domain)
         colnames(tmp) <- c("a", "b")
         rownames(tmp) <- paste("dim ", 1:x$dim, ": ", sep="")
         print(tmp)
      }

      if (x$features$move == 2){
         cat("\n Grid rescaled.")
         cat("\n mean-vector: ", x$features$m)
         cat("\n covariance matrix: \n")
         print(x$features$C)
      }
   }
   cat("\n")
}

#' @rdname plot.NIGrid
#' @name plot (plot.NIGrid)
#' @title plots an NIGrid-object
#' @description
#' Plots the grid points of an NIGrid-object
#'
#' @param x a grid of type \code{NIGrid}
#' @param plot.dimension vector of length 1, 2 or 3. with the dimensions to be plotted
#' (see examples)
#' @param ... arguments passed to the default plot command
#'
#' @examples
#' myGrid <- createNIGrid(dim=4, type=c("GHe", "cNC1", "GLe", "oNC1"),
#'                        level=c(3,4,5,6))
#' plot(myGrid) ## dimension 1-min(3,dim(myGrid)) are plotted
#' ## Free arranged plots
#' plot(myGrid, plot.dimension=c(4,2,1))
#' plot(myGrid, plot.dimension=c(1,2))
#' plot(myGrid, plot.dimension=c(3))
#' @export
plot.NIGrid <- function(x, plot.dimension=NULL, ...){
   p.dim <- 1:min(x$dim, 3)

   arg.list <- list(...)

   if (!is.null(plot.dimension)){
      if (length(plot.dimension) >= 4) {
         cat("only first, second and third stated dimension is used")
         plot.dimension <- as.numeric(plot.dimension[1:3])
      }

      if (any(!(plot.dimension %in% (1:x$dim)))){
         stop("wrong dimensions stated to be plotted")
      }
      p.dim <- plot.dimension
   }

   if (length(p.dim)==1){
      x.tmp <- as.numeric(getNodes(x)[, p.dim])
      y.tmp <- rep(1,length(x.tmp))
      if (length(which("xlab" == names(arg.list))) == 0) arg.list = c(arg.list, xlab=paste("dim", p.dim))
      if (length(which("axes" == names(arg.list))) == 0) {
         arg.list = c(arg.list, axes=FALSE)
         createXaxis = TRUE
      }
      do.call(plot, c(list(x=x.tmp, y=y.tmp, ylab=""), arg.list))
      if (createXaxis) axis(side = 1)
   }
   if (length(p.dim)==2){
      tmp <- getNodes(x)
      x.tmp <- as.numeric(tmp[, p.dim[1]])
      y.tmp <- as.numeric(tmp[, p.dim[2]])

      if (length(which("xlab" == names(arg.list))) == 0) arg.list = c(arg.list, xlab=paste("dim", p.dim[1]))
      if (length(which("ylab" == names(arg.list))) == 0) arg.list = c(arg.list, ylab=paste("dim", p.dim[2]))
      do.call(plot, c(list(x=x.tmp, y=y.tmp), arg.list))
   }
   if (length(p.dim)==3){
      tmp <- getNodes(x)
      x.tmp <- as.numeric(tmp[, p.dim[1]])
      y.tmp <- as.numeric(tmp[, p.dim[2]])
      z.tmp <- as.numeric(tmp[, p.dim[3]])
      if (length(which("xlab" == names(arg.list))) == 0) arg.list = c(arg.list, xlab=paste("dim", p.dim[1]))
      if (length(which("ylab" == names(arg.list))) == 0) arg.list = c(arg.list, ylab=paste("dim", p.dim[2]))
      if (length(which("zlab" == names(arg.list))) == 0) arg.list = c(arg.list, zlab=paste("dim", p.dim[3]))
      do.call(rgl::plot3d, c(list(x=x.tmp, y=y.tmp, z=z.tmp), arg.list))
   }
}

#' @rdname rescale.NIGrid
#' @name rescale (rescale.NIGrid)
#' @title moves, rescales and/or rotates a multidimensional grid.
#' @description
#' \code{rescale.NIGrid} manipulates a grid for more efficient numerical integration with
#' respect to a given domain (bounded integral) or vector
#' of means and covariance matrix (unbounded integral).
#'
#' @param object an initial grid of type \code{NIGrid}
#' @param ... further arguments passed to or from other methods
#' @return This function modifies the "support-attribute" of the grid. The
#' recalculation of the nodes and weights is done when the \code{\link{getNodes}} or \code{\link{getWeights}}
#' are used.
#' @references Peter Jaeckel (2005): A note on multivariate Gauss-Hermite quadrature
#' @seealso
#' \code{\link{quadrature}}, \code{\link{createNIGrid}}
#' @export
rescale <- function(object, ...) UseMethod("rescale")

#' @rdname rescale.NIGrid
#' @param domain a (d x 2)-matrix with the boundaries for each dimension
#' @param m vector of means
#' @param C covariance matrix
#' @param dec.type type of covariance decomposition (\cite{Peter Jaeckel (2005)})
#' @examples
#' C = matrix(c(2,0.9,0.9,2),2)
#' m = c(-.5, .3)
#' par(mfrow=c(3,1))
#'
#' myGrid <- createNIGrid(dim=2, type="GHe", level=5)
#'
#' rescale(myGrid, m=m, C=C, dec.type=0)
#' plot(myGrid, col="red")
#'
#' rescale(myGrid, m=m, C=C, dec.type=1)
#' plot(myGrid, col="green")
#'
#' rescale(myGrid, m=m, C=C, dec.type=2)
#' plot(myGrid, col="blue")
#' @export
rescale.NIGrid <- function(object, domain=NULL, m=NULL, C=NULL, dec.type=0, ...){
   # check plausibility --> bounded vs. unbounded integral (dimension wise)

   if (is.null(domain) & is.null(m) & is.null(C)) {
      object$features[["move"]] <- 0L
      object$features[["domain"]] <- NULL
      object$features[["m"]] <- NULL
      object$features[["C"]] <- NULL
      object$features[["dec.type"]] <- NULL
      return(invisible(NULL))
   }

   if (!is.null(domain) & (!is.null(m) | !is.null(C))) stop("rescaling not appropriate defined")
   if (is.null(domain) & (is.null(m) | is.null(C))) stop("rescaling not appropriate defined")

   if (!is.null(domain)){
      if (any(is.infinite(object$features[["initial.domain"]]))) stop("rescaling not appropriate defined (bounded domain for unbounded rule)")
      object$features[["move"]] <- 1L
      object$features[["domain"]] <- matrix(domain, ncol=2)
      return(invisible(NULL))
   }

   if (is.null(domain)) {
      C <- as.matrix(C)
      C.dim <- dim(C)
      if (C.dim[1]==C.dim[2]) {
         C.dim <- C.dim[1]
      } else stop("(in rescale) covariance matrix C is not quadratic")

      if (!(all(C == t(C)) & all(eigen(C)$values > 0))) {
         # check for symmetry and positive semi definitness
         stop("(in rescale) no appropriate covariance matrix C")
      }

      if (!length(m)==C.dim) stop("(in rescale) vector of mean not appropriate")

      if (!object$dim==C.dim) stop("(in rescale) initial grid not appropriate")
      if (C.dim==1) dec.type<-0


      object$features[["move"]] <- 2L
      object$features[["m"]] <- m
      object$features[["C"]] <- C
      object$features[["dec.type"]] <- dec.type
      return(invisible(NULL))
   }
}


#' @title computes the approximated Integral
#' @description
#'  \code{quadrature} computes the integral for a given function based on an NIGrid-object
#'
#' @param f a function which takes the x-values as a (n x d) matrix as a first argument
#' @param grid a grid of type \code{NIGrid}
#' @param ... further arguments for the function \code{f}
#' @return The approximated value of the integral
#' @seealso \code{\link{createNIGrid}}, \code{\link[=rescale.NIGrid]{rescale}}
#' @examples
#' myGrid <- createNIGrid(dim=2, type="GLe", level=5)
#' rescale(myGrid, domain=rbind(c(-1,1),c(-1,1)))
#' plot(myGrid, col="blue")
#' myFun <- function(x){
#'    1 - x[,1]^2 * x[,2]^2
#' }
#' quadrature(myFun, myGrid)
quadrature <- function(f, grid=NULL, ...){
## to-do:
#    - add adaptivity (CW; 2015-01-05)
#    - add error estimate (CW; 2015-01-05)
   f <- match.fun(f)
   ff <- function(x) f(x, ...)

   if (is.null(grid)) stop("quadrature not appropriate specified")

   if (!is.null(grid) & !is(grid, "NIGrid")) stop("'Grid' isn't of class NIGrid")

   # place to implement adaptivity
   y <- ff(getNodes(grid))
   w <- getWeights(grid)
   if (length(y) != length(w)) stop("function 'f' fits not to the supported Grid")
   return(sum(y*w))
}
