## 'value' is a unique named value (compatible with unlist())
## 'simplicesIndices' can be precomputed if they are to be used repeatedly
# uses convhulln for redunddant elim; does not depend on rcdd
redundant.addVeq <- function(vertices, value,
                             simplicesIndices=NULL ## a matrix that describes the facets of the convex hull of the vertices
                             ) { ## a function that adds constraints in an elementary way, then reduces the vertices using qhull
  if (is.null(vertices)) {
    resu <- NULL
  } else if (nrow(vertices)<=ncol(vertices)) { ## convhulln crashes!
    resu <- addSimplexEq(vertices, value)
  } else if (ncol(vertices)==1L) {
    stop.redef("(!) redundant.addVeq called for one-dimensional space.")
    ##if (value>min(vertices) && value<max(vertices)) {resu <- array(value, dim=c(1, 1))} else resu <- NULL
  } else {
    if (is.null(simplicesIndices)) {
      if (sessionInfo()$R.version$`svn rev` < "69993") {
        ## FR->FR should be obsolete some day... load def of capture.output from R devel, future 3.3.0
        capture.output <- temp_capture.output
      } ## else R already has the right capture.output
      abyss <- capture.output(simplicesIndices <- try(convhulln(vertices, "Pp")),type="message")
      if (inherits(simplicesIndices,"try-error")) simplicesIndices <- delaunayn(vertices,"Pp")
      ##                                          redundant, but appears numerically more robust
    }
    blub <- lapply(seq(nrow(simplicesIndices)), function(ii){addSimplexEq(vertices[simplicesIndices[ii, ], ], value)})
    if(length(blub)>1L) { ## list length
      blub <- blub[ ! unlist(lapply(blub,is.null))]  ## removes null elements of list
      resu <- do.call("rbind", blub) ## general case
    } else resu <- blub[[1L]]
  }
  if( ! is.null(resu)) {
    if( ! is.matrix(resu)) dim(resu) <- c(1,length(resu)) ## cf ?as.matrix
    colNames <- colnames(vertices);varName <- names(value) ## hence has original number of columns
    diffNames <- setdiff(colNames, varName)
    colnames(resu) <- colNames;resu <- resu[, diffNames, drop=FALSE] ## now one fewer column
    if (nrow(resu)>ncol(resu)) { ##
      if (ncol(resu)>1L) {
        resu <- resu[unique(as.numeric(convhulln(resu, "Pp"))), , drop=FALSE] ## redundant vertex elimination
      } else resu <- array(c(min(resu), max(resu)), dim=c(2L, 1L))
    }
    colnames(resu) <- diffNames
  }
  return(resu)
}
