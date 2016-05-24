#######################################################################
####
#### querygraph provides unified interface to graph operations.
####
#### Works on graphNEL objects, igraph objects, and adjacency matrices
####
#### Notice: when a graph is returned it is always a graphNEL object
####
#######################################################################

querygraph <-function(object, op, set=NULL, set2=NULL, set3=NULL){

  ## From RBGL / graph packages
  graph.RBGL <-
    c("maxClique",
      "connectedComp",
      "separates",
      "is.triangulated",
      "adj",
      "subgraph",
      "nodes",
      "edges")
  ## From gRbase
  gRbase <-
    c("ancestors",
      "ancestralGraph",
      "ancestralSet",
      "children",
      "closure",
      "edgeList",
      "is.decomposition",
      "is.complete",
      "is.simplicial",
      "parents",
      "simplicialNodes",
      "vpar")

  op <- match.arg(op, choices=c(graph.RBGL, gRbase))
  object <- coerceGraph(object, "graphNEL")

  switch(op,
         ## Functions from graph/RBGL package here.
         "maxClique"=        { RBGL::maxClique(object)$maxCliques               },
         "connectedComp"=    { RBGL::connectedComp(object)                      },
         "separates"=        { RBGL::separates(set, set2, set3, object)         },
         "is.triangulated"=  { RBGL::is.triangulated(object)                    },

         "adj"=              { graph::adj(object, set)                          },
         "subgraph"=         { graph::subGraph(set, object)                     },
         "nodes"=            { graph::nodes(object)                             },
         "edges"=            { graph::edges(object)                             },
         ## gRbase functions
         "ancestors"=,"an"=  { gRbase::ancestors(set, object)		                },
         "ancestralGraph"=   { gRbase::ancestralGraph(set, object)	            },
         "ancestralSet"=     { gRbase::ancestralSet(set, object)                },
         "children"=         { gRbase::children(set, object)         	          },
         "closure"=          { gRbase::closure(set, object)          	          },
         "edgeList"=         { gRbase::edgeList(object)	       		              },
         "is.decomposition"= { gRbase::is.decomposition(set, set2, set3, object) },
         "is.complete"=      { gRbase::is.complete(object, set)         	      },
         "is.simplicial"=    { gRbase::is.simplicial(set, object)         	    },
         "parents"=          { gRbase::parents(set, object)         		        },
         "simplicialNodes"=  { gRbase::simplicialNodes(object)         	        },
         "vpar"=             { gRbase::vpar(object)         			              }
         )
}


########################################################################
###
### Functions which return vectors
###
########################################################################

## adjmat based
ancestors <- function(set, object){
  amat  <- graphNEL2M(object)
  ##if (isUndirectedMAT(amat))
  if (isugMAT_(amat))
    return(NULL)

  An <- setorig <- set
  amat  <- amat[-match(set, rownames(amat)),]

  repeat{
      set2 <- rowSums(amat[,set, drop=FALSE])
      set  <- names(which(set2>0))
      if (!length(set))
          break()
      An <- c(An, set)
      amat  <- amat[set2 == 0,,drop=FALSE]
  }
  setdiff(An, setorig)
}

## adjmat based -- Must be very slow !!!
ancestralSet <- function(set, object){

  amat  <- graphNEL2M(object)
  ##if (isUndirectedMAT(amat))
  if (isugMAT_(amat))
    return(NULL)

  if (missing(set))
    stop("'set' must be given..\n")
  vn   <- colnames(amat)
  an   <- rep(0, length(vn))
  names(an) <- vn
  an[set]   <- 1

  A0 <- set
  repeat{
    x <- amat[,A0,drop=FALSE]
    B <- rownames(x)[apply(x,1,sum)>0]
    if (!length(B))
      break()
    an[B] <- 1
    idx   <- match(A0, colnames(amat))
    amat  <- amat[-idx,-idx,drop=FALSE]
    vn    <- colnames(amat)
    A0    <- intersect(B,vn)
    if (!length(A0))
      break()
  }
  names(an[an>0])
}

## adjmat based
parents <- function(set, object){
  amat  <- graphNEL2M(object)
  ##if (isUndirectedMAT(amat))
  if (isugMAT_(amat))
    return(NULL)

  pa   <- names(which(amat[,set]>0))
  pa   <- setdiff(pa,set)
  if (length(pa)) pa else NULL
}

## graphNEL based
children <- function(set, object){
  if (graph::edgemode(object)=="undirected")
    return(NULL)
  ch <- structure(unlist(graph::edges(object)[set]), names=NULL)
  if (length(ch)) ch else NULL
}

## graphNEL based
closure <- function(set, object){
  uniquePrim(c(set, unlist(graph::adj(object, set))))
}

## graphNEL based
simplicialNodes <- function(object){
  nodes <- graph::nodes(object)
  b     <- unlistPrim(lapply(nodes, function(s) is.simplicial(s, object)))
  sim   <- nodes[b]
  return(sim)
}



########################################################################
###
### Functions which return graphs
###
########################################################################

ancestralGraph <- function(set, object){
  graph::subGraph(ancestralSet(set, object), object)
}

########################################################################
###
### Boolan graph funcions (is.something)
###
########################################################################

is.complete <- function(object, set=NULL){
  if (is.null(set))
    submat <- graphNEL2M(object)
  else
    submat <- graphNEL2M(object)[set,set]
  all(submat[upper.tri(submat)]>0)
}

is.decomposition <- function(set, set2, set3, object){
  vn <- uniquePrim(c(set, set2, set3))
  if (setequal(vn, graph::nodes(object))){
    RBGL::separates(set, set2, set3, object) & is.complete(object, set3)
  } else {
    FALSE
  }
}

is.simplicial <- function(set, object){
  x <- unlist(graph::adj(object,set))
  is.complete(graph::subGraph(x, object), x)
}
