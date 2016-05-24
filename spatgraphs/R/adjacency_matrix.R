##############################################################################
#' Transpose sg object
#'
#' This will transpose the adjacency matrix underlying the graph. Will transform
#' to and from sgadj-object (see 'sg2adj')
#'
#' @param x sg-object.
#'
#' @export

t.sg<-function(x)
{
  z<-sg2adj(x)
  z$matrix<-t(z$matrix)
  adj2sg(z)
}


##############################################################################
#' Transpose sgadj object
#'
#' This will transpose the adjacency matrix underlying the graph.
#'
#' @param x sgadj object
#'
#' @export

t.sgadj<-function(x)
{
  x$matrix<-t(x$matrix)
  x
}


###################################################################
#' Symmetrisation of sg adjacency matrix
#' wrapper for 1way and 2way symmetrisation
#'
#' @param x sg object
#' @param way 1: OR rule, 2: AND rule for keeping edges.
#'
#' @export
sg2sym<-function(x, way=1)
{
  if(way==1)# symmetrisize with OR : one direction link is enough
  {
    for(i in 1:length(x$edges) )
      if(length(x$edges[[i]])>0)
        for(j in x$edges[[i]])
          x$edges[[j]]<-union(i,x$edges[[j]])
        x$symmetric<-TRUE
  }
  else # symmetrisize with AND: remove one direction links
  {
    for(i in 1:length(x$edges) )
      if(length(x$edges[[i]])>0)
        for(j in x$edges[[i]])
          if(! (i%in%x$edges[[j]]) )
            x$edges[[i]]<-setdiff(x$edges[[i]],j)
          x$symmetric<-TRUE
  }
  return(x)
}


##############################################################################
#' Make a sparse adjacency matrix from sg-object
#'
#' @param x sg-object
#'
#' @import Matrix
#' @export

sg2sparse<-function(x) {
  ij<-NULL
  for(i in 1:x$N)
    if(length(x$edges[[i]])>0)
      ij<-rbind(ij, cbind(i, x$edges[[i]]))
  sparseMatrix(i=ij[,1], j=ij[,2], dims=c(x$N, x$N))
}

#' Make an sg-object from adjacency matrix
#'
#' @param x square matrix. non-0 elements are taken as edge presence.
#' @export

sparse2sg<-function(x) {
  if(ncol(x)!=nrow(x)) stop("parse2sg: adjacency matrix needs to be a square matrix.")
  edges<-vector("list", ncol(x))
  for(i in 1:ncol(x)){
    edges[[i]]<-which(x[i,]!=0)
  }
  as.sg(edges=edges, type="?", pars=NULL, note="Converted from an unknown matrix")
}


##############################################################################
#' sg to sgadj
#' @param x sg object
#'
#' @export

sg2adj<-function(x)
{
  if(!is(x,"sg")) stop("'x' not of class 'sg'.")
  A<-sg2sparse(x)
  as.sgadj(A, type=x$type, pars=x$parameters)
}

#' sgadj to sg
#'
#' @param x sgadj object
#'
#' @export
adj2sg<-function(x)
{
  if(!is(x,"sgadj")) stop("'x' not of class 'sgadj'.")
  A<-list()
  for(i in 1:x$N)
  {
    A[[i]]<-(1:x$N)[x$matrix[i,]==1]
  }
  as.sg(A, type=x$type, pars=x$parameters, note = "from sgadj-object" )
}

# ##############################################################################
# ## what is this...
# sg2wadj<-function(x)
# {
#   verifyclass(x,"sg")
#   if(is.null(x$weights)) stop("No weights. Run weight.sg(x,...) .")
#   W<-diag(0,x$N)
#   for(i in 1:x$N)
#   {
#     W[i,x$edges[[i]]]<-x$weights[[i]]
#   }
#   sgadj(W, type=x$type, pars=x$parameters, sym=x$symmetric, other="weighted")
# }
##############################################################################
#' Creator for sgadj-class
#' @param edges edge list-of-lists
#' @param type of the graph
#' @param pars parameters for the graph
#' @param other other comments
#'
#' @export
as.sgadj<-function(edges=NULL,type="?",pars=NULL, other="")
{
  e<-list(matrix=edges)
  e$N<-dim(edges)[1]
  e$type<-type
  e$parameters<-pars
  e$other<-other
  class(e)<-"sgadj"
  e
}
##############################################################################
#' print method for sgadj
#'
#' @param x sgadj object
#' @param ... ignored
#'
#' @export
print.sgadj<-function(x,...)
{
  nam<-names(x$parameters)
  p<-"?"

  p<-paste(", par=(",paste(x$parameters,collapse=","),")",sep="")
  cat(paste("'Spatgraphs' ",x$other," adjacency matrix:",
            "\ngraph type '",x$type,"'",p,", for ",x$N," points.\n",sep=""))
  if(!is.null(x$note))cat(paste("Note: ", x$note,".\n",sep=""))

}


##############################################################################
#' plot sgadj
#'
#' @param x sgadj object
#' @param ... passed to plot.sg
#'
#' converts to sg and plots that.
#'
#' @export
plot.sgadj<-function(x, ...)
{
  plot.sg(adj2sg(x),...)
}


##############################################################################
#' weighted sg to weighted adjacency matrix
#'
#' @param x weighted sg object
#' @export

sg2wadj<-function(x)
{
  is_sg(x)
  if(is.null(x$weights)) stop("No weights. Run weight.sg(x,...) .")
  W<-diag(0,x$N)
  for(i in 1:x$N)
  {
    W[i, x$edges[[i]]]<-x$weights[[i]]
  }
  as.sgadj(W, type=x$type, pars=x$parameters, other="weighted")
}




# eof



