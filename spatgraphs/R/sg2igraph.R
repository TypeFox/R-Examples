#' sg to igraph
#'
#' @param x sg object
#'
#' @export

sg2igraph <- function(x) {
  stop("Not implemented. You can use the 'graph_from_adj_list'-function in 'igraph'-package on the edges-element of the graph.\n Example: graph_from_adj_list(x$edges)")
}

# sg2igraph<-function(g, x=NULL)
# {
#   elist<-NULL
#
#   actors <- data.frame(label=as.character((1:g$N)-1))
#
#   for(i in 1:g$N)
#   {
#     if(length(g$edges[[i]]>0))
#     {
#       a<-cbind(i-1,g$edges[[i]]-1)
#       elist<-rbind(elist, a)
#     }
#   }
#
#   elist<-data.frame(from=as.character(elist[,1]), to=as.character(elist[,2]))
#
#   if(!is.null(x))
#   {
#     x <- sg_parse_coordinates(x)
#     actors<-data.frame(names=as.character((1:nrow(x))-1), x=x[,1], y=x[,2])
#     if(ncol(x)==3) actors$z <- x[,3]
#   }
#
#   a<-graph.data.frame(elist, vertices=actors, directed=FALSE)
#   a$N<-g$N
#   a$parameters<-g$parameters
#   a$type<-g$type
#   a
# }
#
# ####################
# igraph2sg<-function(g)
# {
#   b<-get.adjlist(g)
#   for(i in 1:length(b))
#     b[[i]]<-union(b[[i]]+1,NULL)
#
#   as.sg(b, type=g$type, pars=g$parameters, note="coverted from igraph-object")
# }
