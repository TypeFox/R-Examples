#' makeAVS
#'
#' This function will create disjoint LD blocks from a GRanges object
#' @param x A GRanges object which is outputted by loadLd function
#' @keywords AVS,Granges
#' @examples
#' ld<-loadLd(file.path(system.file("extdata", "ld_BCa_raggr.csv", package="VSE")), type="raggr")
#' avs<-makeAVS(ld)
#' @importFrom igraph graph.data.frame V clusters
#' @import GenomicRanges
#' @export
makeAVS <- function(x){
  foo <- list()
  ldx <- data.frame(tag=as.character(x$idTag), ld=as.character(x$idLd))
  tags <- unique(ldx$tag)
  for (i in 1:length(tags)){
    foo[[i]] <- unique(as.character(ldx[ldx$tag==tags[i],2]))
  }
  edges <- do.call(
    rbind,
    lapply(foo, function(x){
        if (length(x) > 1) cbind (head(x, -1), tail(x,-1))
        else NULL
    })
  )
  g <- igraph::graph.data.frame(edges, directed=FALSE)
  ldblocks <- split(V(g)$name, clusters(g)$membership)
  avs.glist <- GRangesList(x[as.character(elementMetadata(x)[,1]) %in% as.character(ldblocks[[1]])])
#  avs.ids <- as.character(elementMetadata(avs.glist[[1]])[1,2]);
  for (i in 2:length(ldblocks)){
    gr <- GRangesList(x[as.character(elementMetadata(x)[,1]) %in% as.character(ldblocks[[i]])])
    avs.glist <- c(avs.glist, gr)
    #avs.ids <- c(avs.ids, as.character(elementMetadata(gr[[1]])[1,2]))
  }
  return(avs.glist)
}
