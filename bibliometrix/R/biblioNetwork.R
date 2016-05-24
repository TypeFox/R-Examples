#' Creating Bibliographic networks
#'
#' \code{biblioNetwork} creates different bibliographic networks from a bibliographic data frame.
#'
#' The function \code{\link{biblioNetwork}} can create a collection of bibliographic networks following the approach proposed by Batagely and Cerinsek (2013).\cr\cr
#' Typical networks output of \code{biblioNetwork} are:\cr\cr
#' #### Collaboration Networks ############\cr
#' -- Authors collaboration (analysis = "collaboration", network = "authors")\cr
#' -- Country collabortion (analysis = "collaboration", network = "countries")\cr\cr
#' #### Co-citation Networks ##############\cr
#' -- Authors co-citation (analysis = "co-citation", network = "authors")\cr
#' -- Reference co-citation (analysis = "co-citation", network = "references")\cr
#' -- Source co-citation (analysis = "co-citation", network = "sources")\cr\cr
#' #### Coupling Networks ################\cr
#' -- Authors coupling (analysis = "coupling", network = "authors")\cr
#' -- Source coupling (analysis = "coupling", network = "sources")\cr
#' -- Keyword coupling (analysis = "coupling", network = "keywords")\cr
#' -- Author-Keyword coupling (analysis = "coupling", network = "author_keywords")\cr
#' -- Country coupling (analysis = "coupling", network = "countries")\cr\cr
#'
#' @param M is a bibliographic data frame obtained by the converting function
#'   \code{\link{convert2df}}. It is a data matrix with cases corresponding to
#'   manuscripts and variables to Field Tag in the original SCOPUS and Thomson Reuters' ISI Web of Knowledge file.
#' @param analysis is a character object. It indicates the type of analysis have to be performed.
#'   \code{analysis} argument can be \code{"collaboration"}, \code{"coupling"} or \code{"co-citation"}.
#'   Default is \code{analysis = "coupling"}.
#' @param network is a character object. It indicates the network typology. The \code{network} aurgument can be
#' \code{"authors"}, \code{"references"}, \code{"sources"}, \code{"countries"},\code{"keywords"} or \code{"author_keywords"}.
#' Default is \code{network = "authors"}.
#' @param sep is the field separator character. This character separates strings in each column of the data frame. The default is \code{sep = ";"}.
#' @return It is a squared network matrix. It is an object of class \code{dgMatrix} of the package \code{\link{Matrix}}.
#' @examples
#' # EXAMPLE 1: Authors collaboration network
#'
#' library(igraph)
#' data(scientometrics)
#'
#' NetMatrix <- biblioNetwork(scientometrics, analysis = "collaboration", 
#' network = "authors", sep = ";")
#' netDegree <- 2
#' diag <- Matrix::diag 
#' NetMatrix <- NetMatrix[diag(NetMatrix) >= netDegree,diag(NetMatrix) >= netDegree]
#' diag(NetMatrix) <- 0
#'
#' bsk.network <- graph.adjacency(NetMatrix,mode = "undirected")
#' plot(bsk.network,layout = layout.fruchterman.reingold, vertex.label.dist = 0.5,
#' vertex.frame.color = 'blue', vertex.label.color = 'black',
#' vertex.label.font = 1, vertex.label = V(bsk.network)$name, vertex.label.cex = 0.7)
#'
#'
#' # EXAMPLE 2: Co-citation network
#'
#' library(igraph)
#' data(scientometrics)
#'
#' NetMatrix <- biblioNetwork(scientometrics, analysis = "co-citation", 
#' network = "references", sep = ";")
#' netDegree=10
#' diag <- Matrix::diag
#' NetMatrix <- NetMatrix[diag(NetMatrix) >= netDegree,diag(NetMatrix) >= netDegree]
#' diag(NetMatrix) <- 0
#'
#' bsk.network <- graph.adjacency(NetMatrix,mode = "undirected")
#' plot(bsk.network,layout = layout.fruchterman.reingold, vertex.label.dist = 0.5,
#' vertex.frame.color = 'blue', vertex.label.color = 'black',
#' vertex.label.font = 1, vertex.label = V(bsk.network)$name, vertex.label.cex = 0.7)
#'
#' @seealso \code{\link{convert2df}} to import and convert a SCOPUS and Thomson 
#'   Reuters' ISI Web of Knowledge export file in a data frame.
#' @seealso \code{\link{cocMatrix}} to compute a co-occurrence matrix.
#' @seealso \code{\link{biblioAnalysis}} to perform a bibliometric analysis.
#'

biblioNetwork <- function(M, analysis = "coupling", network = "authors", sep = ";"){
  
  crossprod <- Matrix::crossprod
  NetMatrix=NA
  if (analysis=="coupling"){
  switch(network,
      authors={
      WA=cocMatrix(M, Field="AU", type = "sparse", sep)
      WCR=cocMatrix(M, Field="CR", type = "sparse", sep)
      CRA = crossprod(WCR, WA)
      NetMatrix = crossprod(CRA, CRA)
      },
    keywords={
      WK=cocMatrix(M, Field="ID", type = "sparse", sep)
      WCR=cocMatrix(M, Field="CR", type = "sparse", sep)
      CRK = crossprod(WCR, WK)
      NetMatrix = crossprod(CRK, CRK)
      },
    author_keywords={
      WK=cocMatrix(M, Field="DE", type = "sparse", sep)
      WCR=cocMatrix(M, Field="CR", type = "sparse", sep)
      CRK = crossprod(WCR, WK)
      NetMatrix = crossprod(CRK, CRK)
      },
    sources={
      WSO=cocMatrix(M, Field="SO", type = "sparse", sep)
      WCR=cocMatrix(M, Field="CR", type = "sparse", sep)
      CRSO = crossprod(WCR, WSO)
      NetMatrix = crossprod(CRSO, CRSO)
    },
    countries={
      WCO=cocMatrix(M, Field="AU_CO", type = "sparse", sep)
      WCR=cocMatrix(M, Field="CR", type = "sparse", sep)
      CRCO = crossprod(WCR, WCO)
      NetMatrix = crossprod(CRCO, CRCO)
    }
  )}
  if (analysis=="co-citation"){
    switch(network,
           authors={
             WA=cocMatrix(M, Field="CR_AU", type = "sparse", sep)
             NetMatrix = crossprod(WA, WA)
             },
           references={
             WCR=cocMatrix(M, Field="CR", type = "sparse", sep)
             NetMatrix = crossprod(WCR, WCR)
           },
           sources={
             WSO=cocMatrix(M, Field="CR_SO", type = "sparse", sep)
             NetMatrix = crossprod(WSO, WSO)
           }
    )}
  if (analysis=="collaboration"){
    switch(network,
           authors={
             WA=cocMatrix(M, Field="AU", type = "sparse", sep)
             NetMatrix = crossprod(WA, WA)
             },
           countries={
             WCO=cocMatrix(M, Field="AU_CO", type = "sparse", sep)
             NetMatrix = crossprod(WCO, WCO)
           })
    }
return(NetMatrix)
}
