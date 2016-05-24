#'Sample network structure
#'
#'An examples of networks together with network, edge and vertex attributes
#'used primarly for testing. The same networks are stored in objects of class
#'\code{network} and \code{igraph}.
#'
#'Vertices and edges has attribute \code{label}. For vertices these are simply
#'letters from "a" to "o". For edges these are two-letter sequences
#'corresponding to the neightboring vertices, i.e. the label for the edges
#'linking nodes "b" and "c" will be "bc". The order is irrelevant.
#'
#'In the \code{exNetwork} object the \code{label} attribute is also copied to
#'the \code{vertex.names} attribute to facilitate plotting.
#'
#'The \code{exIgraph} object has additional graph attribute \code{layout} so
#'that by default Fruchterman-Reingold placement is used for plotting.
#'
#'@name exNetwork
#'@aliases exNetwork exIgraph exNetwork2 exIgraph2
#'@docType data
#'@format \describe{ \item{exNetwork,exNetwork2}{is of class \code{network}}
#'\item{exIgraph,exIgraph2}{is of class \code{igraph}} } Objects
#'\code{exNetwork} and \code{exIgraph} store directed version of the network.
#'Objects \code{exNetwork2} and \code{exIgraph2} store the undirected version:
#'all direction information from the edges is removed.
#'
#'The network consists of 15 vertices and 11 edges. \itemize{ \item Vertex 1 is
#'an isolate.  \item Vertices 2-6 constitute a star with vertex 2 as a center.
#'\item Vertices 7-8 and 9-10 make two dyads \item Vertcies 11, 12, 13,14 and
#'15 make a stem-and-leaf network. }
#'@keywords datasets
#'@example examples/exNetwork.R
NULL

