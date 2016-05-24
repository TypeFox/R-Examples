#### Graph  <-->  {from-to}-matrix
#### =============================  in addition to things in 'graph' package

stopifnot(require(graph),
          require(Matrix))
if(packageVersion("Matrix") < "1.1.1")
    stop("need Matrix >= 1.1.1 (with improved graph -> matrix)")


##' @title The "edgeList" of a graph object as a {from -> to} 2-column matrix
##' @param g a 'graph' object
##' @return "ftM", the 2-column matrix (i, j) of integer indices: one row <=> one edge
##' @seealso \code{\link{ftM2graphNEL}} etc (but I don't find a graph2ftM())
##' @author Martin Maechler
graph2ftM <- function(g) {
    A <- graph2T(g) # <- from "Matrix"
    if(edgemode(g) == "undirected") ## symm.matrix: keep only upper triangle
        A <- triu(A)
    cbind(i = A@i + 1L, j = A@j + 1L)
}

##' From of a graph, given by its edgeList ({from -> to} matrix), 
##' compute the edgeList of the "permuted graph", i.e., where the nodes are permuted
##'
##' @title Update the edgeList of graph when the nodes are permuted
##' @param ftM a {from -> to}-matrix, i.e. 2-column, with rows == (i, j) for i--j edges
##' @param perm a permutation vector, i.e. a permutation of 1:p
##' @param edgemode "directed" or "undirected" must be specified
##' @return a new {from -> to} matrix (of the same dimension)
##' @author Martin Maechler, Dec.12, 2013
perm.ftM <- function(ftM, perm, edgemode = "undirected") {
    stopifnot(is.matrix(ftM), (d <- dim(ftM))[2] == 2) # 1 <= perm <= max(ftM)
    ne <- d[1L] # nrow(ftM) = #{edges}
    if(ne == 0) return(ftM)
    ## else ne >= 1 :
    p <- length(perm)# = #{nodes}
    ## Conceptually: ip <- invPerm(perm); m <- ftM; m[] <- ip[ftM]
    ## ==> m = (i,j) ftMatrix and now "reorder" m[.,.]:
    ## (1) [only if undirected]: for each row swap 1st and 2nd columns such that i <= j;
    ## (2) sort rows of m[] wrt j, i.e., first increasing j (and then increasing i)
    switch(edgemode,
	   "directed" = {
	       m <- ftM; m[] <- invPerm(perm)[ftM]
	       ## (2) sort:
	       m[order(m[,2], m[,1]) , , drop=FALSE]
	   },
	   "undirected" = {
	       ## more directly:  m[,1] = first half of ip[ftM], and m[,2] = 2nd half:
	       ij <- invPerm(perm, zero.res=TRUE)[ftM]
	       graph2ftM(T2graph(new("nsTMatrix", Dim=c(p,p),
				     i = ij[1L:ne], j = ij[(ne+1L):(ne+ne)]),
				 edgemode="undirected"))
	   },
	   stop("invalid 'edgemode': ", edgemode))
	   
}
