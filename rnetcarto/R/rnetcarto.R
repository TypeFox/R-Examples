#' @useDynLib rnetcarto, netcarto_binding
#' @importFrom stats runif
NULL
#' Compute modularity and modularity roles for graphs using simulated
#' annealing
#'
#' @docType package
#' @name rnetcarto

#' @title Computes modularity and modularity roles from a network.
#'
#' @param web network either as a square adjacency matrix or a list
#' describing E interactions a->b: the first (resp. second) element is
#' the vector of the labels of a (resp. b), the third (optional) is
#' the vector of interaction weights.
#' @param seed Seed for the random number generator: Must be a
#' positive integer.
#' @param iterfac At each temperature of the simulated annealing
#' (SA), the program performs fN^2 individual-node updates (involving the
#' movement of a single node from one module to another) and fN
#' collective updates (involving the merging of two modules and the split
#' of a module). The number "f" is the iteration factor.
#' @param bipartite If True use the bipartite definition of modularity.
#' @param coolingfac Temperature cooling factor.
#' @return A list. The first element is a dataframe with the name,
#' module, z-score, and participation coefficient for each row of the
#' input matrix. The second element is the modularity of this
#' partition.
#' @examples
#' # Generate a simple random network
#' a = matrix(as.integer(runif(100)<.3), ncol=10) 
#' a[lower.tri(a)] = 0
#' # Find an optimal partition for modularity using netcarto.
#' netcarto(a)
#' 
#' @export
netcarto <- function(web,
                     seed=as.integer(floor(runif(1, 1,100000001))),
                     iterfac=1.0,
                     coolingfac=0.995,
                     bipartite=FALSE)
{  
    # Read matrix of list input..
    if(is.matrix(web)){
        # Sanity checks...
        if(!bipartite){
            if(ncol(web) != nrow(web)){
                stop("Input matrix must be square for non bipartite networks.")
            }
            
            if(!isSymmetric(web)){
                if (sum(web[lower.tri(web)]!=0,diag=FALSE)!=0 &&  sum(web[upper.tri(web,diag=FALSE)]!=0)!=0){
                    warning("Input matrix should be symmetric or triangular for non bipartite networks. \n (max of web(ij) et web(ji) was taken).\n")
                }
                web[upper.tri(web)] = pmax(t(web)[upper.tri(web)],web[upper.tri(web)])
            }
            
            if (any(rownames(web) != colnames(web))){
                warning("Columns and row names are not matching, are you sure this is an adjacency matrix of a non bipartite network ?")
            }
            if (is.null(rownames(web))){
                rownames(web) = 1:nrow(web)
            }
            
            # Removing the upper part of the matrix.
            web[lower.tri(web)] <- 0
            # Removing empty (lines and colums).
            mask = colSums(web==0)+rowSums(web==0)!=2*nrow(web)
            web = web[mask,mask]
        } else{
            # Removing (empty columns) and (empty lines).
            web = web[rowSums(web==0)!=ncol(web), colSums(web==0)!=nrow(web)]
        }

        # Get non zero positions.
        non_zero <- which(!web == 0)
        
        # Parameters
        nodes1 = row(web)[non_zero]-1L
        nodes2 = col(web)[non_zero]-1L
        weights = web[cbind(row(web)[non_zero], col(web)[non_zero])]
        names = rownames(web)

    } else if (is.list(web)){

        E = length(web[[1]]) # Number of edges

        if (length(web[[1]]) != length(web[[2]])){
            stop("Bad labels number: all elements of the input list should have the same length.")}

        # Read the weight if they are supplied
        if (length(web)==3){
            weights = web[[3]]
            if (length(web[[3]])!=length(web[[1]])){
                stop("Bad weight number: all elements of the input list should have the same length.")
            }
        } else if (length(web)==2){
            weights = numeric(E) + 1
        } else{
            stop("Input edge list should be of length 2 (unweighted edges) or 3 (weighted edges)");
        }

        web[[1]] = as.character(web[[1]])
        web[[2]] = as.character(web[[2]])

        # Convert the species names to integer
        if(bipartite==FALSE){
          fct = factor(c(web[[1]], web[[2]]))
          idx = as.integer(fct) - 1L
          names = levels(fct)
          nodes1 = idx[1:E]
          nodes2 = idx[(E+1):(2*E)]
        } else {
          fct_par1= factor(web[[1]])
          fct_par2= factor(web[[2]])
          nodes1 = as.integer(fct_par1) - 1L
          nodes2 = as.integer(fct_par2) - 1L
          names = levels(fct_par1)
        }
    }
    N = length(names)
    E = length(nodes1)
    roles = 1
    clustering = 1
    diagonal_term = ifelse(bipartite, 0,1)

    if (N<=2 || E<=1){
        stop("Trivial graph (less than 2 nodes or 1 edge).")
    }
    
    # Call rgraphlib
    ans <- .Call("netcarto_binding",
                 as.integer(nodes1),
                 as.integer(nodes2),
                 as.numeric(weights),
                 as.integer(N),
                 as.integer(bipartite),
                 as.integer(clustering),
                 as.integer(roles),
                 as.integer(diagonal_term),
                 as.numeric(coolingfac),
                 as.integer(seed),
                 as.numeric(iterfac),
                 PACKAGE="rnetcarto")

    # Build the dataframe
    df = data.frame(names, ans[[1]], ans[[2]], ans[[3]])
    names(df) <- c("name","module","connectivity","participation")
    df = df[with(df, order(module,connectivity,participation,name)), ]

    # Assign Roles
    zlimits = c(2.5, 2.5, 2.5, 2.5, Inf, Inf, Inf)
    plimits = c(0.050,0.620,0.800,Inf,.03,0.75,Inf)
    levels = c(0,1,2,3,4,5,6)
    df$role = -1
    labels = c("Ultra peripheral","Peripheral","Connector","Kinless","Peripheral Hub","Connector Hub","Kinless Hub")
    for (i in length(levels):1){
        df$role[df$connectivity<zlimits[i]  & df$participation<plimits[i]] = levels[i]
    }
    df$role = factor(df$role,levels=levels,labels,ordered=TRUE)
    return(list(df,ans[[4]]))
}
