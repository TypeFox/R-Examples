mixer<-function( x, qmin=2, qmax=NULL, method="variational",
                 directed=NULL, nbiter=10, fpnbiter=5, improve=FALSE, verbose=TRUE )
  {


  ## How the graph is coded ?
  if (is.character(x) & (length(x) == 1) ){
    ## spm file
    g <- new("spmgraph",x)
    m.save <- getEdges(g)
    m <- m.save
    NodeName <- g@nodenames
    NbrNodes <- max( m ) 
    if( is.null(directed) ) {
      directed = ! is.symmetric(m)
      if (verbose) {
        cat("Mixer: the edge list has been transformed to ")
        if (directed) 
          cat("a directed adjacency matrix\n")
        else
          cat("an undirected adjacency matrix\n")
      }
    } else  if( !(directed) & !(is.symmetric( m )) ) {
      cat("Mixer: unsymmetric matrix not suitable with directed=FALSE \n")
      return(NULL)
    }
      
    # Only connected nodes are in "m" 
    # The nodes are renumbered
  } else if (dim(x)[1]==dim(x)[2]){

    ## Adjacency matrix
    ##
    if( is.null(directed) ) {
      directed <- (! is.symmetric( x ) )
      if (verbose) {
        cat("Mixer: the adjacency matrix has been transformed in a ")
        if (directed) 
          cat("directed edge list\n")
        else
          cat("undirected edge list\n")
      }
    } else if( !(directed) & !(is.symmetric( x )) ) {
      cat("Mixer: unsymmetric matrix not suitable with directed=FALSE \n")
      return(NULL)
    }
    
    NbrNodes <- dim(x)[1]
    NodeName <- dimnames(x)[1]
    m <- AdjMat2Edges(x, directed=directed, verbose=verbose)
    m.save <- m
  } else if( dim(x)[1] == 2) {
    ## Edge list
        m <- x
    NodeName <- getNodeName( m ) 
    # To avoid index gap
    NbrNodes <- length( NodeName )
    m <- renumberNodes( m )
    m.save <- m

    if( is.null(directed) ) {
      directed <- (! is.symmetric( m ) )
      if( verbose) {
        cat("Mixer: the edge list has been transformed in a ")
        if (directed) 
          cat("directed one\n")
        else
          cat("undirected one\n")
      }
    } else if( !(directed) & !(is.symmetric( m )) ) {
      cat("Mixer: unsymmetric matrix not suitable with directed=FALSE \n")
      return(NULL)
    }
  } else {
    cat("Mixer: not an adjacency matrix or bad edge list\n")
  }


  ## Get the mixnet node order
  # Invalid : readingOrder<-unique(as.numeric(m));
  
  # Get the mapping Mixnet -> initial Graph
  Mixnet2Graph <- getMixnetNumbering( m )
  m <- removeLoops( m )
  m <- renumberNodes( m ) 

  # NbrConnectedNodes : number of connected nodes
  NbrConnectedNodes <- length( Mixnet2Graph )  

  ## prepare the arguments
  undirected<- !( directed) 
  loop<-FALSE
  kmeansnbclass<- 0   # Accelerate the initialization (used to start the HAC
                      # (should be between NbrConnectedNodes and qmax)
  kmeansnbiter<-30    #  
  emeps<-1e-10        # tolerance for em (compare the likelihood)
  fpeps<-1e-4         # tolerance for the fixed point internal loop
                      # (on the taus...)
  nokmeans<-TRUE      # No acceleration via kmeans
  
  silent<-TRUE        # no verbose
  initnbv<-0          # size of the initial network for the online version
  improvenbiter<-3    # number of iteration for the improvment phase


  ## ensure the options compatibility
  if  (method=="classification") {
    classif<-TRUE; stochastique<-FALSE; online<-TRUE}
  else    {
    stochastique<-FALSE; classif<-FALSE; online<-FALSE}


  if (undirected==TRUE){
    symetrize<-TRUE}
  else {
    symetrize<-FALSE} 
    
  ## Ensure number of classes coherence
  if (is.null(qmax)){
    qmax<-qmin
  }

  if( NbrConnectedNodes < qmax ) {
    stop("q-class value greater than the number of nodes.")
  } 
  

  ## compute the size of the returned array from the .C call
  
  nbrClasses <- qmax - qmin + 1
  span       <- qmin:qmax
  nbrICL     <- nbrClasses;         elts <- c(nbrICL)
  nbrAlphas  <- sum(span);          elts <- c(elts, nbrAlphas)
  nbrPis     <- sum(span*span);     elts <- c(elts, nbrPis)
  nbrTaus    <- NbrConnectedNodes*nbrAlphas; elts <- c(elts, nbrTaus)
  nbrValues  <- sum(elts)

  ##Chose the method for the parameter estimation
  if (method=="bayesian"){
    bout <- VariationalBayes(m, qmin, qmax, nbiter, fpnbiter,
                             emeps, fpeps, directed )
  }
  else if (method=="classification"){
    xout <- .C("main_ermgo",
          as.integer(loop),
          as.integer(silent),
          as.integer(initnbv),
          as.integer(improvenbiter),
          as.integer(nbiter),
          as.integer(improve),
          as.integer(classif),
	  as.integer(stochastique),
          as.integer(qmax),
          as.integer(qmin),
          nbrEdges = as.integer(length(m)/2),# size of the given array
          size = as.integer(nbrValues),  # size of the returned array
          lNodes = as.integer(m),       # given array
          res = double(nbrValues))  # returned array
     y <- vector("list", length(span))
   } else {
         xout <- .C("main_ermg",
          as.integer(symetrize),
          as.integer(loop),
          as.integer(undirected),
          as.integer(silent),
          as.integer(improvenbiter),
          as.integer(kmeansnbclass),
          as.integer(kmeansnbiter),
          as.integer(nbiter),
          as.double(emeps),
          as.integer(fpnbiter),
          as.double(fpeps),
          as.integer(improve),
          as.integer(classif),
          as.integer(nokmeans),
          as.integer(qmax),
          as.integer(qmin),
          nbrEdges = as.integer(length(m)/2),# size of the given array
          size = as.integer(nbrValues),  # size of the returned array
          lNodes = as.integer(m),       # given array
          res = double(nbrValues))  # returned array
         
         y <- vector("list", length(span))
       }

  if (method != "bayesian") {
    j <- 1
    cur <- 1
    for (i in span){
      ## format : y[[j]]$name <- dataFormat(x$res[cur:end]);
      ##          cur <- (offset equal to the size of dataFormat)
      y[[j]]$criterion    <- xout$res[cur]; cur <- cur+1
      y[[j]]$alphas <- xout$res[cur:(cur-1+i)]; cur <- cur+i
      y[[j]]$Pis    <- matrix(xout$res[cur:(cur-1+(i*i))], i,i);
                       cur <- cur+(i*i)

      tmp <- matrix(xout$res[cur:(cur-1+(i*NbrConnectedNodes))], i,
                              NbrConnectedNodes,byrow=TRUE);

      # Invalid : y[[j]]$Taus[,readingOrder] <- y[[j]]$Taus
      # replaced by :
      
      y[[j]]$Taus <- matrix( 0, i, NbrNodes ) 
      y[[j]]$Taus[ , Mixnet2Graph[ ]] <- tmp[ ,  ]
      
      # y[[j]]$Taus   <- matrix(xout$res[cur:(cur-1+(i*NbrConnectedNodes))], i,
      #                        NbrConnectedNodes,byrow=TRUE);
      cur <- cur+(i*NbrConnectedNodes)
      j <- j+1
    }
    result<-list(method=method,nnames=NodeName, nnodes=NbrNodes,
                 map=Mixnet2Graph, edges=m.save,qmin=qmin,qmax=qmax,output=y,
                 directed=directed)
  } else {
    j <- 1
    for (i in span){
      tmp <- bout[[j]]$Taus 
      bout[[j]]$Taus <- matrix( 0, i, NbrNodes ) 
      bout[[j]]$Taus[ , Mixnet2Graph[ ]] <- tmp[ ,  ]
      j <- j+1
    }
    
    result<-list(method=method, nnames=NodeName, nnodes=NbrNodes,
                 map=Mixnet2Graph, edges=m.save, qmin=qmin, qmax=qmax, output=bout,
                 directed=directed)
  }
  class(result)<-"mixer"
  return(result)
}
############################################################
# Plot the icl criterion
############################################################

ploticl<-function(x,q,...)
  {
    if (x$method == "bayesian" ){
      title = "Bayesian criterion vs class number"
      y.lab = "Bayesian criterion"
    } else {
      title = "Integrated Classification Likelihood"
      y.lab = "ICL"
    }
    Q<-unlist(lapply(x$output,ICL<-function(x) length(x$alphas)))
    ICL<-unlist(lapply(x$output,ICL<-function(x) x$criterion))
    plot(Q,ICL,xlab="Number of classes",ylab=y.lab,main=title)
    lines(Q,ICL)
    abline(v=q,col="red",lty=2)
  }

############################################################
# Plot the reorganized adjacency matrix
############################################################

plotam<-function(edges,cluster)
  {
    neworder<-order(cluster)
    max(edges)->n
    m<-t(matrix(order(neworder)[as.numeric(edges)],2))
    plot(1, 1, xlim = c(0, n + 1), ylim = c(n + 1, 0), type = "n", axes= FALSE,xlab="classes",ylab="classes",main="Reorganized Adjacency matrix")
    rect(m[,2]-0.5,m[,1]-0.5,m[,2]+0.5,m[,1]+0.5,col=1)
    rect(m[,1]-0.5,m[,2]-0.5,m[,1]+0.5,m[,2]+0.5,col=1)
    table(cluster)->limits # find the class limits
    cumsum(limits)[1:(length(limits)-1)]+0.5->limits
    abline(v=c(0.5,limits,n+0.5),h=c(0.5,limits,n+0.5),col="red")
  }

############################################################
# Plot the Pis matrix and alphas vector using spectral decomposition
############################################################
plotparam<-function(Pis,alphas,q=NULL){
length(alphas)->q
if (q==1) {D<-list(vector=data.frame(1,1)); a<-b<-1} else {
if (q==2) {a<-b<-1} else {a<-2; b<-3}
D<-colSums(Pis)
L<-diag(rep(1,q)) -  diag(D^(-1/2)) %*% Pis %*% diag(D^(-1/2))
eigen(L)->D
}
plot(D$vector[,a],D$vector[,b],cex=1/min(alphas^(1/2))*alphas^(1/2)*3,axes=FALSE,xlab="",ylab="",main="Specral view of the connection matrix",pch=19,col="red")
points(D$vector[,a],D$vector[,b],cex=1/min(alphas^(1/2))*alphas^(1/2)*3)

text(D$vector[,a],D$vector[,b],label=1:q)  
#gplot((Pis>median(Pis))*Pis,vertex.cex=1/min(alphas^(1/2))*alphas^(1/2)*3,edge.lwd=(Pis>median(Pis))*Pis*1/min(median(Pis)),label=1:length(alphas),label.pos=6)
}


############################################################
# Plot the reorganized adjacency matrix
############################################################


mixture<-function(x,alphas,lambdaq){
  fx<-0; for (q in 1:length(alphas)) {
    fx<-fx+alphas[q]*dpois(x,lambda=lambdaq[q])
  }
  return(fx)
}

plotmixture<-function(degrees,Pis,alphas,n, directed=FALSE){
  if( directed )
    colSums(Pis*alphas)*(2*n-2)->lambdaq
  else
    colSums(Pis*alphas)*(n-1)->lambdaq

  # Remove unconnected nodes
  degrees <- degrees[ which( degrees != 0) ]
  min(degrees):max(degrees)->x
  mixture(x,alphas,lambdaq)->y
  histo<-hist(degrees,plot=FALSE)
  plot(histo,ylim=c(0,max(histo$density,y)),freq=FALSE,col=7,main="Degree distribution",)
  lines(x,y,lwd=2,col="blue")
  points(x,y)
  }
  

############################################################
# Plot the estimated degree distribution
############################################################
is.mixer<-function(x){if (class(x)=="mixer") TRUE else FALSE}

plot.mixer<-function(x, q=NULL, frame=1:4, classes=NULL, classes.col=NULL, quantile.val=0.1, ...){  

  # Test x
  if (!is.mixer( x ))
    stop("Not a mixer object")
  x->mixer.res
  
  # Test q
  if( ! is.null(q)) {
    if( ! (q %in% x$qmin:x$qmax)  )
      stop("Bad value of 'q'")
  }

  # Test frame
  if( ! (is.numeric(frame) &  all( frame %in% 1:5) ) )
      stop("Bad frame number")
  
  # Test classes
  if( ! ( is.factor( classes ) | is.null( classes ) ))
      stop("'classes' not factor")
  
  if( ! is.null( classes) & length(classes) != x$nnodes )
      stop("Bad 'classes' length ")
    
  if( ! is.numeric(quantile.val) )
      stop("Bad 'quantile.val' value")

  #
  # Frames
  #
  
  # Remove bad Frames numbers
  index <- which( frame > 5 )
  if( length( index ) != 0 )
    frame <- frame[ - index ]
  frame <- unique( frame )
  
  nb.frame <-  length( frame )

  # Tool large number of frames : remove the last frames
  if( nb.frame > 4 ) {
      frame <- frame[ -5:-nb.frame]
      nb.frame <- 4
  }
  
  n<-dim(mixer.res$x)[1]
 
  if ( is.null(q) ) {
    # find the best number of classes according ICL
    ICL<-unlist( lapply( mixer.res$output,
                         ICL<-function(x) x$criterion))
    which.max(ICL)->i
    q<-length(mixer.res$output[[i]]$alphas)
  }

  index <- q-mixer.res$qmin+1
  apply(mixer.res$output[[index]]$Taus,2,which.max)->cluster
  # Not connected nodes
  cluster[ (colSums( mixer.res$output[[index]]$Taus ) == 0 ) ] <- -1
  
  Pis    <- mixer.res$output[[q-mixer.res$qmin+1]]$Pis
  alphas <- mixer.res$output[[q-mixer.res$qmin+1]]$alphas

  # Frames to display
  nb.col <- 2
  nb.lin <- 2
  if ( nb.frame == 1) {
    nb.col <- 1
    nb.lin <- 1
  } else if ( nb.frame == 2) {
    nb.col <- 2
    nb.lin <- 1
  }
  
  par(mfrow=c(nb.lin, nb.col))

  if( 1 %in% frame ){
    ploticl(mixer.res,q)
  }
  if( 2 %in% frame ) {
    plotam(mixer.res$edges,cluster)
  }
  if( 3 %in% frame ){
    Degrees <- rep( 0, mixer.res$nnodes )
    for ( i in 1:dim(mixer.res$edges)[2] ) {
      # Remove loops
      if( mixer.res$edges[1, i] != mixer.res$edges[2, i] ) {
        node = mixer.res$edges[1, i]
        Degrees[ node ] = Degrees[ node ] + 1
        node = mixer.res$edges[2, i]
        Degrees[ node ] = Degrees[ node ] + 1
      }
    }
    plotmixture( Degrees, Pis, alphas, length( mixer.res$map ), mixer.res$directed  ) 
  }
  if( 4 %in% frame ){
    if( !is.null( classes ) ) {
      pie.coef = table( factor( cluster, levels=1:q ) , classes )

      # Normalization
      for ( i in 1:dim(pie.coef)[1] ) {
        max =  max( pie.coef[i, ] )
        if ( max != 0 )
          pie.coef[i, ] = pie.coef[i, ] / max
      }
    } else {
      pie.coef = NULL
    }
    Gplot( Pis, type="pie.nodes", node.weight=alphas, node.pie.coef=pie.coef,
           quantile.val = quantile.val, colors=classes.col,
           main="Inter/intra class probabilities",
          ... )
  }
  if( 5 %in% frame )
    Gplot( mixer.res$edges, class=cluster, colors=classes.col,
          main="Graph", directed=x$directed,
          ... )
    
  par( mfrow=c(1, 1) )
}


 
############################################################
# Simulation of an affiliation graph
############################################################

class.ind<-function (cl)
{ 
    n <- length(cl)
    cl <- as.factor(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (unclass(cl) - 1)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    x
}


graph.affiliation<-function( n=100,
                             alphaVect=c(1/2,1/2), lambda=0.7, epsilon=0.05,
                             directed=FALSE) {
      # INPUT  n: number of vertex
      #           alphaVect : vecteur of class proportion
      #           lambda: proba of edge given  same classe
      #           epsilon: proba of edge given two different classes
      # OUTPUT x: adjacency matrix
      #              cluster: class vector
      #           
     
      x<-matrix(0,n,n);
      Q<-length(alphaVect);
      NodeToClass <- vector(length=n) 
      rmultinom(1, size=n, prob = alphaVect)->nq;
      Z<-class.ind(rep(1:Q,nq));
      Z<-Z[sample(1:n,n),];
      for (i in 1:n) {
        NodeToClass[i] <- which.max( Z[i,] )
      }
      for (i in 1:n) {
        if ( i != n) {
          for (j in (i+1):n) {
            # if i and j in same class
            if ( NodeToClass[i] ==  NodeToClass[j]) p<-lambda else  p<-epsilon
            if ( (rbinom(1,1,p) )) { x[i,j] <- 1 }
          }
          if ( directed ) {
            if ( i != 1) {
              for (j in 1:(i-1)) {
                if ( NodeToClass[i] ==  NodeToClass[j]) p<-lambda else  p<-epsilon
                if ( (rbinom(1,1,p) )) { x[i,j] <- 1 }
              }
            }
          }
        }
      }
      if ( ! directed ) {
        x <- x + t(x)
      }
      return(list(x=x,cluster=apply(Z,1,which.max)) )   
  }


##############################################################
#  Spectral Clustering using normalized Laplacian
##############################################################
spectralkmeans<-function(x,q=2){
  #INPUT:
  #    x is an adjacency matrix
  #OUTPUT:
  #    An object of class "kmeans" which is a list with components:
  n<-dim(x)[1]
  D<-colSums(x)
  L<-diag(rep(1,n)) -  diag(D^(-1/2))%*% x %*% diag(D^(-1/2))
  eigen(L)->D
  kmeans(as.matrix(D$vectors[,max(1,(n-q)): (n-1)]),q)->res         
}

##############################################################
#  Compute the rand index between two partition
##############################################################
randError<-function(x, y) {
  # function to calculate the adjusted rand statistic
  # x and y are vectors containing the two partitions to be compared
  # first, get crosstabs
  ctab <- table(x,y);

  # now calculate 4 intermediary sums
  cellsum <- sum(ctab*(ctab-1)/2)
  totsum <- sum(ctab)*(sum(ctab)-1)/2

  # use matrix multiplication to get row and column marginal sums
  rows <- ctab %*% rep(1,ncol(ctab))
  rowsum <- sum(rows*(rows-1)/2)
  cols <- rep(1,nrow(ctab)) %*% ctab
  colsum <- sum(cols*(cols-1)/2)
  # now put them together
  adj.rand <- (cellsum - (rowsum*colsum/totsum))/(.5*(rowsum +colsum)-(rowsum*colsum/totsum))
  return (adj.rand);
}

##############################################################
#  transform of an adjacency matrix  into an array of edges  
##############################################################
AdjMat2Edges<-function(x, directed=FALSE, verbose=TRUE ) {
  
  if (dim(x)[1] == dim(x)[2]){
    
    # Adjacency matrix case
    nbrNodes<-dim(x)[1]
    ConnectedNodes <- getConnectedNodes( x )
    # if ( length(ConnectedNodes) > 0 ) {
    #  x <- x[ ConnectedNodes, ConnectedNodes ]
    # } 
    if ( directed ) {
      m <- t( which( (x==1) , arr.ind=TRUE) )
    } else {
      m <- t( which( (x==1) & (upper.tri(x, diag=TRUE)), arr.ind=TRUE) )
    }
  }

  return( m )
}

##############################################################
#  Return the edge list or adjacency matrix without loops
##############################################################
removeLoops<-function(x, adj=FALSE) {
   if (adj){
     ## Adjacency matrix
     diag(x) <- 0
   } else if ( dim(x)[1] == 2) {
     ## Edge list
     ilist <- which( x[1,] != x[2,])
     if( length(ilist) != 0) {
       x <-  as.matrix( x[ , ilist] )
     }
   } else {
     cat("Mixer: removeLoops not an adjacency matrix nor edge list\n")
   }

   return( x )
 }

##############################################################
#  Return the index list of connected nodes 
##############################################################
getConnectedNodes<-function( x ) {

  
  if (dim(x)[1] == dim(x)[2]){
    
    # Adjacency matrix case
    nbrNodes<-dim(x)[1]
    ConnectedNodes <- which( (rowSums( x ) + colSums(x)) != 0)

  }
  else if ( dim(x)[1] == 2) {
    ConnectedNodes = unique( as.vector( x ) )
  } else {
    cat("Mixer: getConnectedNodes not an adjacency matrix nor edge list\n")
  }
  return( ConnectedNodes )
}

##############################################################
# Set the seed of random functions of C/C++ part
##############################################################
setSeed <- function( seed=1 ) {
  #invisible( .C("srand_stdlib",
  #                 as.integer(seed)
  #              ) )
  set.seed(seed)
}

##############################################################
# Declare a generic function
##############################################################
getModel <- function( object, ... )
{
  UseMethod("getModel", object)
}

##############################################################
# Return model parameters.
##############################################################
getModel.mixer <- function( object, ...) {
  
  # Test x
  if ( !is.mixer( object ) )
    stop("Not a mixer object")
  x <- object
  
  # Get optional parameter q
  q <- sub.param("q", NULL	, ...)
  
  # Test q
  if( ! is.null(q)) {
    if( ! (q %in% x$qmin:x$qmax)  )
      stop("Bad value of 'q'")
  }
  if ( is.null(q) ) {
    # find the best number of classes according ICL
    ICL <- unlist( lapply( x$output, ICL<-function(x) x$criterion))
    i   <- which.max(ICL)
    q   <- length(x$output[[i]]$alphas)
  }

  i <- q - x$qmin + 1
  res <- list(q         = q,
              criterion = x$output[[i]]$criterion ,
              alphas    = x$output[[i]]$alphas,
              Pis       = x$output[[i]]$Pis,
              Taus      = x$output[[i]]$Taus
              )
  return( res ) 

}

##############################################################
# Return tke mapping.
# 'x[2, nbedges]' edge list
##############################################################
getMixnetNumbering <- function(x) {

  NodeList <- vector()
  NbEdges <- dim(x)[2]
  n.nodes <- 0
  for ( i in 1:NbEdges ) {

    if(  x[1,i] != x[2,i] ) {

      # Add in the node list if a new one
      if ( !( x[1, i] %in% NodeList ) ) {
        n.nodes <- n.nodes+1
        NodeList[n.nodes] <-  x[1, i]  
      }

      if ( !( x[2,i] %in% NodeList ) ) {
        n.nodes <- n.nodes+1
        NodeList[n.nodes] <-  x[2,i]  
      }
    }
  }
  
  # Return the mapping Mixnet to the initial graph 'x'
  return( NodeList)
}

##############################################################
# Renumber the nodes.
# 'x[2, nbedges]' edge list
##############################################################
renumberNodes <- function( x ) {

  NodeList <- vector()
  NbEdges <- dim(x)[2]
  n.nodes <- 0
  res <- matrix( 0, dim(x)[1],  dim(x)[2])
  for ( i in 1:NbEdges ) {

    # Add in the node list if a new one
    if ( !( x[1, i] %in% NodeList ) ) {
      n.nodes <- n.nodes+1
      NodeList[n.nodes] <-  x[1, i]
       res[1, i] <- n.nodes
    } else {
       res[1, i] <- which( NodeList == x[1, i] )
    }

    if ( !( x[2,i] %in% NodeList ) ) {
      n.nodes <- n.nodes+1
      NodeList[n.nodes] <-  x[2,i]  
      res[2, i] <- n.nodes
    } else {
      res[2, i] <- which( NodeList == x[2, i] )
    }
  }
  
  # Return the mapping Mixnet to the initial graph 'x'
  return( res )

}

##############################################################
# getNodeNames
# 'x[2, nbedges]' edge list
##############################################################
getNodeName <- function( x ) {

  NodeName <- vector()
  NbEdges <- dim(x)[2]
  n.nodes <- 0
  for ( i in 1:NbEdges ) {

    # Add in the node list if a new one
    if ( !( x[1, i] %in% NodeName ) ) {
      n.nodes <- n.nodes+1
      NodeName[n.nodes] <-  x[1, i]
    }

    if ( !( x[2,i] %in% NodeName ) ) {
      n.nodes <- n.nodes+1
      NodeName[n.nodes] <-  x[2,i]  
    }
  }
  
  # Return the mapping Mixnet to the initial graph 'x'
  return( NodeName )

}

