hc <- function( data, node.sizes, scoring.func = 0, cpc, cont.nodes = c(), ess = 1, tabu.tenure = 100, 
                init.net = NULL )
{
  n.nodes <- ncol(data)
  n.cases <- nrow(data)
  
  
  # just to be sure
  storage.mode( node.sizes ) <- "integer"
  
  # quantize data of continuous nodes 
  levels <- rep( 0, n.nodes )
  levels[cont.nodes] <- node.sizes[cont.nodes]
  
  # data <- quantize.with.na.matrix( data, levels )
  data <- quantize.matrix( data, levels )
  
  # start with init.net if not NULL, otherwise with empty matrix
  if( !is.null(init.net) )
  {
    curr.g <- init.net
    # just to be sure 
    storage.mode( curr.g ) <- "integer" 
    # add init edges out of cpc
    cpc <- cpc | init.net | t(init.net)
  }
  else
    curr.g <- matrix(0L,n.nodes,n.nodes) # integers!  
  
  curr.score.nodes <- array(0,n.nodes)
  for( i in 1L:n.nodes )
    curr.score.nodes[i] <- .Call( "score_node", data, node.sizes, i-1L, which(curr.g[,i]!=0)-1L, scoring.func,
                                  ess, PACKAGE = "bnstruct" )

  # global best solution
  global.best.g <- curr.g
  global.best.score <- sum(curr.score.nodes)
  
  # tabu list
  tabu <- array(0L, c(n.nodes,n.nodes,tabu.tenure))
  tabu.pt <- 1
  
  # worsening moves
  wm.count <- 0
  wm.max <- 15
  while( wm.count < 15 )
  {
    next.score.diff <- rep(-Inf,n.nodes)
    next.pert <- rep(-1L,n.nodes)
    # print(tabu[,,1:10])
    # try all possible perturbations
    for( node in 1L:n.nodes )
    {
      for( par in 1L:n.nodes )
      {
        if( cpc[par,node] == 1L )
        {
          next.g <- curr.g;
          s.diff <- -Inf
        
          if( curr.g[par,node] == 1L ) # edge removal
          {
            next.g[par,node] = 0L;
            if( .Call("is_acyclic", next.g, PACKAGE = "bnstruct") &!.Call("in_tabu", next.g, tabu, PACKAGE = "bnstruct"))
            {
#               cat(node.sizes,'\t',node-1L,'\t',which(next.g[,node]!=0)-1L,'\n')
              s.diff <- .Call( "score_node", data, node.sizes, node-1L, which(next.g[,node]!=0)-1L, 
                               scoring.func, ess, PACKAGE = "bnstruct" ) - curr.score.nodes[node];
            }
          }
          else if( curr.g[node,par] == 0L ) # edge addition
          {
            next.g[par,node] = 1L;
            # print(c(node,par,.Call("is_acyclic", next.g, PACKAGE = "bnstruct"),!.Call("in_tabu", next.g, tabu, PACKAGE = "bnstruct")))
            if( .Call("is_acyclic", next.g, PACKAGE = "bnstruct") & !.Call("in_tabu", next.g, tabu, PACKAGE = "bnstruct"))
            {
              # print("here\n");
              
#               cat(node.sizes,'\t',node-1L,'\t',which(next.g[,node]!=0)-1L,'\n')
              s.diff <- .Call( "score_node", data, node.sizes, node-1L, which(next.g[,node]!=0)-1L,
                               scoring.func, ess, PACKAGE = "bnstruct" ) - curr.score.nodes[node];
            }
          }
          else # edge reversal
          {
            next.g[par,node] = 1L;
            next.g[node,par] = 0L;
            if( .Call("is_acyclic", next.g, PACKAGE = "bnstruct") & !.Call("in_tabu", next.g, tabu, PACKAGE = "bnstruct"))
            {
#               cat(node.sizes,'\t',node-1L,'\t',which(next.g[,node]!=0)-1L,'\t',par-1L,'\t',which(next.g[,par]!=0)-1L,'\n')
              s.diff <- .Call( "score_node", data, node.sizes, node-1L, which(next.g[,node]!=0)-1L,
                               scoring.func, ess, PACKAGE = "bnstruct" ) + 
                        .Call( "score_node", data, node.sizes, par-1L, which(next.g[,par]!=0)-1L,
                               scoring.func, ess, PACKAGE = "bnstruct" ) -
                        ( curr.score.nodes[node] + curr.score.nodes[par] )
            }
          }
          # test for local improvement
          if( s.diff > next.score.diff[node] )
          {
            next.score.diff[node] <- s.diff
            next.pert[node] <- par
          }
        }
      }
    }
    
    best.node <- which.max(next.score.diff)
    
    # no possible improvements given the tabu list
    if( next.score.diff[best.node] == -Inf )
      break
    
    # update current graph and scores
    curr.g[next.pert[best.node],best.node] = 1L - curr.g[next.pert[best.node],best.node] # flip
    if( curr.g[best.node,next.pert[best.node]] == 1L) # need to reverse
    {
      curr.g[best.node,next.pert[best.node]] = 0L;
      # cat(node.sizes,'\t',best.node-1L,'\t',which(curr.g[,best.node]!=0)-1L,'\t',next.pert[best.node]-1L,'\t',which(curr.g[,next.pert[best.node]]!=0)-1L,'\n')
      curr.score.nodes[best.node] <- .Call( "score_node", data, node.sizes, best.node-1L, 
                                            which(curr.g[,best.node]!=0)-1L, scoring.func, ess, PACKAGE = "bnstruct" )
      curr.score.nodes[next.pert[best.node]] <- .Call( "score_node", data, node.sizes, next.pert[best.node]-1L, 
                                            which(curr.g[,next.pert[best.node]]!=0)-1L, scoring.func, ess, PACKAGE = "bnstruct" )
    }
    else
      curr.score.nodes[best.node] = curr.score.nodes[best.node] + next.score.diff[best.node]
    
    # print(c(tabu.pt,best.node,next.pert[best.node],sum(curr.score.nodes)))
    
    if( global.best.score < sum(curr.score.nodes) ) # check for global best
    {
      wm.count <- 0
      global.best.g <- curr.g
      global.best.score <- sum(curr.score.nodes)
    }
    else
      wm.count <- wm.count + 1
    
    # update tabu list
    tabu[,,tabu.pt] <- curr.g
    tabu.pt <- (tabu.pt)%%tabu.tenure + 1
    # print(curr.g)
  }
  
  return(global.best.g)
}


mmpc <- function( data, node.sizes, cont.nodes = NULL, chi.th = 0.05,
						layering = NULL, layer.struct = NULL)
						
{
  n.nodes <- ncol(data)
  n.cases <- nrow(data)
  
  # just to be sure
  storage.mode( node.sizes ) <- "integer"
  
  # quantize data of continuous nodes 
  levels <- rep( 0, n.nodes )
  levels[cont.nodes] <- node.sizes[cont.nodes]
  
  # data <- quantize.with.na.matrix( data, levels )
  data <- quantize.matrix( data, levels )
    
  # default values for layering
  if( is.null(layering) )
  {
    layering <- rep.int( 1, n.nodes )
    if( !is.null(layer.struct) )
      stop( "Argument layer.struct without layering\n" )
  }
  n.layers <- length(unique(layering))
  if( is.null(layer.struct) )
  {
    layer.struct <- matrix(1,n.layers,n.layers)
    if(n.layers > 1)
    {
      layer.struct[lower.tri(layer.struct)] <- 0
      layer.struct[1,1] <- 0 # default: no edges between nodes at level 1
    }
  }
  # print(n.layers)
  # print(layering)
  # print(layer.struct)
  
  # apply layering
  layer.mat <- matrix(1, n.nodes, n.nodes) 
  for( i in 1:n.layers )
    for( j in 1:n.layers )
      layer.mat[ layering==i, layering==j ] <- layer.struct[i,j]
  diag(layer.mat) <- 0
  
  cpc.mat <- layer.mat | t(layer.mat) # constrain the cpc search
  allowed <- cpc.mat
  
  # print(cpc.mat)
  
  storage.mode(data) <- "integer"
  
  # forward addition of nodes
  for( i in 1:n.nodes )
  {
    cpc.mat[i,] <- mmpc.fwd( data, node.sizes, allowed, i, chi.th )
    # cat("cpcMat ",i,": ",cpc.mat[i,],"\n")
    allowed[,i] <- allowed[,i] & t(cpc.mat[i,])
  }
  
  # print(cpc.mat)
  
  # backwards removal of nodes
  for( i in 1:n.nodes )
  {
    cpc.mat[i,] <- mmpc.bwd( data, node.sizes, cpc.mat[i,], i, chi.th )
    # cat("cpcMat ",i,": ",cpc.mat[i,],"\n")
    cpc.mat[,i] <- cpc.mat[,i] & t(cpc.mat[i,])
  }
  
  # print(cpc.mat)
  
  # symmetry enformcement
  cpc.mat = cpc.mat * t(cpc.mat)
  
  # further filter with layering
  cpc.mat = cpc.mat * layer.mat
  
  # print(cpc.mat)
  
  return( cpc.mat )
}

mmpc.fwd <- function( data, node.sizes, allowed, x, chi.th )
{
  # cat("\n",x,":\n");
  n.nodes <- ncol(data)
  n.cases <- nrow(data)

  # test without conditioning
  # print(chi.th)
  minAssoc <- rep(0,n.nodes)
  for( y in 1:n.nodes )
    if( allowed[x,y] )
      minAssoc[y] <- g2( data, node.sizes, x, y, chi.th )
  allowed[x,minAssoc==0] <- 0 # remove already independent nodes
  
  m <- max( minAssoc )
  if( m == 0 )
    return(rep(0,n.nodes))
  
  m.ind <- which.max( minAssoc )
  cpc <- m.ind
  # cat(cpc,",\t",which( !allowed[x,] ),"\n")
  allowed[x,m.ind] <- 0
  minAssoc[m.ind] <- 0
  
  while( sum(allowed[x,]) > 0 )
  {
    # try to add one node to the cpc
    for( y in 1:n.nodes )
    {
      if( allowed[x,y] )
      {
        # condition on all possible combinations of cpc elements, 
        # from smaller to larger
        n <- as.integer( length(cpc) )
        s <- sort( node.sizes[cpc], index.return=T )$ix # useful for early stopping
        for( zsize in seq_len( n ) )
        {
          # test for possible early stopping, no indep test can be performed
          if( prod(node.sizes[c(x,y,cpc[s[1:zsize]])]) > n.cases / 5 )
            break
          
			    comb <- as.integer(1:zsize)
			    while( comb[1]!=0 ) # check for the end of combinations
			    {
            assoc = g2( data, node.sizes, x, y, chi.th, cpc[comb] )
            if( assoc < minAssoc[y] )
              minAssoc[y] <- assoc
            if( assoc == 0 ) # do not try with y anymore
            {
              allowed[x,y] = 0
              break
            }
				    comb <- .Call( "next_comb", comb, n, PACKAGE = "bnstruct" )
          }
          if( allowed[x,y] == 0 ) # came out from the break
            break
        }
      }
    }
    # cat(minAssoc,"\n")
    m <- max( minAssoc )
    if( m == 0 )
      break
    
    m.ind <- which.max( minAssoc )
    minAssoc[m.ind] <- 0
    cpc <- union( cpc, m.ind )
    # cat(cpc,",\t",which( !allowed[x,] ),"\n")
    allowed[x,m.ind] <- 0
  }
  
  cpc.vec <- rep(0,n.nodes)
  cpc.vec[cpc] <- 1
  return(cpc.vec)
}

mmpc.bwd <- function( data, node.sizes, cpc.vec, x, chi.th )
{
#  cat("\n",x,":\n");
#  cat(which(cpc.vec!=0),"\n");
  # not worth if cpc has less than 2 elmts
  if( sum(cpc.vec) < 2 )
    return( cpc.vec )
  
  n.nodes <- ncol(data)
  n.cases <- nrow(data)
  
  for( y in 1:n.nodes )
  {  
    if( cpc.vec[y]!=0 )
    {
      # condition on all possible combinations of cpc elements, 
      # from smaller to larger
      cpc <- setdiff( which( cpc.vec > 0 ), y )
		  n <- as.integer( length(cpc) )
      s <- sort( node.sizes[cpc], index.return=T )$ix # useful for early stopping
      for( zsize in seq_len( n )  )
      {
        # test for possible early stopping, no indep test can be performed
        if( prod(node.sizes[c(x,y,cpc[s[1:zsize]])]) > n.cases / 5)
          break
        
        comb <- as.integer(1:zsize)
        while( comb[1]!=0 ) # check for the end of combinations
        {  
          assoc = g2( data, node.sizes, x, y, chi.th, cpc[comb] )
          if( assoc == 0 ) # do not try with y anymore
          {
            cpc.vec[y] = 0
#            cat(which(cpc.vec!=0),"\n");
            break
          }
          comb <- .Call("next_comb", comb, n, PACKAGE = "bnstruct" )
        }
        if( cpc.vec[y] == 0 ) # came out from the break
          break
      }
    }
  }
  
  return( cpc.vec )
}
  
g2 <- function( data, sizes, x, y, chi.th = 0.05, z=c() )
{
  # if less than 5 counts per cell on average, cannot exclude dependence
  if( dim(data)[1]/prod(sizes[c(x,y,z)]) < 5 )
    return( chi.th )

  # tab <- compute.counts( data[,c(x,y,z)], sizes[c(x,y,z)] )
#   tab <- .Call("compute_counts", data = data[,c(x,y,z)], node_sizes = sizes[c(x,y,z)], 
#                PACKAGE = "bnstruct" )
#   
#   dim(tab) <- c(sizes[x], sizes[y], prod(sizes[z]))
#   sz <- array( rep(apply(tab, 3, sum), each=sizes[x]*sizes[y]), dim(tab) )
#   syz <- array( rep(apply(tab, c(2,3), sum), each=sizes[x]), dim(tab) )
#   sxz <- aperm( array(rep(apply(tab, c(1,3), sum), each=sizes[y]), 
#                       c(dim(tab)[2],dim(tab)[1],dim(tab)[3])), c(2,1,3) )
#   
#   s <- tab * log((tab*sz)/(syz*sxz)) # there can be NaNs for cells with count 0
#   s[is.na(s)] <- 0
#   
#   print(tab)
#   print(sz)
#   print(syz)
#   print(sxz)
#   print(s)
#   
#   df <- sum(tab!=0) + sum(sz!=0)/(sizes[x]*sizes[y]) - 
#     (sum(sxz!=0)/sizes[y] + sum(syz!=0)/sizes[x])
# 
#   # sanity check
#   if( df < 0 )
#   {
#     cat("Negative df\n")
#     return( chi.th )
#   }
#   
#   print(c(2*sum(s),df))
#   
#   return( max(pchisq(2*sum(s),df) - 1+chi.th, 0) )
  
  stat <- .Call("g2_stat", data = data[,c(x,y,z)], node_sizes = sizes[c(x,y,z)], 
                 PACKAGE = "bnstruct")
  
  # print(stat)
  
  # sanity check
  if( stat[2] < 0 )
  {
    # cat("Negative df\n")
    return( chi.th )
  }

  return( max(pchisq(stat[1],stat[2]) - 1+chi.th, 0) )
}
