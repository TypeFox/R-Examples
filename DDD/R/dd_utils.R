brts2phylo <- function(times,root=FALSE,tip.label=NULL)
{
# This code is taken from the package TESS by Sebastian Hoehna, where the function is called tess.create.phylo
# It takes a set of branching times and adds a random topology.
  times = sort(times)
  n <- as.integer(length(times))+1
  if ( root ) {
    n <- n-1
  }
  nbr <- 2*n - 2

  # create the data types for edges and edge-lengths
  edge <- matrix(NA, nbr, 2)
  edge.length <- numeric(nbr)

  h <- numeric(2*n - 1) # initialized with 0's
  pool <- 1:n
  # VERY VERY IMPORTANT: the root MUST have index n+1 !!!
  nextnode <- 2L*n - 1L
  if ( n > 1) {
    for (i in 1:(n - 1)) {
      # sample two nodes that have no parent yet
      y <- sample(pool, size = 2)
      # compute the edge indices (we just order the edges from 1 to 2n-2)
      ind <- (i - 1)*2 + 1:2
      # set the source node of the new edges (i.e. the new internal node)
      edge[ind, 1] <- nextnode
      # set the destination of the new edges (i.e. the two sampled nodes)
      edge[ind, 2] <- y
      # compute the edge length from the difference between the node heights (child <-> parent)
      edge.length[ind] <- times[i] - h[y]
      # store the node height of this new internal node
      # we cannot use times because then we would get into trouble with the indices and we would need to check for tip nodes ...
      h[nextnode] <- times[i]
      # reset the pool of available nodes to merge
      pool <- c(pool[! pool %in% y], nextnode)
      # increase the node index counter
      nextnode <- nextnode - 1L
    }
  }

  phy <- list(edge = edge, edge.length = edge.length)
  if (is.null(tip.label))
    tip.label <- paste("t", 1:n, sep = "")
  phy$tip.label <- sample(tip.label)
  phy$Nnode <- n - 1L

  if ( root ) {
    phy$root.edge <- times[n] - times[n-1]
    phy$root <- times[n] - times[n-1]
  }

  class(phy) <- "phylo"

  phy <- reorder(phy)
  ## to avoid crossings when converting with as.hclust:
  phy$edge[phy$edge[, 2] <= n, 2] <- 1:n

  return(phy)
}

conv = function(x,y)
{
   lx = length(x)
   ly = length(y)
   lxy = length(x) + length(y)
   x = c(x,rep(0,lxy - lx))
   y = c(y,rep(0,lxy - ly))
   cvxy = rep(0,lxy)
   for(i in 2:lxy)
   {
      cvxy[i] = crossprod(x[(i-1):1],y[1:(i-1)])
   }
   return(cvxy[2:lxy])
}

flavec = function(ddep,la,mu,K,r,lx,kk,n0)
{
   if(ddep == 1 | ddep == 5)
   {
       lavec = pmax(rep(0,lx),la - 1/(r + 1) * (la - mu)/K * ((0:(lx - 1)) + kk))
   }
   if(ddep == 1.3)
   {
       lavec = pmax(rep(0,lx),la * (1 - ((0:(lx - 1)) + kk)/K))
   }
   if(ddep == 2 | ddep == 2.1 | ddep == 2.2)
   {
       x = -(log(la/mu)/log(K + n0))^(ddep != 2.2)
       lavec = pmax(rep(0,lx),la * (((0:(lx - 1)) + kk) + n0)^x)
   }
   if(ddep == 2.3)
   {
       x = K
       lavec = pmax(rep(0,lx),la * (((0:(lx - 1)) + kk) + n0)^x)
   }
   if(ddep == 3 | ddep == 4 | ddep == 4.1 | ddep == 4.2)
   {
       lavec = la * rep(1,lx)
   }
   return(lavec)
}

L2phylo = function(L,dropextinct = T)
# makes a phylogeny out of a matrix with branching times, parent and daughter species, and extinction times
{
   L = L[order(abs(L[,3])),1:4]
   age = L[1,1]
   L[,1] = age - L[,1]
   L[1,1] = -1
   notmin1 = which(L[,4] != -1)
   L[notmin1,4] = age - L[notmin1,4]
   if(dropextinct == T)
   {
      sall = which(L[,4] == -1)
      tend = age
   } else {
      sall = which(L[,4] >= -1)
      tend = (L[,4] == -1) * age + (L[,4] > -1) * L[,4]
   }
   L = L[,-4]
   linlist = cbind(L[sall,],paste("t",abs(L[sall,3]),sep = ""),tend)
   done = 0
   while(done == 0)
   {
      j = which.max(linlist[,1])
      daughter = linlist[j,3]
      parent = linlist[j,2]
      parentj = which(parent == linlist[,3])
      parentinlist = length(parentj)
      if(parentinlist == 1)
      {
         spec1 = paste(linlist[parentj,4],":",as.numeric(linlist[parentj,5]) - as.numeric(linlist[j,1]),sep = "")
         spec2 = paste(linlist[j,4],":",as.numeric(linlist[j,5]) - as.numeric(linlist[j,1]),sep = "")
         linlist[parentj,4] = paste("(",spec1,",",spec2,")",sep = "")
         linlist[parentj,5] = linlist[j,1]
         linlist = linlist[-j,]
      } else {
         #linlist[j,1:3] = L[abs(as.numeric(parent)),1:3]
         linlist[j,1:3] = L[which(L[,3] == as.numeric(parent)),1:3]
      }
      if(is.null(nrow(linlist))) { done = 1 }
   }
   linlist[4] = paste(linlist[4],":",linlist[5],";",sep = "")
   phy = read.tree(text = linlist[4])
   tree = as.phylo(phy)
   return(tree)
}

roundn = function(x, digits = 0)
{
    fac = 10^digits
    return(trunc(fac * x + 0.5)/fac)
}

sample2 = function(x,size,replace = FALSE,prob = NULL)
{
    if(length(x) == 1)
    { 
        x = c(x,x)
        prob = c(prob,prob)
    }
    return(sample(x,size,replace,prob))
}

simplex = function(fun,trparsopt,optimpars,...)
{
  numpar = length(trparsopt)
  reltolx = optimpars[1]
  reltolf = optimpars[2]
  abstolx = optimpars[3]
  maxiter = optimpars[4]

  ## Setting up initial simplex
  v = t(matrix(rep(trparsopt,each = numpar + 1),nrow = numpar + 1))
  for(i in 1:numpar)
  {
      parsoptff = 1.05 * trparsopt[i]/(1 - trparsopt[i])
      trparsoptff = parsoptff/(1 + parsoptff)
      fac = trparsoptff/trparsopt[i]
      if(v[i,i + 1] == 0)
      {
         v[i,i + 1] = 0.00025
      } else {
         v[i,i + 1] = v[i,i + 1] * min(1.05,fac)
      }
  }
  
  fv = rep(0,numpar + 1)
  for(i in 1:(numpar + 1))
  {
     fv[i] = -fun(trparsopt = v[,i], ...)
  }
  
  how = "initial"
  itercount = 1
  string = itercount
  for(i in 1:numpar)
  {
     string = paste(string, v[i,1]/(1 - v[i,1]), sep = " ")
  }
  string = paste(string, -fv[1], how, "\n", sep = " ")
  cat(string)
  flush.console()
  
  tmp = order(fv)
  if(numpar == 1)
  {
     v = matrix(v[tmp],nrow = 1,ncol = 2)
  } else {
     v = v[,tmp]
  }
  fv = fv[tmp]
  
  ## Iterate until stopping criterion is reached
  rh = 1
  ch = 2
  ps = 0.5
  si = 0.5
  
  v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))
  
  while(itercount <= maxiter & ( ( is.nan(max(abs(fv - fv[1]))) | (max(abs(fv - fv[1])) - reltolf * abs(fv[1]) > 0) ) + ( (max(abs(v - v2) - reltolx * abs(v2)) > 0) | (max(abs(v - v2)) - abstolx > 0) ) ) )
  { 
     ## Calculate reflection point
  
     if(numpar == 1)
     {
         xbar = v[1]
     } else {
         xbar = rowSums(v[,1:numpar])/numpar
     }
     xr = (1 + rh) * xbar - rh * v[,numpar + 1]
     fxr = -fun(trparsopt = xr, ...)
   
     if(fxr < fv[1])
     {
         ## Calculate expansion point
         xe = (1 + rh * ch) * xbar - rh * ch * v[,numpar + 1]
         fxe = -fun(trparsopt = xe, ...)
         if(fxe < fxr)
         {
             v[,numpar + 1] = xe
             fv[numpar + 1] = fxe
             how = "expand"
         } else {
             v[,numpar + 1] = xr
             fv[numpar + 1] = fxr
             how = "reflect"
         }
     } else {
         if(fxr < fv[numpar])
         {      
             v[,numpar + 1] = xr
             fv[numpar + 1] = fxr
             how = "reflect"
         } else {
             if(fxr < fv[numpar + 1])
             {
                ## Calculate outside contraction point
                xco = (1 + ps * rh) * xbar - ps * rh * v[,numpar + 1]
                fxco = -fun(trparsopt = xco, ...)
                if(fxco <= fxr)
                {
                   v[,numpar + 1] = xco
                   fv[numpar + 1] = fxco            
                   how = "contract outside"
                } else {
                   how = "shrink"
                }
             } else {
                ## Calculate inside contraction point
                xci = (1 - ps) * xbar + ps * v[,numpar + 1]
                fxci = -fun(trparsopt = xci, ...)
                if(fxci < fv[numpar + 1])
                {  
                   v[,numpar + 1] = xci
                   fv[numpar + 1] = fxci
                   how = "contract inside"
                } else {
                   how = "shrink"
                }
             }
             if(how == "shrink")
             {
                 for(j in 2:(numpar + 1))
                 {
  
                     v[,j] = v[,1] + si * (v[,j] - v[,1])
                     fv[j] = -fun(trparsopt = v[,j], ...)
                 }
             }
         }
     }
     tmp = order(fv)
     if(numpar == 1)
     {
        v = matrix(v[tmp],nrow = 1,ncol = 2)
     } else {
        v = v[,tmp]
     }
     fv = fv[tmp]
     itercount = itercount + 1
     string = itercount;
     for(i in 1:numpar)
     {
         string = paste(string, v[i,1]/(1 - v[i,1]), sep = " ")
     }
     string = paste(string, -fv[1], how, "\n", sep = " ")
     cat(string)
     flush.console()
     v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))
  }
  if(itercount < maxiter)
  {
     cat("Optimization has terminated successfully.","\n")
  } else {
     cat("Maximum number of iterations has been exceeded.","\n")
  }
  out = list(par = v[,1], fvalues = -fv[1], conv = as.numeric(itercount > maxiter))
  invisible(out)
}

optimizer = function(optimmethod = 'simplex',optimpars = c(1E-4,1E-4,1E-6,1000),fun,trparsopt, ...)
{
    if(optimmethod == 'simplex')
    {
        out = simplex(fun = fun,trparsopt = trparsopt,optimpars = optimpars,...)
    }
    if(optimmethod == 'subplex')
    {
        minfun = function(fun,trparsopt,...)
        {           
           return(-fun(trparsopt = trparsopt,...))
        }
        out = subplex::subplex(par = trparsopt,fn = minfun,control = list(abstol = optimpars[3],reltol = optimpars[1],maxit = optimpars[4]),fun = fun,...)
        out = list(par = out$par, fvalues = -out$value, conv = out$convergence)
    }
    return(out)
}