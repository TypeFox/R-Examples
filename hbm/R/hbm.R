mcl <- function(m, infl, iter=1000, remove.self.loops = FALSE,
                prune = FALSE, thresh = 1e-06, pruning.prob = 1e-06, 
                use.sparse = NULL, verbose = FALSE)
{
  if (nrow(m) != ncol(m))
  {
    stop("input matrix must be a square matrix")
  }
  
  if (remove.self.loops)
  {
    diag(m) = 0
  }
  n = nrow(m)^2
  m <- sweep(m,2,colSums(m),`/`)
  m[is.na(m)] <- 0  
  if (prune)
  {
    m[m < pruning.prob] = 0   
  }
  force.sparse = FALSE
  if (length(use.sparse) == 0)
  {
    force.sparse = ((sum(m == 0)/n) >= 0.5)
    use.sparse = force.sparse
  }
  if (use.sparse || force.sparse)# actually made sparse if half the entries are 0 
  {
    {
      m = Matrix(m)
      if (verbose)
      {
        print("sparse matrices will be used")
      }
    }
  }

  m0 <- m
  # make sure not circular
  m <- m %*%m 
  m <- m^infl
  m <- sweep(m,2,colSums(m),`/`)
  m[is.na(m)] <- 0   
  i = 1
  if (sum(m0-m) != 0)
  {
    for (i in 2:iter){
      
      m <- m %*%m # expansion 
      m <- m^infl
      m <- sweep(m,2,colSums(m),`/`)
      m[is.na(m)] <- 0 
      if ((sum( m > 0 & m < 1) == 0) || (sqrt((sum((m-m0)^2)/n)) < thresh))
      {
        break 
      }
      if (prune)
      {
        m[m < pruning.prob] <- 0
      }
      m0 = m
      
    }
  }
  if (verbose)
  {
    print(paste("mcl converged after",i, "iterations"))
  }
  
  if (class(matrix) != "matrix")
  {
    m = as.matrix(m)
  }
  
  nrow <- nrow(m)
  ncol <- ncol(m)
  clusters <- vector(length = nrow, mode = "numeric")
  csums = colSums(m)
  lonely = which(csums == 0)
  clustered = which(csums > 0)
  clusters[lonely] = lonely
  attractors = sort(which(rowSums(m) > 0))
  j = 1
  lcc = length(clustered)
  unlabelled = lcc
  while (unlabelled > 0)
  {
    i = attractors[j]
    if (clusters[i] == 0)
    {
      attracts <- which(m[i,1:ncol] > 0)
      clusters[attracts] <- i
    }
    unlabelled = sum(clusters[clustered] == 0)
    j = j + 1
    
  } 
  
  return(clusters)
}

add.noise<- function(m, ...) {
  x = m
  v = x[upper.tri(x)]
  x[upper.tri(x)] = jitter(v, ...)
  v = diag(x)
  diag(x) = jitter(v, ...)
  x[x < 0] = 0
  x[lower.tri(x)] = t(x)[lower.tri(t(x))]
  return(x)
}


hbm.features<-function(m, noise.factor, ncores = 1, ref = NULL, ...)
{
  if (length(ref) == 0)
  {
    ref = hbm(m, ...)$hm
  }
  n = nrow(ref)
  if (length( noise.factor) == 0)
  {
    stop("factor must be a non empty vector")
  }else if (sum( noise.factor <= 0) > 0)
  {
    stop("factor must be a positive numeric vector")
  }

  l = list()
  if (ncores > 1)
  {
    ff = NULL
    registerDoParallel(ncores)
    on.exit( stopImplicitCluster())
    l <- foreach(ff=noise.factor, .export = c("hbm", "add.noise", "mcl"), .packages = c("Matrix")) %dopar% {
      x = add.noise(m, factor=ff)
      hbm(x, 2, ...)$hm  }
  }else{
    l = lapply(noise.factor, function(f) {
    x = add.noise(m, factor=f)
    hbm(x, 2, ...)$hm})
  }  
  target = Reduce('+', l)
  
  target = target/(length(noise.factor))
  maxscale = floor(max(as.vector(target)))
  target[target > maxscale] = maxscale
  features = target
  features[which(ref-target != 0)] = NA
  
  return(list("noisy.hm" = target, "features" = features))
}

hbm<-function(m, infl=2, ...)
{
  
  N = nrow(m)
  if (N != ncol(m))
  {
    stop("input matrix must be a square matrix")
  }
  hm = matrix(0, N, N)
  scale = 1
  previds = c()
  diff = 1
  scales = list()
  while(N > 1)
  {
    clusters = mcl(m, infl, ...)
    clustid = unique(clusters)
    ids = lapply(clustid, function(id) which(clusters == id))
    Ncurr = length(clustid)
    if (N - Ncurr == 0)
    {
      break
    }
    
    mtag = matrix(0, Ncurr, Ncurr)
    scales[[scale]] = ids
    for (i in 1:Ncurr)
    {
      for (j in i:Ncurr)
      {
        mtag[i,j] = mean(m[ids[[i]], ids[[j]]])
      }
    }
    mtag[lower.tri(mtag)] = t(mtag)[lower.tri(t(mtag))]
    scale = scale+1
    m = mtag
    N = Ncurr
  }
  
  if (length(scales) > 1)
  {
    for (scale in 2:length(scales))
    {
      currscale = scales[[scale]]
      prevscale = scales[[(scale-1)]]  
      for (i in 1:length(currscale))
      {
        indices = currscale[[i]]
        if (scale > 1)
        {
          indices = unlist(lapply(indices, function(x) prevscale[[x]]))
          currscale[[i]] = indices
        }
        scales[[scale]] = currscale
      }
    }
  }
  
  for (scale in length(scales):1)
  {
    currscale = scales[[scale]]
    for (i in 1:length(currscale))
    {
      indices = currscale[[i]]
      hm[indices,indices] = scale
    }

  }
  hm[hm == 0] = length(scales)+1
  diag(hm) = 0
  return(list("hm"=hm, "scales"=scales))
}

generate.random.conf<-function(n, k = 3, perturb = NULL, 
                               scale = T, mean = 0, sd = 1)
{
  m = matrix(0, n, k)
  m[1, 1:k] = rnorm(k)
  for (i in 2:n)
  {
    for (j in 1:k)
    {
      m[i,j] = m[i-1, j] + rnorm(1, mean = mean, sd = sd)
    }
  }
  
  
  if (scale)
  {
    scalevec = rep(1, k)
    for (i in 1:k)
    {
      scalevec[i] = abs(m[1,i]-m[n,i])
    }
    m = scale(m, center = FALSE, scale = scalevec)
  }
  
  
  if (length(perturb) > 0)
  {
    if (max(perturb) > n || length(perturb) > n || min(perturb) < 1 )
    {
      stop("perturb must be a vector of size > 1 and <=n with integers in the range 1 to n (inclusive)")
    }
    perturb = sort(perturb)
    indices = sample(perturb, length(perturb), replace = FALSE)
    trials = 5
    i = 0
    while (sum(indices != perturb) > 0  && i < trials ) # edge case
    {
      indices = sample(perturb, length(perturb), replace = FALSE)
      i = i + 1
    }
    m[perturb, 1:k] = m[indices,1:k]
  }
  return(m)
}

detect.movement<-function(ref, target, ref.res, target.res, 
                  motion.prop.thresh = 0.75, siglevel = 0.05, verbose = FALSE)
{
  
  n = nrow(ref)
  if (sum(dim(ref) != dim(target)) > 0)
  {
    stop("reference and target matrices must be of the same dimension")
  }else if (ncol(ref) != n)
  {
    stop("input matrices must be square matrices")
  }
  
  movement = matrix(0, n, n)
  
  
  ## transform contact maps to probability matrices
  ref = sweep(ref,2,colSums(ref),`/`)
  ref[is.na(ref)] = 0
  diag(ref) = NA
  target = sweep(target,2,colSums(target),`/`)
  target[is.na(target)] = 0
  diag(target) = NA
  
  # extract scales and hierarchical matrix from the data structure
  ref.scales = ref.res$scales
  target.scales = target.res$scales
  ref.hm = ref.res$hm
  target.hm = target.res$hm
  
  maxscale = length(ref.scales) 
  
  # iterate through the scales and the clusters of the refernce
  # and look for changes, then map them to directed movement
  for (scale in 1:maxscale)
  {
    clusters = ref.scales[[scale]]
    for (i in 1:length(clusters))
    {
      sig.pack = FALSE
      sig.unfold = FALSE
      ids = clusters[[i]]
      
      # only check clusters with size > 1 that 
      # their maximal scale is the scale under consideration
      if (max(ref.hm[ids,ids]) == scale & length(ids) > 1) 
      {
        
        # check if the cluster has become significantly more packed or unfolded
        # this acchived by testing the relative within and
        #  outer connectivity in the ref and target
        notids = which(!(1:n %in% ids))
        ref.conn = unlist(lapply(ids, function(x) mean(ref[ids,x], na.rm = T)))/unlist(lapply(ids, function(x) mean(ref[notids,x], na.rm = T)))
        target.conn = unlist(lapply(ids, function(x) mean(target[ids,x], na.rm = T)))/unlist(lapply(ids, function(x) mean(target[notids,x], na.rm = T)))        
        answer = tryCatch(
          {
            pval.g = t.test(ref.conn, target.conn, alternative = "greater")$p.value
            
            if (!is.na(pval.g) && pval.g < siglevel)
              #target cluster is less compact than ref cluster
            {
               sig.unfold = TRUE
               if (verbose)
               {
                 cat(paste("cluster", i, "at scale", scale, 
                           "is less compact in target (p-value <",pval.g,"\n",sep=" "))
               }
            }
            else 
            {
              pval.l = t.test(ref.conn, target.conn, alternative = "less")$p.value
              if (!is.na(pval.l) && pval.l < siglevel)
              {
                #target cluster is more compact than ref cluster
                sig.pack = TRUE
                if (verbose)
                {
                  cat(paste("cluster", i, "at scale", scale, 
                            "is more compact in target (p-value <",pval.l,"\n",sep=" "))
                }
                # we cannot detect orientation so we make it symmetric
                movement[ids, ids] = 1
               
              }
            }
          },error = function(w){}
        )
        
        # now inspect for local movements 
        c1 = ref.hm[ids,ids]
        c2 = target.hm[ids,ids]
        
        # look only within the scale of interest to avoid reps.
        c1[which(c1 != scale)] = 0
        c2[which(c1 == 0)] = 0
        
        # get the scales in the target cluster (and remove 0 that is always on the diagonal)
        unique.target.scales = sort(unique(as.vector(c2)))
        unique.target.scales = unique.target.scales[2:length(unique.target.scales)]
        
        # detect possible motion: indicated by multiple scales at the target cluster 
        min.scale = unique.target.scales[1]
        target.c2 = target.hm
        target.c2[ids,ids] = c2
        if (length(unique.target.scales) > 1)
        {
          prop.moved = unlist(lapply(ids, function(x) sum(target.c2[ids,x] > min.scale)/length(ids)))
          if (sum(prop.moved >= motion.prop.thresh) > 0) # possible movement
          {
            if (sig.pack) # movement because of significant within-cluster packing
            {
              clustered = ids[which(prop.moved < motion.prop.thresh)]
              moved = ids[which(prop.moved >= motion.prop.thresh)]
              movement[clustered, moved] = -1
              # we cannot be sure if any others have moved away additionally so we put 0.5
              for (region in moved)
              {
                attractor = which(target[,region] == max(target[,region], na.rm = T))[1]
                movement[region, attractor] = 0.5
              }
             
            }else 
            {
              moved = ids[which(prop.moved >= motion.prop.thresh)]
              for (region in moved)
              {
                attractor = which(target[,region] == max(target[,region], na.rm = T))[1]

                if (attractor %in% ids)
                {
                  if (sig.unfold) # this can happen becuase of movements away within the cluster
                  {
                    movement[region, attractor] = 1
                    others = ids[!(ids %in% c(region, attractor))]
                    vals = ref[region, others]
                    repulsor = others[which(vals == max(vals, na.rm = T))]
                    movement[region, repulsor] = -1
                    movement[region, others[others!=repulsor]] = -0.5
                  }
                  
                }else
                {
                  movement[region, attractor] = 1
                  movement[region, which(target.hm[attractor,] == 1)] = 0.5
                  others = ids
                  vals = ref[region, others]
                  repulsor = others[which(vals == max(vals, na.rm = T))]
                  movement[region, repulsor] = -1
                  movement[region, others[others!=repulsor]] = -0.5
                }
              }
            }
            
          }
        } # end if of local movements
        
      }# end big if
      
    }#end cluster for
  }#end scale for
  diag(movement) = 0
  return(movement)
}

get.movements<-function(movement, hm, features = NULL)
{
  
  movement.summary = which(movement != 0, arr.ind = T)
  indices = which(movement != 0)
  type = movement[indices]
  scale = hm[indices]
  robust = c()
  if (length(features) > 0)
  {
    if (is.matrix(features))
    {
      robust = features[indices]
      robust[is.na(robust)] = 0
      robust[robust > 0] = 1
    }else if (is.list(features))
    {
      for (i in 1:length(features))
      {
        robust.i = features[indices]
        robust.i[is.na(robust.i)] = 0
        robust.i[robust.i > 0] = 1
        robust = cbind(robust, robust.i)
      }
      colnames(robust) = unlist(lapply(1:length(features), function(x) paste("robust", x, sep = "")))
    }
   
  }
  
  movement.summary = data.frame("from" = movement.summary[,1], "to" = movement.summary[,2], "type" = type, "scale" = scale) 
  if (length(robust) > 0)
  {
    movement.summary = cbind(movement.summary, "robust" = robust)
  }
  return(movement.summary)
  
}


hierarchy<-function(hm)
{

  n = nrow(hm)
  pi = unlist(lapply(1:n, function(i)
    {
    x = hm[i,]
    if (i < 3)
    {
      x2 = x[(i+1):n]
      diff2 = x2[2:(length(x2)-1)]-x2[1:(length(x2)-2)]
      s2 = sum(diff2 < 0)
      s2
    }else if (i > n-2)
    {
      x1 = x[1:i]
      x1 = x1[x1 > 0]
      diff1 = x1[2:(length(x1)-1)]-x1[1:(length(x1)-2)]
      s1 = sum(diff1 > 0) 
      s1
     
    }else
    {
      x1 = x[1:i]
      x1 = x1[x1 > 0]
      diff1 = x1[2:(length(x1)-1)]-x1[1:(length(x1)-2)]
      s1 = sum(diff1 > 0)
      x2 = x[(i+1):n]
      diff2 = x2[2:(length(x2)-1)]-x2[1:(length(x2)-2)]
      s2 = sum(diff2 < 0) 
      s1+s2
      
    }
  }))
  
  return(-mean(pi))

}

communicability<-function(hm)
{
  n = nrow(hm)
  ns = max(hm)-1
  comm = matrix(0, n, n)
  for (s in 1:ns)
  {
    A = hm
    A[A > s] = 0
    comm = comm + as.matrix(expm(A))
  }
  comm = comm/ns
  return(comm)
}




