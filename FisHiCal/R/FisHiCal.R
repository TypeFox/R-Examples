# find matching indices in hic genomic Coordinates
findMatchingIndices<-function(fishCoord, hicCoord)
{
  n = nrow(fishCoord)
  indices = vector(length = n, mode = "integer")
  midPoints = (fishCoord$start + fishCoord$end)/2 
  for (i in 1:n)
  {
    mid = midPoints[i]
    ind = which(as.character(hicCoord$chr) == as.character(fishCoord[i,]$chr) & 
                  (hicCoord$start <= mid & hicCoord$end > mid))
    if (length(ind) > 0){indices[i] = ind}
    
  }
  return(indices)  
}

prepareData<-function(fish, fishCoord, hic, hicCoord)
{
  indices = findMatchingIndices(fishCoord, hicCoord)
  uindices = unique(indices)  
  un = length(uindices)
  frequencies = c()
  distances = c()
  for (i in 1:un)
  {
    ind1 = uindices[i]
    if (ind1 > 0)
    {
      f1 = which(indices == ind1)# matching FISH indices
      for (j in i:un)
      {
        ind2 = uindices[j] 
        if (ind2 > 0)
        {
          freq = as.numeric(hic[ind1, ind2])
          if (freq > 0)
          {
            f2 = which(indices == ind2)# matching FISH indices
            m = fish[f1, f2]
            d = sum(m)
            if (d > 0)
            {
              d = as.numeric(min(m[m>0]))
              distances = c(distances, d)
              frequencies = c(frequencies, freq)
            }
          }
          
            
        }
      }
    }
   
  } 

  data = as.data.frame(cbind(distances, frequencies))
  data = unique(data)
  data = data[order(data$distances),]
  return(data)
}

prepareCalib<-function(data, npoints, threshold = NULL, useMax = TRUE, delta = 0.05, buffer = 1)
{
  data = data[order(data$distances),]
  xpoints = c()
  if (length(npoints) > 1)
  {
	xpoints = npoints
  }else{
	xpoints = 1:npoints
  }
  fit = lm(log(data[xpoints,]$frequencies)~log(data[xpoints,]$distances))
  a = fit$coefficients[[2]] # slope
  b = fit$coefficients[[1]]
  
  # prepare f - the calibration function 
  if (length(threshold) == 0)
  {
    if (!useMax)
    {
      dev = abs(exp((log(data$frequencies) - b)/a) - data$distances)
      devd = which(dev < delta)
      if (length(devd) > 0 )
      {
        threshold = max(data$distances[devd])
      }else
      {
        useMax = T
      }
    }
    if(useMax)
    {
      threshold = max(data$distances) + buffer
    } 
  }
  params = list(a, b, threshold)
  f <- function(m, params)
  {
    a = params[[1]]
    b = params[[2]]
    threshold = params[[3]]
    m = exp((log(m) - b)/a)
    m[m == Inf] = 0 
  	m[m > threshold] = 0
    return(m)
  }
  calib = list("f" = f, "params" = params)
  return(list("calib" = calib, "fit" = fit))
}

updateCalib<-function(calib, paramVal, paramIndex)
{
  calib$params[[paramIndex]] = paramVal
  return (calib)
}

calibrate<-function(hic, calib)
{
  calibMat = calib$f(hic, calib$params)
  return(calibMat)
}

getInfoLevelForChr<-function(calibHiC, hicCoord, chr)
{
  indices = which(hicCoord$chr == chr)
  infoLevel = NULL
  if (length(indices) > 0)
  {
    m = calibHiC[indices,indices]
    ut = m[upper.tri(m, diag = F)]
    infoLevel = sum(ut > 0)/length(ut)
  }
  return(infoLevel)
}

lsmacof<-function(diss, infD, itermax = 10000, eps = 1e-06, init = NULL, k = 3, 
                          verbose = FALSE, infW = NULL)
{
  if (length(init) == 0)
  {
    init = cmdscale(diss, k = k)
  }
  # updates the dissimilarity matrix and builds the weight matrix for local stress
  n = nrow(diss)
  w = matrix(1, nrow = n, ncol = n)
  x = which(diss <= 0)
  if (length(infW) > 0)
  {
    w[x] = infW
  }
  else if (infD <= 1)
  {
    w[x] = 0.05
  }
  else
  {
    w[x] = 1/infD
  }
  diss[x] = infD
  diag(w) = 0
  diag(diss) = 0
  conf = smacof(diss, w, init, itermax, eps, verbose) # cpp impl
  return(conf)
}



searchInc<-function(calibHiC, hicCoord)
{
  diag(calibHiC) = 0
  colnames(calibHiC) = 1:ncol(calibHiC)
  g = graph.adjacency(calibHiC, mode="undirected", weighted=TRUE, diag=TRUE)
  chrs = unique(as.vector(hicCoord$chr))
  splits = list(length = nrow(calibHiC)) # for each locus we save its neighborhood if it's
  # splitted - i.e. not a connected graph
  
  sIndices = c() # indices of loci that have a splitted neigborhood 
  j = 1
  for (refChr in chrs)
  {
    indices = which(hicCoord$chr == refChr)
    i1 = indices[1]
    iN = indices[length(indices)]
    # iterate through all the neigborhoods of the loci in teh current chromosome
    neighbors = neighborhood(g, order=1, nodes=V(g)[i1:iN], mode="all")
    for (i in 1:length(indices))
    {
      gNgr = NULL
      myNgr = neighbors[[i]]
      # get set of neighbors in other chromosomes (trans) and check if connected
      transNgr = myNgr[which(hicCoord[myNgr,]$chr != refChr)]
      if (length(transNgr) > 1)
      { 
        gNgr = induced.subgraph(g, transNgr) # induced graph of trans neighbors
        if (!is.connected(gNgr))
        { 
          v = sort(c(indices[i], transNgr))
          myIndex = which(v == indices[i]) 
          comp = clusters(gNgr)
          n = length(v)
          membership = vector(mode = "integer", length = n)
          if (myIndex == 1)
          {
            membership[2:n] = comp$membership 
            membership[1] = 0
          }else if (myIndex == n)
          {
            membership[n] = 0
            membership[1:(n-1)] = comp$membership
          }else{
            membership[1:(myIndex-1)] = comp$membership[1:(myIndex-1)]
            membership[(myIndex+1):n] = comp$membership[myIndex:(n-1)]
            membership[myIndex] = 0
          }
          
          
          gNgr = induced.subgraph(g, v)
          fullNames = unlist(lapply(v,  function(i) paste(hicCoord[i,]$chr, ":", hicCoord[i,]$start, "-", hicCoord[i,]$end, sep = "")))
          chrNames = unlist(lapply(v,  function(i) hicCoord[i,]$chr))
          
          gNgr = set.vertex.attribute(gNgr, name = "membership", index=V(gNgr), value = membership)
          gNgr = set.vertex.attribute(gNgr, "fullName", index=V(gNgr), fullNames)
          gNgr = set.vertex.attribute(gNgr, "chr", index=V(gNgr), chrNames)
          sIndices = c(sIndices, j)
          
        }
      }
      splits[[j]] = gNgr
      j = j + 1
    }
    
  }
  return(list("neighborhoods" = splits, "indices" = sIndices))
  
}

summaryInc<-function(indices, neighborhoods)
{
  inconsistency = c()
  if (length(indices) > 0)
  {
    for (i in indices)
    {
      g = neighborhoods[[i]]
      if (length(g) > 0)
      {
          chr = (get.vertex.attribute(g, "chr", index = V(g)))
          component = (get.vertex.attribute(g, "membership", index = V(g)))
          index = (get.vertex.attribute(g, "name", index = V(g)))
          fullName = (get.vertex.attribute(g, "fullName", index = V(g)))
          splitIndex = index[which(component == 0)]
          splitChr = chr[which(component == 0)]
          
          myInc = as.data.frame(cbind(chr, component, index, fullName, splitChr, splitIndex))
          myInc = myInc[order(component),] 
          inconsistency = rbind(inconsistency, myInc) 
      }
    }
  }
  return(inconsistency)
}

plotInc<-function(i, neighborhoods, label = "fullName", size = 20, interactive = F)
{
  
  g = neighborhoods[[i]]
  if (length(g) > 0)
  {
    if (!(label %in% c("name", "fullName", "chr")))
    {
      label = "name"
    }
    memberAtt = get.vertex.attribute(g, "membership", index = V(g))
    cols = rainbow(length(unique(memberAtt)) -1 )
    vertexcols = unlist(lapply(memberAtt, function(x){ if (x == 0){"white"}else{cols[x]}} ))
    if (!interactive)
    {
      plot(g, vertex.label = get.vertex.attribute(g, label, index = V(g)), vertex.color = vertexcols, vertex.size = size) 
    }else
    {
      tkplot(g, vertex.label = get.vertex.attribute(g, label, index = V(g)), vertex.color = vertexcols, vertex.size = size) 
      
    }
    
  }
  
}

