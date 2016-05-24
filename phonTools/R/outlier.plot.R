# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

outlier.plot = function (x, y, category, xsampa = TRUE, logaxes = TRUE, ellipsesd = 2, borders = c(2,3), select = 0, nearest = 1){

  ffs = as.matrix(cbind(x,y))
  if (logaxes) ffs = log(ffs)
  
  temp = createtemplate (ffs,category)
  cs = as.factor(category)
  lcs = levels (cs)
  ncs = length (lcs)
  type = as.numeric(cs)
  
  tmp.env = environment()
  
  dist = zeros(x)
  for (i in 1:nrow(temp$means)){ 
    use = (category == lcs[i])
    dist[use] = mahalanobis (ffs[use,], temp$means[i,], cov (ffs[use,]))
  }
  dist = sqrt(dist)
  
  cols = rep('', length(dist))
  cols[dist<borders[1]] = 'forestgreen';
  cols[dist>borders[1] & dist < borders[2]] = 'gold3';
  cols[dist>borders[2]] = 'firebrick';
  
  sizes = zeros(dist)
  sizes[dist<borders[1]] = .5;sizes[dist>borders[1] & dist < borders[2]] = 1.2;sizes[dist>borders[2]] = 1.7;

  #vns = c('x','y','category', 'xsampa', 'logaxes', 'sizes')
  #vs = list(x,y,category, xsampa, logaxes, sizes)
  #for (i in 1:6) assign (vns[i],vs[[i]],envir=.GlobalEnv)
  
  oldpar = par()
  par (mar = c(4.1,4.1,1,1))
  vplot (x,y,category, xsampa = xsampa, logaxes = logaxes, cex = sizes, colors = cols,
  xlab='Dimension 1',ylab='Dimension 2')
  vplot (x,y,category, xsampa = xsampa, logaxes = logaxes, colors = 1, 
         add = TRUE, meansonly = TRUE, cex = 2.5)
  
  for (i in 1:nrow(temp$means)){  
      if (!logaxes) sdellipse (cbind (x[cs==lcs[i]],y[cs==lcs[i]]), 
                              stdev = ellipsesd, col = 1,lwd=2,lty='dotted') 
      
      if (logaxes){ tmp = sdellipse (log(cbind (x[cs==lcs[i]],y[cs==lcs[i]])), 
                          stdev = ellipsesd, show = F); lines (exp(tmp), col = 1,
						  lwd=2,lty='dotted')}
  }  
  suppressWarnings (par (oldpar))
  if (select > 0){
    coords = locator(select) 
    coords = log(as.matrix(cbind(coords$x,coords$y)))
    siginv = solve (cov (ffs))  
  
    index = NULL
    selection = NULL
    closest = NULL
  
    for (i in 1:select){
      dists = mahalanobis (ffs, coords[i,], cov (ffs))
      index = c(index, order(dists)[1:nearest])
      selection = c(selection, rep(i, nearest))
      closest = c(closest, 1:nearest)
    }  
    return (data.frame (index, selection, closest))  
  }
}

