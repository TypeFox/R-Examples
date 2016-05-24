
# .estimateWindow <- function(y){
#   #Get middle 60%
#   m = quantile(1:length(y), c(0.25, 0.75))
#   y = y[m[1]:m[2]]
#   w = length(y) / (which.max(Re(fft(y-mean(y)))[1:(length(y)/2)]) - 1)
#   w
# }

.estimateWindow2 <- function(x, bw=15){
  m = quantile(1:length(x), c(0.25, 0.75))
  x = x[m[1]:m[2]]
  xo = x
  x = ksmooth(1:length(x), x, kernel='normal',bandwidth=bw)$y
  pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0) + 2
  
#   plot(xo, type='l')
#   lines(x, col='blue')
#   abline(v=pks, col='red')
  
  pks = sort(diff(pks), decreasing=T)
  pks = pks[1:(length(pks)/2)]
  return(median(pks, na.rm=T))
}

.sinCor <- function(x,w){
  # Smooth
  x = ksmooth(1:length(x), x, kernel='normal',bandwidth=15)$y
  
  # Reference curve
  ref = sin( base::seq(-pi, 2*pi, length.out=w) )
  # Correlate
  # cr = mclapply
  cr = lapply(1:length(x), function(i){
    s = i - (w/2)
    e = i + (w/2)
    if(s < 1) s = 1
    if(e > length(x)) e = length(x)
    z = x[s:e][1:length(ref)]
    suppressWarnings( cor(ref, z, method="spearman") )
  })
  cr = simplify2array(cr)
  # Remove less then 0.3 correlations
  cr[is.na(cr) | cr < 0.3] = 0
  cr
}


.colonyPeaks <- function(x, n, format, plot=T){
  bw = 15
  if(format < 500) bw = 40
  w = .estimateWindow2(x, bw)
  loginfo('Estimated window size: %s', w)
  
  x = as.vector(x*runmed(x,.roundOdd(w/3)))
  
  max.win = (length(x) / n)
  if(w > max.win){
    w = floor(max.win)
    loginfo('Window too large, changing to %s', w)
  }
  
  x = c(rep(0, w), x, rep(0, w))
  cr = .sinCor(x, w)
  cr = cr[(w+1): (length(cr)-w)]
  
  rn = 1:length(cr)
  pe = .getPeaks(cr, floor( w/3 ))
  rp = rn[pe]
  pv = cr[pe]
  
  if(length(rp) < n)
    stop('Not enough peaks found')
  
  peak.dist = .shift(rp, 1) - rp
  peak.height = .shift(pv, 1) - pv
  #peak.height = .shift(x[pe], 1) - x[pe]
  
  num.s = (length(rp) - n) + 1
  P = sapply(1: num.s, function(s.i){
    dist = peak.dist[s.i:(s.i + (n-2))]
    w = peak.height[s.i:(s.i + (n-2))]
    sum(log( dnorm(dist, median(dist), sd(dist)) * dnorm(w, median(w), sd(w)) ) )
  })
  
  s = which.max(P)
  peaks = rp[s:(s + (n-1))]
  
  # Compute new delta
  delta = median( .shift(peaks, 1) - peaks, na.rm=T )
  
  if(plot){
    plot(cr/max(cr), lwd=1, ylab='Sum of pixel intensities', xlab='Index', bty='o', type='l')
    abline(v=rp, col="red")
    abline(v=peaks, col='#2e756d', lwd=2)
    lines(x/max(x), col='blue', lwd=2)
    text(peaks, quantile(cr,.1, na.rm=T), as.character(1:length(peaks)))
  }
  return(list(peaks=peaks, all.peaks=rp, window=delta/2))
}

.getWindow <- function(data, pos, window){
  win.left = window
  win.right = window
  if(pos - window < 0) win.left = pos
  if(pos + window > length(data)) win.right = length(data) - pos
  
  return(data[(pos-win.left):(pos+win.right)])
}

# Shift vector by an offset
.shift <- function(x, offset, na.pad=NA){
  t = (1:length(x)) + offset
  t[t<1] = na.pad
  return(x[ t ])
}


.getPeaks <- function(x, halfWindowSize, type="max") {
  if(type=="min")
    x = 1/x
  
  windowSize <- halfWindowSize * 2 + 1
  windows <- embed(x, windowSize)
  localMaxima <- max.col(windows, "first") == halfWindowSize + 1
  
  return(c(rep(FALSE, halfWindowSize), localMaxima, rep(FALSE, halfWindowSize)))
}


.splitHalf <- function(vec){
  t = ceiling(length(vec)/2)#middle
  return(list(left=vec[1: (t-1)], right=vec[(t+1):length(vec)]))
}

.fixOuterPeaks <- function(peaks){
  n = length(peaks)
  left = seq(peaks) < n/2
  p1 = peaks[left]
  p2 = peaks[!left]
  
  dist = .shift(peaks, 1) - peaks
  d1 = .shift(p1, 1) - p1
  d2 = p2 - .shift(p2, -1)
  
  s = 4*sd(dist[(0.2*length(dist)):(0.8*length(dist))], na.rm=T)
  m = median(dist, na.rm=T)
  
  f = 1
  for(i in f){
    if(d1[i] < m-s | d1[i] > m+s)
      peaks[1] = round(peaks[2] - m)
    j = length(d2) - (i-1)
    if(d2[j] < m-s | d2[j] > m+s)
      peaks[length(peaks)] = round(peaks[n-i] + m)
  }
  return(peaks)
}
