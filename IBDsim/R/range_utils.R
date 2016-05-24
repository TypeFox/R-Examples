rangeUnion <-
function(x) {
  nrws = lengths(x, use.names=FALSE)/2
  if (sum(nrws)==0) return(matrix(numeric(0), ncol=2))
  
  n = length(nrws)
  vec = unlist(x, recursive=TRUE, use.names=FALSE)
  startstop = rep.int(rep.int(c(1,-1), n), rep.int(nrws, rep.int(2,n))) #1 = interval start; -1 = interval stop
  extvec = vec - startstop*1e-6 #extend all regions by 1bp in both directions (to count consecutive regions as overlapping) 
  
  ord = order(extvec, startstop, na.last=TRUE, decreasing=FALSE) #NB: in rangeIntersect this must be '-startstop'
  m = extvec[ord]
  startstop.ord = startstop[ord]
  
  #print(cbind(m,startstop.ord,cumsum(startstop.ord)))  #look at this to understand.
  ind.end = seq_along(ord)[cumsum(startstop.ord) == 0]
  ind.start = c(1, ind.end[-length(ind.end)] + 1)
  cbind(m[ind.start]+1e-6, m[ind.end]-1e-6, deparse.level=0)
}

rangeIntersect <-
function(x) {
  nrws = lengths(x, use.names=FALSE)/2
  if (sum(nrws)==0) return(matrix(numeric(0), ncol=2))
  
  n = length(nrws)
  vec = unlist(x, recursive=TRUE, use.names=FALSE)
  startstop = rep.int(rep.int(c(1,-1), n), rep.int(nrws, rep.int(2,n)))

  ord = order(vec, -startstop, na.last=TRUE, decreasing=FALSE)
  m = vec[ord]
  startstop.ord = startstop[ord]

  ind.start = seq_along(vec)[cumsum(startstop.ord)==n]
  ind.end = ind.start+1
  cbind(m[ind.start], m[ind.end], deparse.level=0)
}


rangeCompl <-
function(m, endpoint) { #m a matrix with columns start & end
  n = nrow(m)
  if (n==0) return(cbind(1, endpoint, deparse.level=0))
  newstart = m[-n, 2]
  newend = m[-1, 1]
  
  if(m[1,1] > 0) { 
	newstart = c(0, newstart); newend = c(m[1, 1], newend)
  }
  if(m[n, 2] < endpoint) {
	newstart = c(newstart, m[n, 2]); newend = c(newend, endpoint)
  }
  res = cbind(newstart, newend, deparse.level=0)
  res[res[, 1] < res[, 2], , drop=F]
}