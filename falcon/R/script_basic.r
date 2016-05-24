## old functions from seqCBS
getAutoGridSize <-
function(nL) {
	index10 = floor(log(nL, base=10))
	if(nL/(10^index10) < 3) {
		index10 = index10-1
	}
	grid.size = 10^(1:index10)
	return(grid.size)
}

## new functions
hardthres = function(v, low=0.9, high=1.1){
  n = length(v)
  for (i in 1:n){ if (v[i]>low && v[i]<high) v[i] = 1 }
  v
}

comb = function(tau, cns1, cns2, low=0.9, high=1.1){
  n = length(cns1)
  if (n==1) return(tau)
  code = rep("",n)
  for (i in 1:n){
    if (cns1[i]<low){
      code[i] = "0"
    }else if (cns1[i]>high){
      code[i] = "2"
    }else{
      code[i] = "1"
    }
    if (cns2[i]<low){
      code[i] = paste(code[i],"0",sep="")
    }else if (cns2[i]>high){
      code[i] = paste(code[i],"2",sep="")
    }else{
      code[i] = paste(code[i],"1",sep="")
    }
  }
  removeid = c()
  for (i in 2:n){
    if (code[i]==code[i-1]){
      removeid = c(removeid, i)
    }
  }
  if (length(removeid)==0) return(tau)
  return(tau[-removeid])
}


plotCN = function(n, tauhat, ascn, pos=NULL, gaincol="red", losscol="blue", neutralcol="green",xlab=NULL, ylab=NULL, pch=".",...){
  tauhat = sort(unique(c(1,tauhat,n)))
  ascn1 = ascn[1,]
  ascn2 = ascn[2,]
  if (is.null(pos)){
    pos = 1:n
    if (is.null(xlab)) xlab = "SNP #"
  }else{
    tauhat = pos[tauhat]
    if (is.null(xlab)) xlab = "Position (bp)"
  }
  if (is.null(ylab)) ylab = "Allele-specific CN"
  poscn1 = poscn2 = rep(1,n) 
  K = length(tauhat)-1
  m = match(tauhat[1:K], pos)
  if (K>1){
    for (i in 1:(K-1)){
      poscn1[m[i]:(m[i+1]-1)] = ascn1[i]
      poscn2[m[i]:(m[i+1]-1)] = ascn2[i]
    }
  }
  poscn1[m[K]:n] = ascn1[K]
  poscn2[m[K]:n] = ascn2[K]
  g1 = l1 = g2 = l2 = c()
  if (K>1){
    for (i in 1:(K-1)){
      if (ascn1[i]>1){
        g1 = c(g1, m[i]:(m[i+1]-1))
      }else if (ascn1[i]<1){
        l1 = c(l1, m[i]:(m[i+1]-1))
      }
      if (ascn2[i]>1){
        g2 = c(g2, m[i]:(m[i+1]-1))
      }else if (ascn2[i]<1){
        l2 = c(l2, m[i]:(m[i+1]-1))
      }
    }
  }
  plot(pos, poscn1, col=neutralcol, ylim = c(0, max(c(ascn1,ascn2))+0.5),xlab=xlab, ylab=ylab, pch=pch, ...)
  points(pos, poscn2, col=neutralcol, pch=pch, ...)
  for (i in 1:K){
    if (i==K){
      ids = m[K]:n
    }else{
      ids = m[i]:(m[i+1]-1)
    }
    nids = length(ids)
    if (ascn1[i]>1){
      points(pos[ids], poscn1[ids], col=gaincol, pch=pch,...)
      if (i>1 && ascn1[i-1]>1){
        a = ascn1[i-1]
      }else{
        a = 1
      }
      points(c(pos[ids[1]], pos[ids[1]]), c(a, ascn1[i]), col=gaincol, type="l")
      if (i<K-1 && ascn1[i+1]>1){
        a = ascn1[i+1]
      }else{
        a = 1
      }
      points(c(pos[ids[nids]], pos[ids[nids]]), c(a, ascn1[i]), col=gaincol, type="l")
    }else if (ascn1[i]<1){
      points(pos[ids], poscn1[ids], col=losscol, pch=pch, ...)
      if (i>1 && ascn1[i-1]<1){
        a = ascn1[i-1]
      }else{
        a = 1
      }
      points(c(pos[ids[1]], pos[ids[1]]), c(a, ascn1[i]), col=losscol, type="l")
      if (i<K-1 && ascn1[i+1]<1){
        a = ascn1[i+1]
      }else{
        a = 1
      }       
      points(c(pos[ids[nids]], pos[ids[nids]]), c(a, ascn1[i]), col=losscol, type="l")
    }
    if (ascn2[i]>1){
      points(pos[ids], poscn2[ids], col=gaincol, pch=pch, ...)
      if (i>1 && ascn2[i-1]>1){
        a = ascn2[i-1]
      }else{
        a = 1
      }
      points(c(pos[ids[1]], pos[ids[1]]), c(a, ascn2[i]), col=gaincol, type="l")
      if (i<K-1 && ascn2[i+1]>1){
        a = ascn2[i+1]
      }else{
        a = 1
      }
      points(c(pos[ids[nids]], pos[ids[nids]]), c(a, ascn2[i]), col=gaincol, type="l")
    }else if (ascn2[i]<1){
      points(pos[ids], poscn2[ids], col=losscol, pch=pch, ...)
      if (i>1 && ascn2[i-1]<1){
        a = ascn2[i-1]
      }else{
        a = 1
      }
      points(c(pos[ids[1]], pos[ids[1]]), c(a, ascn2[i]), col=losscol, type="l")
      if (i<K-1 && ascn2[i+1]<1){
        a = ascn2[i+1]
      }else{
        a = 1
      }
      points(c(pos[ids[nids]], pos[ids[nids]]), c(a, ascn2[i]), col=losscol, type="l")
    }
  }
}


