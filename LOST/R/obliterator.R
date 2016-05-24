obliterator <- function(x,remperc,landmarks,expo=1) {
  
  
  
  probability.generator<-function(newsp,distances,removes,expo=1,sa){
    zeros<-0
    if (sum(removes)==0){nrcur<-1;current=0} else
    {current<-setdiff(removes,zeros); nrcur<-length(current)}
    probmatrix<-rep(1,2*length(newsp))
    for (m in 1:(nrcur)){
      distancesx1<-distances[1,]
      distancesx2<-distances[2,]
      distancesy1<-distances[3,]
      distancesy2<-distances[4,]
      distancesz1<-distances[5,]
      distancesz2<-distances[6,]
      distancesx<-rbind(distancesx1,distancesx2)
      distancesy<-rbind(distancesy1,distancesy2)
      distancesz<-rbind(distancesz1,distancesz2)
      if (sum(removes)==0){anchor<-0} else
      {anchor<-current[m]}
      if (anchor==0){sss<-1} else
      {sss<-sa[m]}
      if (anchor==0){basex<-0} else
      {basex<-distancesx[sss,anchor]}
      if (anchor==0){basey<-0} else
      {basey<-distancesy[sss,anchor]}
      if (anchor==0){basez<-0} else
      {basez<-distancesz[sss,anchor]}
      distsx1<-sqrt((distancesx1-basex)^2)
      distsx2<-sqrt((distancesx2-basex)^2)
      #
      distsy1<-sqrt((distancesy1-basey)^2)
      distsy2<-sqrt((distancesy2-basey)^2)
      #
      distsz1<-sqrt((distancesz1-basez)^2)
      distsz2<-sqrt((distancesz2-basez)^2)
      #
      distsstart<-sqrt(distsx1^2+distsy1^2+distsz1^2)
      distsstop<-sqrt(distsx2^2+distsy2^2+distsz2^2)
      dists<-c(distsstart,distsstop)
      #
      nozeros<-ifelse(dists==0,(max(dists)*10),dists)
      dists<-ifelse(dists==0,(min(nozeros)/2),dists)
      inv<-1/(dists^expo)
      inv[anchor]<-0
      ll<-length(inv)/2
      inv[(anchor+ll)]<-0
      sums<-sum(inv)
      ones<-rep(1,length(inv))
      if (anchor==0){probs<-ones} else
      {probs<-inv/sums}
      checker<-c(newsp,newsp)
      probs<-ifelse(is.na(checker),0,probs)
      probmatrix<-rbind(probmatrix,probs)
    }
    probs<-apply(probmatrix,2,prod)
    sums<-sum(probs)
    probs<-probs/sums
    return(probs)
  }
  
  
  
  
  remove.points<-function(specimen,r,distances,expo){
    removes<-rep(0,r)
    newsp<-specimen
    l<-length(specimen)
    nl<-1:(2*l); anchor<-0
    site<-c(1:l,1:l)
    startorstop<-c(rep(1,l),rep(2,l))
    sa<-removes
    for (k in 1:r){
      probs<-probability.generator(newsp,distances,removes,expo,sa)
      a<-sample(nl,1,prob=probs)
      b<-site[a]
      cc<-startorstop[a]
      newsp[b]<-NA
      removes[k]<-b
      sa[k]<-cc
    }
    return(newsp)
  }
  
  
    newx1<-as.matrix(x)
    distances<-as.matrix(landmarks)
    totaldata<-nrow(x)*ncol(x)
    n<-round(totaldata*remperc)
    ndat<-1:totaldata
    remove<-sample(ndat,n,replace=FALSE)
    for (k in 1:n) {
      i<-remove[k]
      newx1[i]<-NA
    }
    binary<-ifelse(is.na(newx1),1,0)
    numberper<-apply(binary,1,sum)
    rows<-1:nrow(x)
    numbersp<-setdiff((ifelse(numberper==0,0,1)*rows),0)
    nsp<-length(numbersp)
    sa<-rep(0,ncol(x))
    newx<-x
    for (k in 1:nsp){
      i<-numbersp[k]
      specimen<-x[i,]
      r<-numberper[k]
      newsp<-remove.points(specimen,r,distances,expo)
      newx[i,]<-newsp
    }
    return(newx)
  }
