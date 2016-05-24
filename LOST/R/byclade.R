byclade <-
function(x,remperc,ngroups,groups) {
    
  remove.dat<-function(specimen,removes){
    ndat<-length(specimen)
    rems<-sample(ndat,removes,replace=FALSE)
    for (k in 1:removes){
      m<-rems[k]
      specimen[m]<-NA
    }
    return(specimen)
  }
  
  
  newx1<-as.matrix(x)
    grouping<-as.factor(groups)
    newx2<-as.matrix(x)
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
    numbersp<-ifelse(numberper==0,0,1)*rows
    nsp<-length(numbersp)
    sorted<-sort(numberper,decreasing=TRUE)
    splitgroups<-split(as.data.frame(x),grouping)
    npergroup<-sapply(splitgroups,nrow,simplify=TRUE)
    counts<-1:nsp
    for (i in 1:nsp){
      m<-groups[i]
      a<-npergroup[m]
      counts[i]<-a
    }
    counts<-counts
    sums<-sum(npergroup)
    ratio<-sums/counts
    probs<-ratio/sum(ratio)
    orders<-sample(1:nsp,nsp,replace=FALSE,prob=probs)
    for (k in 1:nsp){
      removes<-sorted[k]
      spnumber<-orders[k]
      specimen<-newx2[spnumber,]
      newsp<-remove.dat(specimen,removes)
      newx2[spnumber,]<-newsp
    }
    return(newx2)
  }
