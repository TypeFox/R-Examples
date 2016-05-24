test.combinations <-
function(SpA, SpB){
  max.alleles<-length(SpA)
  SpA<-SpA[!is.na(SpA)]
  SpB<-SpB[!is.na(SpB)]
  total.alleles<-length(SpA) ## number of !is.na alleles
  max.delta<-0
  obs.delta<-delta(SpA, SpB)
  search.complete<-FALSE
  
  ## gives the alleles with a frequency greater than zero in each population
  spa.alleles<-which(SpA > 0)
  spb.alleles<-which(SpB > 0)

  ## sets working.alleles to the alleles from the population with fewer
  ## alleles, this number becomes nalleles
  if(length(spa.alleles) <= length(spb.alleles)){
    nalleles <- length(spa.alleles)
    working.alleles<-spa.alleles
  }
  else{
    nalleles<- length(spb.alleles)
    working.alleles<-spb.alleles
  }

  if(nalleles == 1){ 
    focalset<-working.alleles
    combo1<-focalset
    combo2<-(1:total.alleles)[-focalset]
    group.a<-c(sum(SpA[focalset]), sum(SpA[-focalset]))
    group.b<-c(sum(SpB[focalset]), sum(SpB[-focalset]))
    max.delta<-delta(group.a, group.b)
  }
  else if (abs(obs.delta - 1) < 0.000001 ){  ## diagnostic alleles from the outset
    focalset<-working.alleles
    combo1<-focalset
    combo2<-(1:total.alleles)[-focalset]
    max.delta<-obs.delta
  }
  else {
    my.levels<-nalleles:1
    for(i in 1:length(my.levels)){
      ## cat(i, "of ", length(my.levels), ", max.delta=", max.delta,", obs.delta=", obs.delta, fill=T)
      if(search.complete==TRUE){
        break	
      }	
      working<-utils::combn(nalleles, my.levels[i])  ## these are now indexes to the working.alleles
      for (x in 1:ncol(working) ){
        focalset<-working.alleles[working[,x]]
        ## note that combn and combinations are transposed relative to
        ## another.  group.a and group.b give the frequency of some
        ## combination of alleles is species a and species b
        group.a<-c(sum(SpA[focalset]), sum(SpA[-focalset]))
        group.b<-c(sum(SpB[focalset]), sum(SpB[-focalset]))
        new.delta<-delta(group.a, group.b)
        if ( abs(new.delta - obs.delta) < 0.00001 ||
            abs(max.delta -obs.delta) < 0.00001 ){
          ## data reduction cant get better than full set
          ## so if this is satisfied, stop
          combo1<-focalset
          combo2<-(1:total.alleles)[-focalset]
          search.complete<-TRUE
          max.delta<- new.delta
          break
        }
        else if(  new.delta > max.delta ){
          max.delta<- new.delta
          combo1<-focalset
          combo2<-(1:total.alleles)[-focalset]
        }	
      }
    }
  }

  class.freq<-matrix(data=c(sum(SpA[combo1]), sum(SpA[combo2]),
                       sum(SpB[combo1]), sum(SpB[combo2])), 
                     dimnames=list(c("spa", "spb"), c("c1", "c2")),
                     byrow=TRUE, ncol=2)
  ## return(list(obs.delta=obs.delta, max.delta=max.delta,
  ## combo1=combo1, combo2=combo2, class.freq=class.freq))
  padding1<-rep(NA,max.alleles-length(combo1))
  padding2<-rep(NA,max.alleles-length(combo2))
  return(c(obs.delta, max.delta, combo1, padding1, combo2, 
           padding2, as.vector(class.freq)))
}

