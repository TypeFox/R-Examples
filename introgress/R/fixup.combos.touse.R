fixup.combos.touse <-
function(working, myset, nallele){
  for(i in 1:length(myset)){
    if(working$spb.c1[i] < working$spa.c1[i]){
                                        # swap c1 and c2
      new.spa.c1<-working$spa.c2[i]
      new.spb.c1<-working$spb.c2[i]
      new.spa.c2<-working$spa.c1[i]
      new.spb.c2<-working$spb.c1[i]
      working$spa.c1[i] <- new.spa.c1
      working$spb.c1[i] <- new.spb.c1
      working$spa.c2[i] <- new.spa.c2
      working$spb.c2[i] <- new.spb.c2

      spb.range<-(3+nallele):(3+nallele*2-1)  ## the 3 refers to the third column onward in working
      spa.range<-3:(3+nallele-1)
      new.c1<-working[i,spb.range]
      new.c2<-working[i,spa.range]
      working[i,spb.range]<-new.c2
      working[i,spa.range]<-new.c1
    }
  }
  return(working)
}

