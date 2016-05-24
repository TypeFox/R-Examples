joinconne<-function(leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index){
#
#
# 1. We find the startpoints for the sets
#
# we will make three vectors of startpoints:
# a. startpoints for the whole sets (concatenate left and right child)
#    (startpointsS)
# b. startpoints for the right boundary of the left child and 
#    left boundary of the right child
#    (startpointsB)
# c.1. startpoints for the left boundary of the left child a
#    (startpointsNewBleft)
# c.2. startpoints for the right boundary of the right child
#    (startpointsNewBright)
#
#
# startpoints are pointers to 
# a. atomsSepaAtoms/atomsSepaNext
# b. atomsSepaAtoms/atomsRBounNext, atomsSepaAtoms/atomsLBounNext
# c.1. atomsSepaAtoms/atomsLBounNext
# c.2.  atomsSepaAtoms/atomsRBounNext
#
# startpoints are found from 
# a. begsSepaBegs
# b. begsLeftBoun, begsRighBoun
# c.1 begsLeftBoun, 
# c.2 begsRighBoun
#
# b. is used to check which touch
# a., c.1. and c.2. are joined together
# Note that some sets of the boundary are empty 
# (we store 0 to the respective location in begsLeftBoun, begsRighBoun )
#
suppi<-length(begsSepaNext)
startpointsS<-matrix(0,suppi,1)
startpointsB<-matrix(0,suppi,1)
startpointsNewBleft<-matrix(0,suppi,1)
startpointsNewBright<-matrix(0,suppi,1)
#new boundary: left bound. of left child, right b. of right child
#
induksi<-1
anfang<-separy[leftbeg]
startpointsS[induksi]<-anfang
startpointsB[induksi]<-begsRighBoun[anfang]
startpointsNewBleft[induksi]<-begsLeftBoun[anfang]
while (begsSepaNext[anfang]>0){
  anfang<-begsSepaNext[anfang]
  induksi<-induksi+1
  startpointsS[induksi]<-begsSepaBegs[anfang]
  startpointsB[induksi]<-begsRighBoun[anfang]
  startpointsNewBleft[induksi]<-begsLeftBoun[anfang]  
}
mleft<-induksi
induksi<-induksi+1
anfang<-separy[rightbeg]
startpointsS[induksi]<-anfang
startpointsB[induksi]<-begsLeftBoun[anfang]
startpointsNewBright[induksi]<-begsRighBoun[anfang]
while (begsSepaNext[anfang]>0){
  anfang<-begsSepaNext[anfang]
  induksi<-induksi+1
  startpointsS[induksi]<-begsSepaBegs[anfang]
  startpointsB[induksi]<-begsLeftBoun[anfang]
  startpointsNewBright[induksi]<-begsRighBoun[anfang]  

}
startpointsS<-startpointsS[1:induksi]
startpointsB<-startpointsB[1:induksi]
startpointsNewBleft<-startpointsNewBleft[1:induksi]
startpointsNewBright<-startpointsNewBright[1:induksi]
m<-induksi
mright<-m-mleft  
#
#
# 2. We make "links" matrix and apply declev
#
# We utilize previous programs
#
linkit<-matrix(0,m,m)
do<-1
while (do <= mleft){
   beg1<-startpointsB[do]    #could be 0
   re<-mleft+1
   while (re <= m){
       beg2<-startpointsB[re]    #could be 0
       conne<-FALSE
       begbeg1<-beg1
       while (begbeg1>0){
            begbeg2<-beg2
            while (begbeg2>0){
                atom1<-atomsSepaAtom[begbeg1]
                inde1<-index[atom1,]
                atom2<-atomsSepaAtom[begbeg2]
                inde2<-index[atom2,]
                if (dotouch(inde1,inde2,direction)){
                   conne<-TRUE
                }
                begbeg2<-atomsLBounNext[begbeg2]
            }
            begbeg1<-atomsRBounNext[begbeg1]
       }                
       if (conne){
           linkit[do,re]<-1
       }
       re<-re+1
   }
   do<-do+1
}
for (do in (mleft+1):m){
   beg1<-startpointsB[do]
   for (re in 1:mleft){
       beg2<-startpointsB[re]
       conne<-FALSE
       begbeg1<-beg1
       while (begbeg1>0){
            begbeg2<-beg2
            while (begbeg2>0){
                atom1<-atomsSepaAtom[begbeg1]
                inde1<-index[atom1,]
                atom2<-atomsSepaAtom[begbeg2]
                inde2<-index[atom2,]
                if (dotouch(inde1,inde2,direction)){
                   conne<-TRUE
                }
                begbeg2<-atomsRBounNext[begbeg2]
            }
            begbeg1<-atomsLBounNext[begbeg1]
       }                
       if (conne){
           linkit[do,re]<-1
      }
   }
} 
# huom ylla on nopeutettu, koska tiedetaan, etta atomit
# 1,...,mleft eivat koske toisiaan ja samoin atomit mleft+1,...,m
#
# next we apply "declev" 
rindeksitB<-seq(1,m)
res<-declevnew(rindeksitB,linkit,m)   #res is sepnum*m-matrix of 0/1
sepnum<-dim(res)[1]
# 
# res is sepnum*m-matrix, 1 in some row indicates that set (atom)
# belongs to this component, 0 in other positions
#
# output could be also a vector which contains pointers
# to a list of elements (in one list we find those sets which
# belong to the same component
#
#compopointer<-matrix(0,sepnum,1) 
#compoSet<-matrix(0,m,1)
#compoNext<-matrix(0,m,1)
#
#
#3. We join the sets 
#
# We join the sets whose startpoints are in 
# startpointsS and startpointsNewBleft, startpointsNewBright
# We have pointers separy[leftbeg] and separy[rightbeg]
# which contain pointers to lists which we can utilize
# to make a new list (these two lists contain together at most as many 
# elements as we need)
# 
# cut first list or (join these two lists and cut second)
#
TotalBeg<-separy[leftbeg]
#
tavoite<-1
hiihtaja<-TotalBeg
while ((begsSepaNext[hiihtaja]>0) && (tavoite<sepnum)){
   hiihtaja<-begsSepaNext[hiihtaja]
   tavoite<-tavoite+1
}  
if (tavoite<sepnum){  #now hiihtaja points to the end of the first list
   #join the lists
   begsSepaNext[hiihtaja]<-separy[rightbeg]
   #we continue
   hiihtaja<-separy[rightbeg]
   tavoite<-tavoite+1
   while ((begsSepaNext[hiihtaja]>0) && (tavoite<sepnum)){
      hiihtaja<-begsSepaNext[hiihtaja]
      tavoite<-tavoite+1
   }    
   begsSepaNext[hiihtaja]<-0
}
else{  #we have reached goal, cut without joining
   begsSepaNext[hiihtaja]<-0
}
#
#
nykyinen<-TotalBeg
i<-1
while (i<= sepnum){
  len<-sum(res[i,])            # number of sets to be joined
  #
  # we find vectors which contain pointer to the beginnings
  # of lists of atoms
  #
  osoittajaS<-matrix(0,len,1)  #make vector of pointers to the begs of sets
  osoittajaNewBleft<-matrix(0,len,1)
  osoittajaNewBright<-matrix(0,len,1)
  laskuri<-1
  for (j in 1:m){
     if (res[i,j]==1){
          osoittajaS[laskuri]<-startpointsS[j]  
          osoittajaNewBleft[laskuri]<-startpointsNewBleft[j]    #could be 0
          osoittajaNewBright[laskuri]<-startpointsNewBright[j]  #could be 0
          laskuri<-laskuri+1
     }    
  }
  #
  # handle separy 
  #
  begsSepaBegs[nykyinen]<-osoittajaS[1]    #always non-zero
  #
  k<-1
  while (k<=(len-1)){    
      curre<-osoittajaS[k]
      while (atomsSepaNext[curre]>0){    #find the end
          curre<-atomsSepaNext[curre]
      }
      atomsSepaNext[curre]<-osoittajaS[k+1]
      k<-k+1
  }
  #
  # handle left boundary
  #
  # set kL=0 if all 0 , otherwise kL first nonzero
  #
  k<-1
  while ((k<=len) && (osoittajaNewBleft[k]==0)){
      k<-k+1
  }
  if (k>len){   # all zero
     kL<-0
     begsLeftBoun[nykyinen]<-0
  }
  else{         # kL is first non zero
     kL<-k
     begsLeftBoun[nykyinen]<-osoittajaNewBleft[kL]  
  #
  # update the list of left boundaries
  # concatenate the lists of atoms
  #
  k<-kL
  while (k<=(len-1)){    
      curre<-osoittajaNewBleft[k]         # curre is not zero
      while (atomsLBounNext[curre]>0){    #find the end
            curre<-atomsLBounNext[curre]
      }
      # find the next non zero
      k<-k+1
      while ((k<=len) && (osoittajaNewBleft[k]==0)){
           k<-k+1
      }
      if (k>len){
          atomsLBounNext[curre]<-0
      }
      else{  # found nonzero
          atomsLBounNext[curre]<-osoittajaNewBleft[k]
      }
  }
  }
  #
  # handle right boundary
  #
  # set kR=0 if all 0 , otherwise kR first nonzero
  #
  k<-1
  while ((k<=len) && (osoittajaNewBright[k]==0)){
      k<-k+1
  }
  if (k>len){
     kR<-0
     begsRighBoun[nykyinen]<-0
  }
  else{
     kR<-k
     begsRighBoun[nykyinen]<-osoittajaNewBright[kR]  
  #
  # update the list of right boundaries
  # concatenate the lists of atoms
  #
  k<-kR
  while (k<=(len-1)){    
      curre<-osoittajaNewBright[k]         # curre is not zero
      while (atomsRBounNext[curre]>0){    #find the end
            curre<-atomsRBounNext[curre]
      }
      # find the next non zero
      k<-k+1
      while ((k<=len) && (osoittajaNewBright[k]==0)){
           k<-k+1
      }
      if (k>len){
          atomsRBounNext[curre]<-0
      }
      else{  # found nonzero
          atomsRBounNext[curre]<-osoittajaNewBright[k]
      }
  }
  }
  #
  # we move to the next sepaset
  nykyinen<-begsSepaNext[nykyinen]
  i<-i+1
}
#
return(list(totbegSepary=TotalBeg,separy=separy,
begsSepaNext=begsSepaNext,begsSepaBegs=begsSepaBegs,
begsLeftBoun=begsLeftBoun,begsRighBoun=begsRighBoun,
atomsSepaNext=atomsSepaNext,atomsSepaAtom=atomsSepaAtom,
atomsLBounNext=atomsLBounNext,atomsRBounNext=atomsRBounNext))
}























