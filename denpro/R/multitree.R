multitree<-function(parents)
{
#Makes sibling-links and child-links

itemnum<-length(parents)
sibling<-matrix(0,itemnum,1)
child<-matrix(0,itemnum,1)
roots<-matrix(0,itemnum,1)
siborder<-matrix(0,itemnum,1)

rootnum<-0
for (i in itemnum:1){
  par<-parents[i]
  if (par==0){   #i is root (does not have parent)
     rootnum<-rootnum+1
     roots[rootnum]<-i
     siborder[i]<-rootnum
  } 
  else{          #i has parent
      if (child[par]==0){  #no childs so far
        child[par]<-i
        siborder[i]<-1
      }
      else{    #go to the end of sibling list
        chi<-child[par]
        sibsi<-1
        while(sibling[chi]>0){
           chi<-sibling[chi]
           sibsi<-sibsi+1
        }
        sibling[chi]<-i    #put to the sibling list
        siborder[i]<-sibsi+1
      }
  }
}
roots<-roots[1:rootnum]
return(list(child=child,sibling=sibling,roots=roots,siborder=siborder))
}
