partition.old<-function(tree,suppo=NULL)
{
#Finds the partition corresponding to an evaluation tree.

#tree is a binary tree, list(val,vec,left,right,...)
#supp is 2*xlkm-vector, support of density

#Result is list(values,recs)
#  values is recnum-vector, values of the estimate
#  recs is recnum*(2*xlkm)-matrix

if (is.null(suppo)) suppo<-tree$support
xlkm<-length(suppo)/2
sollkm<-length(tree$left)           #solmujen lkm

# Find parents and number of leaves:

leafloc<-findleafs(tree$left,tree$right)

N<-matrix(0,sollkm,1)   #number of leaves in the tree whose root is i
p<-matrix(0,sollkm,1)   #parent
t<-sollkm
while (t>=1){
  if ((!is.na(leafloc[t])) && (leafloc[t]==1)){  #l(t)=0 eli ollaan lehdessa
   N[t]<-1
  }
  else if ((!is.na(leafloc[t])) && (leafloc[t]==0)){ #non-terminal node
     p[tree$left[t]]<-t       #p[t+1]<-t
     p[tree$right[t]]<-t      #p[endpoint(tree,t)+1]<-t
     N[t]<-N[tree$left[t]]+N[tree$right[t]]
  }
  t<-t-1
}
                     
leafnum<-N[1]
recs<-matrix(0,leafnum,2*xlkm)
values<-matrix(0,leafnum,1)
frekv<-matrix(0,leafnum,1)

if (leafnum==1){ 
  recs<-suppo
  values<-tree$mean[1]
}
else{
 ind<-1
 for (i in 1:sollkm){
   if ((!is.na(leafloc[i])) && (leafloc[i]==1)){  #i is leaf
     values[ind]<-tree$mean[i]
     recs[ind,]<-suppo
     frekv[ind]<-tree$nelem[i]
     j<-i
     while (p[j]>0){  #we are not in the root
         pare<-p[j]
         vari<-tree$direc[pare]
         split<-tree$split[pare]
         if (tree$left[pare]==j){  #i is left child
           if (split<recs[ind,2*vari]){ #if we have new restriction 
               recs[ind,2*vari]<-split
           }
         }
         else{  #i is right child
           if (split>recs[ind,2*vari-1]){ #if we have new restriction 
               recs[ind,2*vari-1]<-split
           }
         }
         j<-pare
     }
     ind<-ind+1
   }
 }
}
# clean zeros out
finrecs<-matrix(0,leafnum,2*xlkm)
finvalues<-matrix(0,leafnum,1)
finfrekv<-matrix(0,leafnum,1)
notzero<-0
for (i in 1:leafnum){
  if (values[i]>0){
      notzero<-notzero+1
      finvalues[notzero]<-values[i]
      finrecs[notzero,]<-recs[i,]
      finfrekv[notzero]<-frekv[i]
  }
}
finvalues<-finvalues[1:notzero]
finrecs<-finrecs[1:notzero,]
finfrekv<-finfrekv[1:notzero]

return(list(values=finvalues,recs=finrecs,frekv=finfrekv))
}













