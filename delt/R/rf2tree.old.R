rf2tree.old<-function(forest,suppo)
{
d<-length(suppo)/2

nr<-length(forest$ndbigtree)     #number of trees in the forest

nrtreemap<-length(forest$treemap)
map<-matrix(0,nrtreemap,1)
infopointer<-matrix(0,nrtreemap,1)
rootinfo<-matrix(0,nr,1)

# create infopointer
laskuri<-1
for (ij in 1:nrtreemap){
     if (forest$treemap[ij]==2) laskuri<-laskuri+1
     if (forest$treemap[ij]!=0){
            infopointer[ij]<-laskuri
            laskuri<-laskuri+1
     }
}
# create rootinfo
rootinfo[1]<-1
cusu<-0
ii<-1
while (ii<=nr){
   rootinfo[ii]<-cusu+1
   cusu<-cusu+forest$ndbigtree[ii]
   ii<-ii+1
}

totalrunner<-0
glob<-1

#####################################################################
while (glob <= nr){

# build a small tree

maxnrnodes<-2*forest$ndbigtree[glob]

left<-matrix(0,maxnrnodes,1)
right<-matrix(0,maxnrnodes,1)
val<-matrix(0,maxnrnodes,1)
vec<-matrix(0,maxnrnodes,1)
mean<-matrix(0,maxnrnodes,1)
nelem<-matrix(0,maxnrnodes,1)
low<-matrix(0,maxnrnodes,d)
upp<-matrix(0,maxnrnodes,d)

# create root

{
if (forest$treemap[1]!=0){
   #left[1]<-2
   #right[1]<-3

   locu<-rootinfo[glob]

   val[1]<-forest$upper[locu]
   vec[1]<-forest$mbest[locu]
   mean[1]<-forest$avnode[locu]
   for (si in 1:d){
      low[1,si]<-suppo[2*si-1]
      upp[1,si]<-suppo[2*si]
   }

   #map[1]<-2
   #map[2]<-3

   nodesleft<-1
   treeind<-1    #3
   
   prevbeg<-totalrunner #0 #1  #beg of previous level
   prevend<-totalrunner #0 #2  #end of previous level
   curbeg<-totalrunner+1  #1  #3
   curend<-totalrunner+2  #2  #6

   #totalrunner<-2
}
else{
   nodesleft<-0
   treeind<-1
}
}

while (nodesleft==1){

   nodesleft<-0
   prevlkm<-0
   curlkm<-curend-curbeg+1
   muisti<-matrix(0,curlkm,1)

   i<-curbeg

   while (i <= (curbeg+(curend-curbeg+1)/2-1)){
       indl<-curbeg+(2*(i-curbeg+1)-1)-1
       indr<-curbeg+2*(i-curbeg+1)-1
       #ind<-i-curbeg+1
       #parind<-prevbeg+ind-1

       if (forest$treemap[indl]!=0){
           
           if (prevbeg==prevend) curparent<-1 
           else{ 
                  parind<-vanh[indl-curbeg+1]    
                  parindmap<-prevbeg+parind-1           
                  curparent<-map[parindmap]
           }

           node<-treeind+1
           left[curparent]<-node

           locu<-infopointer[indl]

           val[node]<-forest$upper[locu]
           vec[node]<-forest$mbest[locu]
           mean[node]<-forest$avnode[locu]
           for (si in 1:d){
              low[node,si]<-low[curparent,si]
              upp[node,si]<-upp[curparent,si]
           }
           split<-vec[curparent]
           upp[node,split]<-val[curparent]

           map[indl]<-node
           treeind<-treeind+1
           prevlkm<-prevlkm+1
           nodesleft<-1
       
        #if (forest$treemap[indr]!=0)

           node<-treeind+1
           right[curparent]<-node

           locu<-infopointer[indr]

           val[node]<-forest$upper[locu]
           vec[node]<-forest$mbest[locu]
           mean[node]<-forest$avnode[locu]
           for (si in 1:d){
              low[node,si]<-low[curparent,si]
              upp[node,si]<-upp[curparent,si]
           }
           split<-vec[curparent]
           low[node,split]<-val[curparent]

           map[indr]<-node 
           treeind<-treeind+1
           prevlkm<-prevlkm+1
           nodesleft<-1
       }
       i<-i+1
   }

   prevbeg<-curbeg
   prevend<-curend
   curbeg<-curend+1
   curend<-curbeg+2*prevlkm-1

   vanh<-matrix(0,curend-curbeg+1,1)
   liuk<-0
   sep<-1
   while (sep <= (prevend-prevbeg+1)){
        if (forest$treemap[prevbeg+sep-1]!=0){
             liuk<-liuk+1
             vanh[2*liuk-1]<-sep
             vanh[2*liuk]<-sep
        }
        sep<-sep+1
   }

}

while ((forest$treemap[curend]==0) && (glob<nr)) curend<-curend+1
totalrunner<-curend-1

left<-left[1:treeind]
right<-right[1:treeind]
val<-val[1:treeind]
vec<-vec[1:treeind]
mean<-mean[1:treeind]
nelem<-nelem[1:treeind]
low<-low[1:treeind,]
upp<-upp[1:treeind,]

trnew<-list(val=val,vec=vec,mean=mean,nelem=nelem,
left=left,right=right,low=low,upp=upp)

# small tree ready
{
if (glob>1) tr<-treeadd(tr,trnew,d)
else tr<-trnew
}

glob<-glob+1

}

return(tr)
}









