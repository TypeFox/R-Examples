modegraph<-function(estiseq,hseq=NULL,paletti=NULL)  #,reverse=F)
{
# we want that the largest h is first (1 mode, oversmoothing)

if (is.null(hseq))
   if (!is.null(estiseq$type)){
       if (estiseq$type=="greedy") hseq<--estiseq$hseq
       if (estiseq$type=="bagghisto") hseq<--estiseq$hseq
       if (estiseq$type=="carthisto")  hseq<--estiseq$leaf
       if (estiseq$type=="kernel")  hseq<-estiseq$hseq    
   }
   else hseq<-estiseq$hseq

hnum<-length(hseq)

treelist<-estiseq$lstseq
d<-dim(treelist[[1]]$center)[1]

if (hseq[1]<hseq[2]){   #(reverse){  
    #if ((hnum>1) && (is.null(hseq))) 
    hseq<-hseq[seq(hnum,1)]
    apuseq<-list(treelist[[hnum]])
    i<-2
    while (i <= hnum){
         apuseq<-c(apuseq,list(treelist[[hnum-i+1]]))
         i<-i+1 
   }
   treelist<-apuseq
}

if (is.null(paletti))
paletti<-c("red","blue","green","turquoise","orange","navy",
"darkgreen","orchid",colors()[50:100])

low<-matrix(0,hnum,1)
upp<-matrix(0,hnum,1)
tot<-moodilkm(treelist[[1]]$parent)$lkm  #tot is the number of modes over all lst:s
low[1]<-1
upp[1]<-tot
i<-2
while (i <= hnum){
  lkmm<-moodilkm(treelist[[i]]$parent)$lkm
  tot<-tot+lkmm
  low[i]<-upp[i-1]+1
  upp[i]<-low[i]+lkmm-1
  i<-i+1
}

xcoor<-matrix(0,tot,d)
ycoor<-matrix(0,tot,1)
parent<-matrix(0,tot,1)
mlabel<-matrix(0,tot,1)
flucpoints<-matrix(0,hnum,1)
nodepointer<-matrix(0,tot,1)
colot<-matrix("",tot,1)

# first we allocate colors for the largest h
colrun<-1  #low[1]
while (colrun<=upp[1]){
   colot[colrun]<-paletti[colrun]
   colrun<-colrun+1
}

laskuri<-1
srun<-1
mlkmpre<-1
flucnum<-0
while (srun<=hnum){  
    mlkm<-moodilkm(treelist[[srun]]$parent)
    if (mlkmpre < mlkm$lkm){
          flucnum<-flucnum+1
          flucpoints[flucnum]<-srun
    }

    for (j in 1:mlkm$lkm){
        loca<-mlkm$modloc[j]
        if (d>1){
           for (dim in 1:d){
              xcoor[laskuri,dim]<-treelist[[srun]]$center[dim,loca]
           }
        }
        else{
              xcoor[laskuri]<-treelist[[srun]]$center[loca]
        }
        ycoor[laskuri]<-hseq[srun]
        mlabel[laskuri]<-j
        nodepointer[laskuri]<-loca

        laskuri<-laskuri+1
    }

    if (srun>1){

       vec1<-matrix(0,mlkmpre,d)
       vec2<-matrix(0,mlkm$lkm,d)
       vec1[1:mlkmpre,]<-xcoor[low[srun-1]:upp[srun-1],]
       vec2[1:mlkm$lkm,]<-xcoor[low[srun]:upp[srun],]
       vm<-vectomatch(vec1,vec2)
       for (jj in low[srun]:upp[srun]){
           parent[jj]<-vm$parent[jj-low[srun]+1]+low[srun-1]-1
           if (vm$newnode[jj-low[srun]+1]==1){ 
                colot[jj]<-paletti[colrun]
                colrun<-colrun+1
           }
           else colot[jj]<-colot[parent[jj]]
      }
    }

    mlkmpre<-mlkm$lkm
    srun<-srun+1 
}

xcoor<-xcoor[1:(laskuri-1),]
ycoor<-ycoor[1:(laskuri-1)]
parent<-parent[1:(laskuri-1)]
colot<-colot[1:(laskuri-1)]
mlabel<-mlabel[1:(laskuri-1)]
nodepointer<-nodepointer[1:(laskuri-1)]
flucpoints<-flucpoints[1:flucnum]

mt<-list(xcoor=xcoor,ycoor=t(ycoor),
parent=parent,colot=colot,hseq=hseq,type=estiseq$type,
upp=upp,low=low,
mlabel=t(mlabel),
flucpoints=t(flucpoints),
nodepointer=t(nodepointer)
)

return(mt)
}




