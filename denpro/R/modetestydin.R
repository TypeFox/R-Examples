modetestydin<-function(lstree,h,N,Q,bootnum,delta,nsimu,minim,kernel){

#we will approximate the estimate with a function whose values are
#levels of level set tree, this estimate has already been 
#normalized to integrate to one

index<-lstree$index
etnum<-dim(index)[1]
d<-length(N)

len<-length(lstree$parent)
mut<-multitree(lstree$parent)
mt<-pruneprof(mut)

#branchnodes<-findbranchMT(mt,len)
branchn<-findbranch(lstree$parent)$indicator
bnumbeg<-length(branchn)
branchnodes<-matrix(0,bnumbeg,1)
bnum<-0
for (i in 1:bnumbeg){
  if (branchn[i]==1){
      bnum<-bnum+1
      branchnodes[bnum]<-i
  }
}
branchnodes<-branchnodes[1:bnum]

testvec<-matrix(0,len,1)    #this is output

i<-1
while (i<=bnum){
   brnode<-branchnodes[i]
   
   #first cut the level set tree
   #3 cases: 1. brnode is linked from parent
   #         2. brnode is linked from previous sibling

   newchild<-mut$child
   newsibling<-mut$sibling
   newroots<-mut$roots

   brpare<-lstree$parent[brnode]
  
   if (brpare>0){  # brnode is not a root
 
   etsi<-mut$child[brpare]
{  if (etsi==brnode){
       newchild[brpare]<-mut$sibling[etsi]
   }
   else{
      while (etsi!=brnode){
         prevetsi<-etsi
         etsi<-mut$sibling[etsi]
      }
      newsibling[prevetsi]<-mut$sibling[etsi]
   }
} 
  
   #normalize the estimate to integrate to one
   
   kerroin<-cintemul(newroots,newchild,newsibling,lstree$volume,lstree$level)
   newlevel<-lstree$level/kerroin
   
   #creat value-vector "newvalue" with cutted values
   #newvalue*volofatom is probablility vector
   
   newvalue<-cutvalue(newroots,newchild,newsibling,
             newlevel,lstree$component,
             lstree$AtomlistAtom,lstree$AtomlistNext,etnum)
   
   #calculate the test statistics

         tstat<-excmas(lstree)[brnode]
          # cuttedlevel<-lstree$level[brpare]
          #cintemul(brnode,mut$child,mut$sibling,
          #lstree$volume,lstree$level-cuttedlevel)
   
   #generate samples from the cutted estimate
   
   overfrekv<-0
   j<-1
   while (j<=bootnum){
      dendatj<-simukern(nsimu,d,seed=j,newvalue,index,delta,minim,h)
      pk<-profkern(dendatj,h,N,Q,compoinfo=T,kernel=kernel)

      #find modes which are in the set associated with node brnode

      mlkm<-moodilkm(pk$parent)
      setissa<-matrix(0,mlkm$lkm,1)
      setissalkm<-0
      for (mrun in 1:mlkm$lkm){
           mloc<-mlkm$modloc[mrun]
           kandi<-pk$center[,mloc]
           torf<-onsetissa(kandi,h,delta,minim,
                           brnode,lstree$component,
                           lstree$index,
                           lstree$AtomlistAtom,lstree$AtomlistNext)
           if (torf){
               setissalkm<-setissalkm+1
               setissa[setissalkm]<-mloc 
           }
      }
      if (setissalkm>0){ 
        
         setissa<-setissa[1:setissalkm]

      #calculate the excess mass over "cuttedlevel" for
      #the modes in setissa 
      #Note that there may be a branching after "cuttedlevel"
      #We find the multitree for pk, then we cut this tree
      #by choosing the root node to be those nodes which are
      #arrived from modes (stop when the cuttedlevel is passed)
      
      cuttedlevel<-lstree$level[brpare]
      pkmut<-multitree(pk$parent)
      pkroots<-matrix(0,setissalkm,1)
      pkrootslkm<-0
      for (mrun in 1:setissalkm){
         node<-setissa[mrun]
         while ((pk$level[node]>cuttedlevel) && (node>0)){
              node<-pk$parent[node]
         }
         if (node>0){
             pkrootslkm<-pkrootslkm+1
             pkroots[pkrootslkm]<-node
         }
      }
      if (pkrootslkm>0){
         pkroots<-pkroots[1:pkrootslkm]

         bootstat<-cintemul(pkroots,pkmut$child,pkmut$sibling,
                            pk$volume,pk$level-cuttedlevel)
      }
      else{    #setissalkm>0, pkrootslkm==0
         bootstat<-0
      }
      
      }
      else{    #setissalkm==0
          bootstat<-0
      }

      if (bootstat<tstat){
          overfrekv<-overfrekv+1
      }
      j<-j+1
   }
   testvec[brnode]<-overfrekv/bootnum    #p-values are returned

   }

   i<-i+1
}

return(testvec)
}













