shapetree<-function(et,lev,bary,ordmet="etaisrec",levmet="proba")
{
# et is an evaluation tree

d<-length(et$step)

# order the atoms for the level set with level "lev"

lenni<-length(et$value)
distat<-matrix(0,lenni,1)
infopointer<-matrix(0,lenni,1)
lkm<-0
for (i in 1:lenni){
  if (et$value[i]>=lev){
     lkm<-lkm+1
     nod<-i  #nod<-et$nodefinder[i]
     if (ordmet=="etaisrec"){
         recci<-matrix(0,2*d,1)
         for (jj in 1:d){
            recci[2*jj-1]<-et$support[2*jj-1]+et$step[jj]*et$low[nod,jj]
            recci[2*jj]<-et$support[2*jj-1]+et$step[jj]*et$upp[nod,jj]
         }
         distat[lkm]<-etaisrec(bary,recci)
     }
     else{
        lowi<-matrix(0,d,1)
        uppi<-matrix(0,d,1)
        for (jj in 1:d){
            lowi[jj]<-et$support[2*jj-1]+et$step[jj]*et$low[nod,jj]
            uppi[jj]<-et$support[2*jj-1]+et$step[jj]*et$upp[nod,jj]
        }
        baryc<-lowi+(uppi-lowi)/2  #et$low[nod,]+(et$upp[nod,]-et$low[nod,])/2  
        distat[lkm]<-etais(baryc,bary)
     }
     infopointer[lkm]<-i
  }
}
distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   #pointe->et$value,et$nodefinder

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
center<-matrix(0,lkm,d)
radius<-matrix(0,lkm,1)

proba<-matrix(0,lkm,1)
ekamome<-matrix(0,lkm,d)
highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0
# volume calculation
vol<-1
k<-1
ip<-infopointer[node]  #et$nodefinder[infopointer[node]]
while (k<=d){
    vol<-vol*(et$upp[ip,k]-et$low[ip,k])*et$step[k]
    k<-k+1
}
volume[node]<-vol
ip2<-infopointer[node]
proba[node]<-et$value[ip2]*vol
radius[node]<-distat[ord[node]]
# ekamome calculation
newcente<-matrix(0,d,1)
for (j in 1:d){
  volmin<-1
  k<-1
  while (k<=d){
      if (k!=j){
         volmin<-volmin*(et$upp[ip,k]-et$low[ip,k])*et$step[k]
      }
      k<-k+1
  }
  ala<-et$support[2*j-1]+et$step[j]*et$low[ip,j]
  yla<-et$support[2*j-1]+et$step[j]*et$upp[ip,j]
  newcente[j]<-volmin*(yla^2-ala^2)/2
}
ekamome[node,]<-newcente

beg<-node             #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-et$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-et$low[note,i]   
  boundrec[node,2*i]<-et$upp[note,i]  #et$index[infopointer[node],i]
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)  #luo sigleton
    note<-infopointer[node]  #note<-et$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-et$low[note,i]  
         rec1[2*i]<-et$upp[note,i] 
    }
    boundrec[node,]<-rec1

    # volume calculation
    vol<-1
    k<-1
    ip<-infopointer[node]    #et$nodefinder[infopointer[node]]
    while (k<=d){
          vol<-vol*(et$upp[ip,k]-et$low[ip,k])*et$step[k]
          k<-k+1
    }
    volume[node]<-vol
    ip2<-infopointer[node]
    proba[node]<-et$value[ip2]*vol
    radius[node]<-distat[ord[node]]
    # ekamome calculation
    newcente<-matrix(0,d,1)
    for (jj in 1:d){
         volmin<-1
         k<-1
         while (k<=d){
            if (k!=jj){
                volmin<-volmin*(et$upp[ip,k]-et$low[ip,k])*et$step[k]
            }
            k<-k+1
         }
         ala<-et$support[2*jj-1]+et$step[jj]*et$low[ip,jj]
         yla<-et$support[2*jj-1]+et$step[jj]*et$upp[ip,jj]
         newcente[jj]<-volmin*(yla^2-ala^2)/2
    }
    ekamome[node,]<-newcente

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        istouch<-touchstep(node,curroot,boundrec,child,sibling,
                           infopointer,et$low,et$upp)
        if (istouch==1){
{
           # paivita parent, child, sibling, volume ekamome
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }
           volume[node]<-volume[node]+volume[curroot]
           proba[node]<-proba[node]+proba[curroot]
           ekamome[node,]<-ekamome[node,]+ekamome[curroot,]
           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]
}
        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter
for (i in 1:lkm){
   for (j in 1:d){
       center[i,j]<-ekamome[i,j]/volume[i]
   }
}

if (levmet=="proba")
level<-taillevel(root,#child,sibling,
parent,volume,proba)
else level<-radius

return(list(
parent=parent,volume=volume,center=t(center),level=level,
root=root,
#child=child,sibling=sibling,  #virhe??
infopointer=infopointer,
proba=proba,radius=radius,
bary=bary,maxdis=distat[ord[length(ord)]]))

}





