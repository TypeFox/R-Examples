leafsfirst.lst<-function(pcf, ngrid=NULL, predictor=NULL, type=NULL)
{
rho<-0

d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]

lkm<-length(pcf$value)
distat<-pcf$value
infopointer<-seq(1,lkm)     # links from nodes to recs

distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   
if (length(rho)==1) rho<-rep(rho,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
ekamome<-matrix(0,lkm,d)
distcenter<-matrix(0,lkm,d)
branchradius<-matrix(0,lkm,1)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]
branchradius[node]<-radius[node]

volume[node]<-1

beg<-node             #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-pcf$down[note,i]   
  boundrec[node,2*i]<-pcf$high[note,i]  
}

found.predictor.node<-FALSE
if ((!is.null(predictor))&&(!found.predictor.node)){
   predictor.rec<-matrix(0,2*d,1)
   for (ii in 1:d){ 
     predictor.rec[2*ii-1]<-floor((predictor[ii]-pcf$support[2*ii-1])/step[ii])
     predictor.rec[2*ii]<-ceiling((predictor[ii]-pcf$support[2*ii-1])/step[ii])
   }
   if (touch(predictor.rec,boundrec[node,])) predictor.node<-node
}
else predictor.node<-NULL

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)  #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-pcf$down[note,i]  
         rec1[2*i]<-pcf$high[note,i] 
    }
    boundrec[node,]<-rec1

    if ((!is.null(predictor))&&(!found.predictor.node)){
       if (touch(predictor.rec,boundrec[node,])) predictor.node<-node
    }

    # radius
    radius[node]<-distat[ord[node]]
    branchradius[node]<-radius[node]
    
    volume[node]<-1

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        rhocur<-rho[infopointer[node]]  
        istouch<-touchstep(node,curroot,boundrec,child,sibling,
                           infopointer,pcf$down,pcf$high,rhocur)
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

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])
           if (branchradius[node]<=branchradius[curroot]) 
                  distcenter[node,]<-distcenter[curroot,]
           branchradius[node]<-max(branchradius[node],branchradius[curroot])

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
          ekamome[i,j]<-ekamome[i,j]/volume[i]
      }
}
bary<-ekamome[root,]

level<-radius
maxdis<-distat[ord[length(ord)]]

lf<-list(
  parent=parent,volume=volume,center=t(ekamome),level=level,
  root=root,
  infopointer=infopointer,
  distcenter=t(distcenter),
  maxdis=maxdis,bary=bary,predictor.node=predictor.node)


# if ngrid given, reduce the lst
if (!is.null(ngrid)){
    stepsi<-maxdis/ngrid
    radii<-seq(0,maxdis,stepsi)
    lf<-treedisc(lf,pcf,r=radii,type=type)
}

return(lf)
}





