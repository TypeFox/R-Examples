leafsfirst.bondary<-function(pcf=NULL,lev=NULL,refe=NULL,type="lst",levmet="radius",
ordmet="etaisrec",ngrid=NULL,
dendat=NULL,rho=0,propor=NULL,lowest="dens",f=NULL)
{
# pcf is a piecewise constant object
# type= "lst"/"shape"
# levmet= "radius"/"proba"

if (lowest=="dens") lowest<-0 else lowest<-min(pcf$value)

if ((!is.null(lev)) || (!is.null(propor))){
    type<-"shape"
    if (!is.null(propor)) lev<-propor*max(pcf$value)
    if (is.null(refe)) refe<-locofmax(pcf)
}
if (!is.null(dendat)) type<-"tail"

if (type=="tail"){
   d<-dim(dendat)[2]
   pcf$high<-dendat
   pcf$down<-dendat
   if (is.null(refe)){
       refe<-matrix(0,1,d)
       for (i in 1:d) refe[1,i]<-mean(dendat[,i])
       refe<-refe[1:d]
   }
}
else{
  d<-length(pcf$N)
  step<-matrix(0,d,1)
  for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
}

if (type=="lst"){
  lkm<-length(pcf$value)
  distat<-pcf$value-lowest
  infopointer<-seq(1,lkm)     # links from nodes to recs
}
else if (type=="shape"){
  lenni<-length(pcf$value)
  distat<-matrix(0,lenni,1)
  infopointer<-matrix(0,lenni,1)
  lkm<-0
  for (i in 1:lenni){
    if (pcf$value[i]>=lev){
       lkm<-lkm+1
       nod<-i  #nod<-pcf$nodefinder[i]
       if (ordmet=="etaisrec"){
           recci<-matrix(0,2*d,1)
           for (jj in 1:d){
              recci[2*jj-1]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
              recci[2*jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
           }
           distat[lkm]<-etaisrec(refe,recci)
       }
       else{
          lowi<-matrix(0,d,1)
          uppi<-matrix(0,d,1)
          for (jj in 1:d){
             lowi[jj]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
             uppi[jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
          }
          baryc<-lowi+(uppi-lowi)/2  
          distat[lkm]<-etais(baryc,refe)
       }
       infopointer[lkm]<-i
    }
  }
}
else{  #type=="tail"
   d<-dim(dendat)[2]
   n<-dim(dendat)[1]
   lkm<-dim(dendat)[1]
   distat<-sqrt(pituus(dendat-t(matrix(refe,d,n))))
   infopointer<-seq(1,lkm)
}

distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   
#if (length(rho)==1) rho<-rep(rho,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
proba<-matrix(0,lkm,1)
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

if (type!="tail"){
  # volume calculation
  vol<-1
  k<-1
  ip<-infopointer[node]  #pcf$nodefinder[infopointer[node]]
  while (k<=d){
      vol<-vol*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
      k<-k+1
  }
  volume[node]<-vol
  ip2<-infopointer[node]
  proba[node]<-pcf$value[ip2]*vol

  # ekamome calculation
  newcente<-matrix(0,d,1)
  for (j in 1:d){
    volmin<-1
    k<-1
    while (k<=d){
       if (k!=j){
          volmin<-volmin*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
       }
       k<-k+1
    }
    ala<-pcf$support[2*j-1]+step[j]*pcf$down[ip,j]
    yla<-pcf$support[2*j-1]+step[j]*pcf$high[ip,j]
    newcente[j]<-volmin*(yla^2-ala^2)/2
  }
  ekamome[node,]<-newcente
  distcenter[node,]<-newcente/vol
}
else{  # type=tail
  if (is.null(f)) volume[node]<-1
  else{ 
       ip<-infopointer[node] 
       volume[node]<-1/(f[ip]*length(f))
  }
}

beg<-node             #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-pcf$down[note,i]   
  boundrec[node,2*i]<-pcf$high[note,i]  
}

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

    # radius
    radius[node]<-distat[ord[node]]
    branchradius[node]<-radius[node]
    if (type!="tail"){
       # volume calculation
       vol<-1
       k<-1
       ip<-infopointer[node]    #pcf$nodefinder[infopointer[node]]
       while (k<=d){
          vol<-vol*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
          k<-k+1
       }
       volume[node]<-vol
       ip2<-infopointer[node]
       proba[node]<-pcf$value[ip2]*vol

       # ekamome calculation
       newcente<-matrix(0,d,1)
       for (jj in 1:d){
            volmin<-1
            k<-1
            while (k<=d){
               if (k!=jj){
                   volmin<-volmin*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
               }
               k<-k+1
            }
            ala<-pcf$support[2*jj-1]+step[jj]*pcf$down[ip,jj]
            yla<-pcf$support[2*jj-1]+step[jj]*pcf$high[ip,jj]
            newcente[jj]<-volmin*(yla^2-ala^2)/2
       }
       ekamome[node,]<-newcente
       distcenter[node,]<-newcente/vol
    }
    else{     #type==tail
       if (is.null(f)) volume[node]<-1
       else{ 
          ip<-infopointer[node] 
          volume[node]<-1/(f[ip]*length(f))
       }
    }

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        rhocur<-rho   #rho[infopointer[node]]  
        istouch<-touchstep.boundary(node,curroot,boundrec,child,sibling,
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
           if (type!="tail"){
              volume[node]<-volume[node]+volume[curroot]
              proba[node]<-proba[node]+proba[curroot]
              ekamome[node,]<-ekamome[node,]+ekamome[curroot,]
           }
           else{  # type == tail
              volume[node]<-volume[node]+volume[curroot]
           }

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
# lf is the level set tree or the shape tree
if (type!="tail"){
   for (i in 1:lkm){
      for (j in 1:d){
          ekamome[i,j]<-ekamome[i,j]/volume[i]
      }
   }
   bary<-ekamome[root,]
}
if (type=="shape"){
  maxdis<-sqrt(distat[ord[length(ord)]])
  if (levmet=="proba")
     level<-taillevel(root,#child,sibling,
            parent,volume,proba)
  else 
     level<-sqrt(radius)
}
else{ #type="lst"
     level<-radius+lowest
     maxdis<-distat[ord[length(ord)]]
}
if (type=="tail"){
   center<-t(dendat[infopointer,])
}

if (type!="tail"){
  lf<-list(
  parent=parent,volume=volume,center=t(ekamome),level=level,
  root=root,
  #child=child,sibling=sibling,  #virhe??
  infopointer=infopointer,
  proba=proba,#radius=radius,
  #branchradius=sqrt(branchradius),
  distcenter=t(distcenter),
  refe=refe,maxdis=maxdis,bary=bary,lev=lev)
}
else{
  lf<-list(
  parent=parent,volume=volume,center=center,level=level,
  root=root,
  #child=child,sibling=sibling,  #virhe??
  infopointer=infopointer,
  #proba=proba,#radius=radius,
  #branchradius=sqrt(branchradius),
  #distcenter=t(distcenter),
  refe=refe,maxdis=maxdis,
  dendat=dendat)
}

# if ngrid given, reduce the lst
if (!is.null(ngrid)){
    stepsi<-maxdis/ngrid
    radii<-seq(0,maxdis,stepsi)
    lf<-treedisc(lf,pcf,r=radii,type=type)
}

return(lf)
}





