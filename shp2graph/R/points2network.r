#mapping.method:
#1 means "Mapped to the nearest node"
#2 means "Mapped to the foot point (as an new node) of the nearest edge"
#3 means "Add virtual edges"
points2network<-function(ntdata,pointsxy,mapping.method=1,ELComputed=FALSE, longlat=F,Detailed=F, ea.prop=NULL)
{
  VElist<-NULL
  if(!is(ntdata,"SpatialLinesDataFrame"))
     stop("Input data is not a proper spatial network data frame, here only SpatialLinesDataFrame is accepted.")
  if (mapping.method%%1!=0)
     stop("Only integer 1,2,3 is accepted to specify a mapping.method.")
  else
      if (mapping.method>4 ||mapping.method<1)
      stop("No matched mapping.method available.")
  if (mapping.method==1)
  {
     res<-readshpnw(data=ntdata, Detailed=Detailed, ELComputed=ELComputed)
     Nodeslist<-res[[2]]
     Edgeslist<-res[[3]]
     edgelength<-res[[4]]
     Nodexlist<-res[[6]]
     Nodeylist<-res[[7]]
     Eadf<-res[[5]]
     CoorespondIDs<-Nearest.nodes(pointsxy, Nodeslist,longlat=longlat)
  }
  else if (mapping.method==2)
  {
    res<-footpoint.nodes(ntdata=ntdata,pointsxy=pointsxy, longlat=longlat,ea.prop=ea.prop, ELComputed=ELComputed)
    Nodeslist<-res[[1]]
    CoorespondIDs<-res[[6]]
    Edgeslist<-res[[4]]
    Nodexlist<-res[[2]]
    Nodeylist<-res[[3]]                    
    Eadf<-res[[5]]
    edgelength<-res[[7]]
  }
  else if (mapping.method==3)
  {
     res<-readshpnw(data=ntdata, ELComputed=T, Detailed=Detailed)
     Nodeslist<-res[[2]]
     Edgeslist<-res[[3]]
     Nodexlist<-res[[6]]
     Nodeylist<-res[[7]]
     Eadf<-res[[5]]
     edgelength<-res[[4]]
     CoorespondIDs<-Nearest.nodes(pointsxy, Nodeslist,longlat=longlat)
     res<-virtualedge.nn(Nodeslist, Edgeslist, Nodexlist, Nodeylist, pointsxy, edgelength, Eadf,CoorespondIDs,longlat=longlat,ea.prop)
     Nodeslist<-res[[1]]
     Edgeslist<-res[[2]]
     Nodexlist<-res[[3]]
     Nodeylist<-res[[4]]
     CoorespondIDs<-res[[6]]
     Eadf<-res[[5]]
     VElist<-res[[7]]
     edgelength<-res[[8]]
  }
  else
  {
    res<-footpoint.nodes(ntdata=ntdata,pointsxy=pointsxy, longlat=longlat,ea.prop=ea.prop, ELComputed=T)
    Nodeslist<-res[[1]]
    CoorespondIDs<-res[[6]]
    Edgeslist<-res[[4]]
    Nodexlist<-res[[2]]
    Nodeylist<-res[[3]]
    Eadf<-res[[5]]
    edgelength<-res[[7]]
    res<-virtualedge.nn(Nodeslist, Edgeslist, Nodexlist, Nodeylist, pointsxy, edgelength, Eadf,CoorespondIDs,longlat=longlat,ea.prop)
    Nodeslist<-res[[1]]
    Edgeslist<-res[[2]]
    Nodexlist<-res[[3]]
    Nodeylist<-res[[4]]
    CoorespondIDs<-res[[6]]
    Eadf<-res[[5]]
    VElist<-res[[7]]
    edgelength<-res[[8]] 
  }
  res<-list(Nodeslist, Edgeslist, CoorespondIDs, Nodexlist, Nodeylist, Eadf, VElist, edgelength)
  #Nodeslist 
}

##################
Nearest.nodes<-function(pointsxy, Nodeslist,longlat=F)
{
  Nn<-length(Nodeslist[,1])
  Np<-length(pointsxy[,1])
  CoorespondIDs<-c()
  for (i in 1:Np)
  {
     pxy<-c(as.double(pointsxy[i,1]),as.double(pointsxy[i,2]))
     tmpdist<-Inf
     id<-0
     for (j in 1:Nn)
     {
       nxy<-c(as.double(Nodeslist[j,2][[1]][1]),as.double(Nodeslist[j,2][[1]][2]))
       distnp<-as.double(0)
       distnp<-.C("edgelength", c(pxy[1],nxy[1]),c(pxy[2],nxy[2]), as.integer(2), distnp, as.integer(longlat))[[4]]
       if (distnp<tmpdist)
       {
         id<-Nodeslist[j,1]
         tmpdist<-distnp
       }       
     }
     CoorespondIDs<-c(CoorespondIDs, id)
  }
  CoorespondIDs
 }

##########################
footpoint.nodes<-function(ntdata,pointsxy, longlat=F, ea.prop, ELComputed=FALSE)
{
  if (longlat)
    warning("If the spatial coordinate is not projected to the Euclidean system, the acquired footpoints won't be accurate")
  Coords<-coordinates(ntdata)
  res<-readshpnw(data=ntdata, longlat<-longlat)
  nodelist<-res[[2]]
  nodexlist<-res[[6]]
  nodeylist<-res[[7]]
  edgelist<-res[[3]]
  Eadf<-res[[5]]
  #edgenum<-max(edgelist[,1])
  Np<-length(pointsxy[,1])
  CoorespondIDs<-c()
  #Changedegdes<-c()
  #newEdgesList<-c()
  #ftpts<-c()
  #numEdges<-length(Coords)
  edgelength<-c()
  #newEL<-c()
  if (is.null(ea.prop))
    stop("If a detailed graph is to be built, the properties of its attributes has to be specified")
  else
    {
     if (length(ea.prop)!=dim(Eadf)[2]-1)
     stop("All the properties of attributs should be specified except EdgeID")
    }
  for (i in 1:Np)
  {
    X<-as.double(pointsxy[i,1])
    Y<-as.double(pointsxy[i,2])
    distp2e<-as.double(1000000000000000000)
    d<-as.double(0)
    fx<-as.double(0)
    fy<-as.double(0)
    tag<-0
    ###Ver 1.1
    numEdges<-length(Coords)
    edgenum<-max(edgelist[,1])###For generating the ID of the new edge
    ##sidx: If the nearest point is located on a road segment(a polyline), sidx is the sequence NO. of the line segment in this polyline
    sidx<-as.integer(0)
    for (j in 1:numEdges)
    {
      M<-dim(Coords[[j]][[1]])[1]
      NodesXL<-as.double(Coords[[j]][[1]][,1])
      NodesYL<-as.double(Coords[[j]][[1]][,2])
      res<-.C("dist_p2e", NodesXL, NodesYL, as.integer(M),sidx, X, Y, d, fx, fy)
      d<-res[[7]]
      if (d<distp2e)
      {
        fx<-res[[8]]
        fy<-res[[9]]
        sidx<-res[[4]]
        distp2e<-d
        tag<-j  ##Save the index of the selected edge        
      }
    }
    ##To check if the found nearest point existed as a node
    res<-Update.nodelist(nodexlist,nodeylist,nodelist, fx, fy)
    if (res[[2]])
    {
       node<-list(res[[1]],c(fx,fy))
       nodelist<-rbind(nodelist,node)
       nodexlist<-c(nodexlist,fx)
       nodeylist<-c(nodeylist,fy) 
       CoorespondIDs<-c(CoorespondIDs, res[[1]])
       #Changedegdes<-c(Changedegdes, edgelist[tag,1])
       #res<-Add.nodes(Coords=Coords, fx=fx, fy=fy, edge.tag=tag,new.nodeID=res[[1]], edge=edgelist[tag,], sidx=sidx, edgenum=edgenum, newEdgesList=newEdgesList, longlat=longlat, newEL)
       ##ver1.1
       res<-Add.nodes(Coords=Coords, fx=fx, fy=fy, edge.tag=tag,new.nodeID=res[[1]], edge=edgelist[tag,], sidx=sidx, edgenum=edgenum,longlat=longlat,edgelist=edgelist,Eadf=Eadf,ea.prop=ea.prop)       
       Coords<-res[[1]]
       edgelist<-res[[2]]
       Eadf<-res[[3]]
    }
    else
    {
       CoorespondIDs<-c(CoorespondIDs, res[[1]]) 
    }     
  }
  #res<-extend.eadf(newEdgesList,Eadf, ea.prop)
  #newEadf<-res[[2]]
  #newE<-res[[1]]
  
  #Eadf<-rbind(Eadf, newEadf)
  #edgelist<-rbind(edgelist, newE)
  #changedidxs<-NULL
  #if (!is.null(Changedegdes))
  #{
   # changedidxs<-which(edgelist[,1]==Changedegdes)
   # edgelist<-edgelist[, -changedidxs]
   # Eadf<-Eadf[, -changedidxs]   
  #}
  
  if (ELComputed)
  {   
     numEdges<-length(Coords)
     for(i in 1:numEdges)
     {
        M<-dim(Coords[[i]][[1]])[1]
        El<-as.double(0)
        El<- .C("edgelength", as.double(Coords[[i]][[1]][,1]),as.double(Coords[[i]][[1]][,2]), as.integer(M), El, as.integer(longlat))[[4]]
        edgelength<-c(edgelength,El)
     }
    # if(is.null(changedidxs))
     #{
      #  edgelength<-edgelength[-changedidxs]
       # edgelength<-c(edgelength, newEL)
    # } 
    #edgelength<-c(edgelength, newEL)
  } 
  rownames(nodelist)<-NULL
  rownames(edgelist)<-NULL
  res<-list(nodelist,nodexlist,nodeylist,edgelist,Eadf,CoorespondIDs, edgelength)
  res  
}

Add.nodes<-function(Coords, fx, fy, edge.tag,new.nodeID, edge, sidx, edgenum, longlat=F,edgelist,Eadf,ea.prop)
{
  #if (is.null(newEdgesList))
#      n<-0
#  else
#      n<-length(newEdgesList)
  M<-dim(Coords[[edge.tag]][[1]])[1]
  El<-as.double(0)
  El1<-as.double(0)
  E1coord<-rbind(Coords[[edge.tag]][[1]][1:sidx,],c(fx, fy))
  E2coord<-rbind(c(fx, fy), Coords[[edge.tag]][[1]][(sidx+1):M,])
  El<- .C("edgelength", as.double(Coords[[edge.tag]][[1]][,1]),as.double(Coords[[edge.tag]][[1]][,2]), as.integer(M), El, as.integer(longlat))[[4]]
  El1<- .C("edgelength", as.double(E1coord[,1]),as.double(E1coord[,2]), as.integer(sidx+1), El1, as.integer(longlat))[[4]]
  ## newEL: lengths of the new edge
  
  #newEL<-c(newEL, El1, El-El1)
  new.edge<-c(edgenum+1, edge[1], edge[2], new.nodeID, El1/El)
  res<-extend.eadf(new.edge,Eadf, ea.prop)
  Eadf<-rbind(Eadf, res[[2]])
  ##ver1.1
  Coords[[edge.tag]]<-NULL
  edgelist<-edgelist[-edge.tag,]
  edgelist<-rbind(edgelist,c(edgenum+1, edge[2], new.nodeID))
  edgelist<-rbind(edgelist,c(edgenum+2, new.nodeID,edge[3]))
  Coords[[length(Coords)+1]]<-list(E1coord)
  Coords[[length(Coords)+1]]<-list(E2coord)
  #print(c(length(Coords),dim(edgelist)[1]))
  #newEdgesList<-rbind(newEdgesList, new.edge)
  new.edge<-c(edgenum+2, edge[1], new.nodeID,edge[3], (El-El1)/El)
  res<-extend.eadf(new.edge,Eadf, ea.prop)
  Eadf<-rbind(Eadf, res[[2]])
  Eadf<-Eadf[-edge.tag,]
  #newEdgesList<-rbind(newEdgesList, new.edge)
  res<-list(Coords, edgelist, Eadf)
  res
}

virtualedge.nn<-function(nodelist, edgeslist, nodexlist, nodeylist, pointsxy, edgelength, Eadf, CoorespondIDs, longlat=longlat, ea.prop)
{
   Np<-length(pointsxy[,1])
   eid<-max(edgeslist[,1])
   NewNIDs<-c()
   newEdgesList<-c()
   for (i in 1:Np)
   {
     res<-Update.nodelist(nodexlist, nodeylist,nodelist, pointsxy[i,1], pointsxy[i,2])
     if (res[[2]])
      {
       node<-list(res[[1]],c(pointsxy[i,1], pointsxy[i,2]))
       nodelist<-rbind(nodelist,node)
       nodexlist<-c(nodexlist,pointsxy[i,1])
       nodeylist<-c(nodeylist,pointsxy[i,2]) 
       NewNIDs<-c(NewNIDs, res[[1]])
      }
   edge<-c(eid+i, CoorespondIDs[i], res[[1]])
   newEdgesList<-rbind(newEdgesList, edge)
   }
   newEadf<-data.frame(matrix(newEdgesList[,1], ncol=1))  
   Nea<-length(ea.prop)
   Nne<-length(newEdgesList[,1])
   Eadf.names<-names(Eadf)
   Eadf1<-Eadf[,-1]
   newEL<-c()
   for (i in 1: Nea)
   {
      newdfi<-c()
      if(ea.prop[i]==0)
      {
        for (j in 1:Nne)
        {
           newdfi<-c(newdfi, NA)
        }
      }
      else
      {
        for (j in 1:Nne)
        {
          nel<-as.double(0)
          nidx<-which(nodelist[,1]==newEdgesList[j,][[2]])
          oidx<-which(nodelist[,1]==newEdgesList[j,][[3]])
          nel<-.C("edgelength", as.double(nodexlist[c(nidx, oidx)]),as.double(nodeylist[c(nidx, oidx)]), as.integer(2), nel, as.integer(longlat))[[4]]
          newEL<-c(newEL, nel)
          Oeidxs<-c(which(edgeslist[,2]==newEdgesList[j,][[3]]), which(edgeslist[,3]==newEdgesList[j,][[3]]))
          oels<-edgelength[Oeidxs]
          w<-nel/oels
          ai<-Eadf1[Oeidxs,i]
          newdfi<-c(newdfi, sum(w*ai)/length(ai))
        }
      }
      newEadf<-data.frame(newEadf, newdfi)
   }
   rownames(newEadf)<-NULL
   rownames(newEdgesList)<-NULL
   rownames(newEL)<-NULL
   names(newEadf)<-Eadf.names
   Eadf<-rbind(Eadf, newEadf)
   edgeslist<-rbind(edgeslist, newEdgesList)
   edgelength<-c(edgelength, newEL)
   CoorespondIDs<-NewNIDs
   res<-list(nodelist, edgeslist, nodexlist, nodeylist, Eadf, CoorespondIDs, newEdgesList, edgelength)
   res 
}
##########################
point.in.bbox<-function(pointxy,bbox)
{
  if (pointxy[1]<=bbox[1,2]&&pointxy[1]>=bbox[1,1])
     if (pointxy[2]<=bbox[2,2]&&pointxy[2]>=bbox[2,1])
         res<-TRUE
     else
         {
         res<-FALSE
         break
         }
  else
         res<-FALSE
  res  
}

ptsinnt.view<-function(ntdata, nodelist, pointsxy, CoorespondIDs, VElist=NULL)
{ 
  ntbbox<-bbox(ntdata)
  ptbbox<-bbox(ntdata)
  margin<-c(ntbbox[1,2]-ntbbox[1,1],ntbbox[2,2]-ntbbox[2,1])/20
  x.lim<-c(min(ntbbox[1,], ptbbox[1,])-margin[1], max(ntbbox[1,], ptbbox[1,])+margin[1])
  y.lim<-c(min(ntbbox[2,], ptbbox[2,])-margin[2], max(ntbbox[2,], ptbbox[2,])+margin[2])
  nodesxy<-Nodes.coordinates(nodelist)
  plot(ntdata, xlim=x.lim, ylim=y.lim)
  points(nodesxy[,1],nodesxy[,2],pch=1,col="grey",cex=0.2)
  points(pointsxy[,1],pointsxy[,2],pch=18,col="green")
  np<-length(CoorespondIDs)
  nn<-length(nodelist[,1])
  Nids<-c()
  for(i in 1:nn)
  {
    Nids<-c(Nids, nodelist[i,1][[1]])
  }
  for (i in 1:np)
  {
    idx<-which(Nids==CoorespondIDs[i])
    points(nodesxy[idx,1],nodesxy[idx,2],pch=1, col="red")
    lines(c(nodesxy[idx,1],pointsxy[i,1]),c(nodesxy[idx,2],pointsxy[i,2]), lty=2,col="blue")
  }
  if (!is.null(VElist))
  {
    Nve<-length(VElist[,1])
    for (i in 1:Nve)
    {
      idx1<-which(Nids==VElist[i,][[2]])
      idx2<-which(Nids==VElist[i,][[3]])
      lines(c(nodesxy[idx1,1],nodesxy[idx2,1]),c(nodesxy[idx1,2],nodesxy[idx2,2]),col="green")
    }
  }
  cols<-c("green","red")
  pchs<-c(18,1)
  legend(x<-ntbbox[1,1], y =ntbbox[2,2], c("Data point", "Nearest point on the network"), col = cols, pch=pchs)
}