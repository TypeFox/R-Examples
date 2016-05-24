readshpnw<-function(data=list(), ELComputed=FALSE, longlat=FALSE, Detailed=FALSE, ea.prop=NULL)
#The function returns a data list composed of "nodelist", "edgelist" and "data frame of edge attributes" 
{
  if(!is(data,"SpatialLinesDataFrame"))
     stop("Input data is not a proper spatial network data frame, here only SpatialLinesDataFrame is accepted.")
  Coords<-coordinates(data)
  numEdges<-length(Coords)
  nodelist<-c()
  edgelist<-c()
  Eadf<-data.frame(data)
  Eadf.names<-names(Eadf)
  edgeID.vec<-1:numEdges
  Eadf<-data.frame(edgeID.vec,data.frame(data))
  names(Eadf)<-c("EdgeID", Eadf.names)
  id.idx<-1
  if (Detailed)
  {
    if (is.null(ea.prop))
    stop("If a detailed graph is to be built, the properties of its attributes has to be specified")
    else
      {
         if (length(ea.prop)!=dim(Eadf)[2]-1)
         stop("All the properties of attributs should be specified except EdgeID")
      } 
  }
    #####################################################
  
  nodexlist<-vector(mode="double",length=0)
  nodeylist<-vector(mode="double",length=0)
  edgelength<-vector(mode="double",length=0)
  fromid<-0
  toid<-0
##################################  
  if (Detailed)
  {
     for (i in 1:numEdges)
    {
      M<-dim(Coords[[i]][[1]])[1]
      SEl<-as.double(10)
      SEl<- .C("edgelength", as.double(Coords[[i]][[1]][,1]),as.double(Coords[[i]][[1]][,2]), as.integer(M), SEl, as.integer(longlat))[[4]]
      edgeid<-edgeID.vec[i]
      nx<-as.double(Coords[[i]][[1]][1,1])
      ny<-as.double(Coords[[i]][[1]][1,2])
      res1<-Update.nodelist(nodexlist,nodeylist, nodelist, nx, ny)
      fromid<-res1[[1]]
      if(res1[[2]])
      {
        node<-list(res1[[1]],c(nx,ny))
        nodelist<-rbind(nodelist,node)
        nodexlist<-c(nodexlist,nx)
        nodeylist<-c(nodeylist,ny)
      }
      for (j in 2:M)
       {            
          nx<-as.double(Coords[[i]][[1]][j,1])
          ny<-as.double(Coords[[i]][[1]][j,2])
          res1<-Update.nodelist(nodexlist,nodeylist, nodelist, nx, ny)
          toid<-res1[[1]]
          if(res1[[2]])
           {
            node<-list(res1[[1]],c(nx,ny))
            nodelist<-rbind(nodelist,node)
            nodexlist<-c(nodexlist,nx)
            nodeylist<-c(nodeylist,ny)
           }
         ####Compute the edgelength
           El<-as.double(0)
           El<- .C("edgelength", as.double(Coords[[i]][[1]][c(j-1, j),1]),as.double(Coords[[i]][[1]][c(j-1, j),2]), as.integer(2), El, as.integer(longlat))[[4]]
           ROL<-El/SEl
           edgelist<-Update.edgelist(edgelist, edgeid, fromid, toid, ROL,Detailed=Detailed)
           if(ELComputed)
           {
             edgelength<-c(edgelength,El)
           }
           fromid<-toid
       }  
    }    
  }
  else
  {
    for (i in 1:numEdges)
    {
      M<-dim(Coords[[i]][[1]])[1]
      nx<-as.double(Coords[[i]][[1]][1,1])
      ny<-as.double(Coords[[i]][[1]][1,2])
      res1<-Update.nodelist(nodexlist,nodeylist, nodelist, nx, ny)
      fromid<-res1[[1]]
      if(res1[[2]])
      {
        node<-list(res1[[1]],c(nx,ny))
        nodelist<-rbind(nodelist,node)
        nodexlist<-c(nodexlist,nx)
        nodeylist<-c(nodeylist,ny)
      }   
      nx<-as.double(Coords[[i]][[1]][M,1])
      ny<-as.double(Coords[[i]][[1]][M,2])
      res1<-Update.nodelist(nodexlist,nodeylist, nodelist, nx, ny)
      toid<-res1[[1]]
      if(res1[[2]])
      {
        node<-list(res1[[1]],c(nx,ny))
        nodelist<-rbind(nodelist,node)
        nodexlist<-c(nodexlist,nx)
        nodeylist<-c(nodeylist,ny)
      }
      edgeid<-edgeID.vec[i]
      edgelist<-Update.edgelist(edgelist, edgeid, fromid, toid, Detailed=Detailed)
       ####Compute the edgelength
           El<-as.double(0)
           if(ELComputed)
           {
             El<- .C("edgelength", as.double(Coords[[i]][[1]][,1]),as.double(Coords[[i]][[1]][,2]), as.integer(M), El, as.integer(longlat))[[4]]
             edgelength<-c(edgelength,El)
           } 
    }    
  }
  if (Detailed)
  {
     res<-extend.eadf(edgelist,Eadf, ea.prop)
     edgelist<-res[[1]]
     Eadf<-res[[2]]
  }
  rownames(edgelist)<-NULL
  rownames(nodelist)<-NULL
  res<-list(Detailed, nodelist, edgelist, edgelength, Eadf, nodexlist, nodeylist)  
}

Update.nodelist<-function(nodexlist,nodeylist, nodelist, nx, ny)
{
   Nid<-length(nodexlist)
   if(Nid==0)
   {
     Nid<-Nid+1
     id<-as.integer(Nid)
     isUpdate<-TRUE
   }
   else
   {
     tag<-as.integer(-1)
     tag <- .C("nodeExisted", nodexlist,nodeylist, as.integer(Nid), nx, ny, tag)[[6]]
     if (tag!=-1)
        {
          id<-tag
          isUpdate<-FALSE
        }
        else
        {
          Nid<-Nid+1
          id<-as.integer(Nid)
          isUpdate<-TRUE
        }
   }
   
   res<-list(id,isUpdate)
   res
}

Update.edgelist<-function(Edgelist, edgeid, fromid, toid, ROL,Detailed)
{
   if (Detailed)
   {
     eid<-as.integer((length(Edgelist)%/%4)+1)
     edge<-c(eid, edgeid, fromid, toid, ROL)  
   }
   else
      edge<-c(edgeid, fromid, toid)
   Edgelist<-rbind(Edgelist, edge)
   Edgelist
}

Nodes.coordinates<-function(nodelist=list())
{
  Nn<-length(nodelist[,1])
  Nodesx<-vector(mode="double",length=0)
  Nodesy<-vector(mode="double",length=0)
  for (i in 1:Nn)
  {
     Nodesx<-c(Nodesx, as.double(nodelist[i,2][[1]][1]))
     Nodesy<-c(Nodesy, as.double(nodelist[i,2][[1]][2]))
  }
  Nodesxy<-cbind(Nodesx, Nodesy)
  colnames(Nodesxy)<-c("X","Y")
  Nodesxy 
}

extend.eadf<-function(newedges,Eadf, ea.prop)
{
  ##ver1.1
  if (is.null(dim(newedges)))
      newedges<-matrix(newedges, nrow=1,byrow=T)
  ###
  id.idx<-1
  ne<-length(newedges[,1])
  edgeIDs<-Eadf[,id.idx]
  Eadf.names<-names(Eadf)
  Eadf<-Eadf[,-id.idx]
  ndf<-length(ea.prop)
  newdf<-newedges[,1]
  for (i in 1:ndf)
  {
    dfi<-Eadf[,i]
    newdfi<-c()
    for (j in 1:ne)
    {
      edgeid<-newedges[j,2]
      edgeidx<-which(edgeIDs==edgeid)
      if (ea.prop[i]==0)
      {
        newdfi<-c(newdfi, dfi[edgeidx])
      }
      else
      {
        w<-newedges[j, 5]
        newdfi<-c(newdfi, as.numeric(dfi[edgeidx])*w)
      } 
    }                                     
    newdf<-cbind(newdf, newdfi)
  }
  newdf<-data.frame(newdf)
  names(newdf)<-Eadf.names
  newedges<-newedges[,-c(2,5)]
  res<-list(newedges,newdf)
  res
}