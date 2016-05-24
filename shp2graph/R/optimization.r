#Self-loop extraction
SL.extraction<-function(nodelist, edgelist, eadf=NULL,Directed=F, DegreeL=NULL, InDegreeL=NULL, OutDegreeL=NULL,Nexception=NULL,Eexception=NULL)
{
  Nn<-length(nodelist[,1])
  Ne<-length(edgelist[,1])
  newEdgelist<-c()
  newEadf<-c()
  if (Directed)
  {
    res<-Degree.list(nodelist, edgelist, Directed=Directed)
    InDegreeL<-res[[1]]
    OutDegreeL<-res[[2]]
  }
  else
  {
     DegreeL<-Degree.list(nodelist, edgelist, Directed=Directed)[[1]]
  }   
  for (i in 1:Ne)
  {
    if(edgelist[i,2]==edgelist[i,3]&&!is.exception(edgelist[i,1],Eexception))
    {
      if (Directed)
      { 
         res<-Minus.DegreeL(nodelist, edgelist[i,], Directed=Directed, InDegreeL=InDegreeL, OutDegreeL=OutDegreeL)
         InDegreeL<-res[[1]]
         OutDegreeL<-res[[2]]
      }
      else
      {
        res<-Minus.DegreeL(nodelist, edgelist[i,], Directed=Directed, DegreeL=DegreeL)
        DegreeL<-res[[1]]
      }
    }
    else
    {
        newEdgelist<-rbind(newEdgelist, edgelist[i,])
        newEadf<-rbind(newEadf, eadf[i,])
    }
  }
  rownames(newEdgelist)<-NULL
  rownames(newEadf)<-NULL
  newEadf<-data.frame(newEadf)
  names(newEadf)<-names(eadf)
  if (Directed)
   {
     De0<-which((InDegreeL+OutDegreeL)==0)
     if (length(De0)>0)
     {
       res<-is.exception(nodelist[De0,1],Nexception)
       De0<-De0[!res]
       if(length(De0)>0)
       {
        InDegreeL<-InDegreeL[-De0]
        OutDegreeL<-OutDegreeL[-De0]
        nodelist<-nodelist[-De0,]
       }
     }
     res<-list(nodelist, newEdgelist, newEadf,InDegreeL, OutDegreeL)
   }
   else
   {
     De0<-which(DegreeL==0)
     if (length(De0)>0)
     {
       res<-is.exception(nodelist[De0,1],Nexception)
       De0<-De0[!res]
       if(length(De0)>0)
       {
         DegreeL<-DegreeL[-De0]
         nodelist<-nodelist[-De0,]
       }
     }
     res<-list(nodelist, newEdgelist,newEadf, DegreeL)
   }
   res   
}

#Multiple-edges simplification
ME.simplification<-function(nodelist, edgelist, eadf=NULL, ea.prop=NULL, Directed=F, DegreeL=NULL, InDegreeL=NULL, OutDegreeL=NULL,Nexception=NULL,Eexception=NULL)
{
  newEadf<-c()
  if (!is.null(eadf))
    {
     if (length(ea.prop)!=dim(eadf)[2])
     stop("All the properties of attributs should be specified")
     newEadf<-c(eadf[1,])
    }  
  Nn<-length(nodelist[,1])
  Ne<-length(edgelist[,1])
  newEdgelist<-c()
  newEdgelist<-rbind(newEdgelist,edgelist[1,]) 
  if (Directed)
  {
    res<-Degree.list(nodelist, edgelist, Directed=Directed)
    InDegreeL<-res[[1]]
    OutDegreeL<-res[[2]]
  }
  else
  {
     DegreeL<-Degree.list(nodelist, edgelist, Directed=Directed)[[1]]
  }    
  for (i in 2:Ne)
  {
    Nne<-length(newEdgelist[,1])
    tag<-F
    edge<-edgelist[i,]
    if (Directed)
    {
      for (j in 1:Nne)
      {
        if(edge[2]==newEdgelist[j,2]&&edge[3]==newEdgelist[j,3]&&!is.exception(edge[1],Eexception))
        {
           tag<-T
           if (!is.null(eadf))
           {
             na<-length(ea.prop)
             eai<-c()
             for (k in 1:na)
             {
                eai<-c(eai,Redef.functions(c(eadf[i,k],newEadf[j,k]),typ=ea.prop[k]))
             }
            newEadf[j,]<-eai 
           }
           res<-Minus.DegreeL(nodelist, edge, Directed=Directed, InDegreeL=InDegreeL, OutDegreeL=OutDegreeL)
           InDegreeL<-res[[1]]
           OutDegreeL<-res[[2]]
           break
        } 
      }
      if(!tag)
      {
        newEdgelist<-rbind(newEdgelist,edge)
        if(!is.null(eadf))
        newEadf<-rbind(newEadf, eadf[i,])
      }      
    }
    else
    {
      for (j in 1:Nne)
      {
        if((edge[2]==newEdgelist[j,2]&&edge[3]==newEdgelist[j,3])||(edge[2]==newEdgelist[j,3]&&edge[3]==newEdgelist[j,2])&&!is.exception(edge[1],Eexception))
        {
           tag<-T
           if (!is.null(eadf))
           {
             na<-length(ea.prop)
             eai<-c()
             for (k in 1:na)
             {
                eai<-c(eai,Redef.functions(c(eadf[i,k],newEadf[j,k]),typ=ea.prop[k]))
             }
            newEadf[j,]<-eai 
           }
           res<-Minus.DegreeL(nodelist, edge, Directed=Directed, DegreeL=DegreeL)
           DegreeL<-res[[1]]
           break
        } 
      }
      if(!tag)
      {
        newEdgelist<-rbind(newEdgelist,edge)
        if(!is.null(eadf))
        newEadf<-rbind(newEadf, eadf[i,])
      } 
    }
  }
  rownames(newEdgelist)<-NULL
  rownames(newEadf)<-NULL
  newEadf<-data.frame(newEadf)
  names(newEadf)<-names(eadf)
  if (Directed)
   {
     res<-list(nodelist, newEdgelist, newEadf, InDegreeL, OutDegreeL)
   }
   else
   {
     res<-list(nodelist, newEdgelist, newEadf, DegreeL)
   }
   res   
}

#Pseudo-node amalgamation
PN.amalgamation<-function(nodelist, edgelist, eadf=NULL, ea.prop=NULL, Directed=F, DegreeL=NULL, InDegreeL=NULL, OutDegreeL=NULL, Nexception=NULL,Eexception=NULL)
{
  newEadf<-c()
  if (!is.null(eadf))
    {
     if (length(ea.prop)!=dim(eadf)[2])
     stop("All the properties of attributs should be specified")
     newEadf<-c(eadf[1,])
    }  
  Nn<-length(nodelist[,1])
  Ne<-length(edgelist[,1])
  if (Directed)
  {
    res<-Degree.list(nodelist, edgelist, Directed=Directed)
    InDegreeL<-res[[1]]
    OutDegreeL<-res[[2]]
  }
  else
  {
     DegreeL<-Degree.list(nodelist, edgelist, Directed=Directed)[[1]]
  }
  newNodelist<-c()
  NIDs<-as.integer(nodelist[,1])
  Nn<-as.integer(length(nodelist[,1]))        
  for (i in 1:Nn)
  {
     nid<-nodelist[i,1]
    if (Directed)
      {
         if(InDegreeL[i]==1&&OutDegreeL[i]==1&&!is.exception(nid,Nexception))
         {
            eidx1<-which(edgelist[,2]==nid)
            eidx2<-which(edgelist[,3]==nid)
            edge1<-edgelist[eidx1,]
            edge2<-edgelist[eidx2,]
            if ((!is.exception(edge1[1],Eexception))&&(!is.exception(edge2[1],Eexception)))
            {
              res<-Minus.DegreeL(nodelist, edge1, Directed=Directed, InDegreeL=InDegreeL, OutDegreeL=OutDegreeL)
              InDegreeL<-res[[1]]
              OutDegreeL<-res[[2]]
              res<-Minus.DegreeL(nodelist, edge2, Directed=Directed, InDegreeL=InDegreeL, OutDegreeL=OutDegreeL)
              InDegreeL<-res[[1]]
              OutDegreeL<-res[[2]]
              edge<-c(edge1[1],edge1[2],edge2[3])
              edgelist[eidx1,]<-edge
              edgelist<-edgelist[-eidx2, ]
              OutDegreeL<-.C("addDegree", NIDs, Nn, as.integer(edge[2]), OutDegreeL)[[4]]
              InDegreeL<-.C("addDegree", NIDs, Nn, as.integer(edge[3]), InDegreeL)[[4]]
              if (!is.null(eadf))
              {
                na<-length(ea.prop)
                eai<-c()
                for (k in 1:na)
                  {
                   eai<-c(eai,Redef.functions(c(eadf[eidx1,k],eadf[eidx2,k]),typ=ea.prop[k]))
                  }
                eadf[eidx1,]<-eai
                eadf<-eadf[-eidx2,] 
              }
            }
            else newNodelist<-rbind(newNodelist, nodelist[i,])
                 
         }
         else
         {
           newNodelist<-rbind(newNodelist, nodelist[i,])
         }
      }
   else
      { 
        if(DegreeL[i]==2&&!is.exception(nodelist[i,1],Nexception))
        {
            eidxs<-which(c(edgelist[,2],edgelist[,3])==nid)
            Endnids<-c()
            Ne<-length(edgelist[,1])
            for (k in 1:2)
            {
              if (eidxs[k]<=Ne) 
              {
                Endnids<-c(Endnids, edgelist[eidxs[k],3])
              }
              else 
              {  
                eidxs[k]<-eidxs[k]-Ne
                Endnids<-c(Endnids, edgelist[eidxs[k],2])
              }
            }
            if ((!is.exception(edgelist[eidxs[1],1],Eexception))&&(!is.exception(edgelist[eidxs[2],1],Eexception)))
            {
              res<-Minus.DegreeL(nodelist, edgelist[eidxs[1],], Directed=Directed, DegreeL=DegreeL)
              DegreeL<-res[[1]]
              res<-Minus.DegreeL(nodelist, edgelist[eidxs[2],], Directed=Directed, DegreeL=DegreeL)
              DegreeL<-res[[1]]
              edge<-c(edgelist[eidxs[1],1],Endnids[1],Endnids[2])
              edgelist[eidxs[1],]<-edge
              edgelist<-edgelist[-eidxs[2], ]
              DegreeL<-.C("addDegree", NIDs, Nn, as.integer(edge[2]), DegreeL)[[4]]
              DegreeL<-.C("addDegree", NIDs, Nn, as.integer(edge[3]), DegreeL)[[4]]
              if (!is.null(eadf))
              {
                na<-length(ea.prop)
                eai<-c()
                if (na==1)
                {
                  eai<-Redef.functions(c(eadf[eidxs[1],1],eadf[eidxs[2],1]),typ=ea.prop)
                  eadf[eidxs[1],1]<-eai
                  eadf<-eadf[-eidxs[2]]
                }
                else
                {
                  for (k in 1:na)
                  {
                   eai<-c(eai,Redef.functions(c(eadf[eidxs[1],k],eadf[eidxs[2],k]),typ=ea.prop[k]))
                  }
                  eadf[eidxs[1],]<-eai
                  eadf<-eadf[-eidxs[2],]
                }   
              }
            }
            else newNodelist<-rbind(newNodelist, nodelist[i,])
        }
        else
        {
          newNodelist<-rbind(newNodelist, nodelist[i,])
        } 
      }
  }
  rownames(edgelist)<-NULL
  rownames(eadf)<-NULL
  rownames(newNodelist)<-NULL
  if (Directed)
   {
     dl<-InDegreeL+OutDegreeL
     idxs<-which(dl==0)
     if (length(idxs)!=0)
     {
       InDegreeL<-InDegreeL[-idxs]
       OutDegreeL<-OutDegreeL[-idxs]
     }
     res<-list(newNodelist, edgelist, eadf,InDegreeL, OutDegreeL)
   }
   else
   {
     idxs<-which(DegreeL==0)
     if (length(idxs)!=0)
     {
       InDegreeL<-InDegreeL[-idxs]
       OutDegreeL<-OutDegreeL[-idxs]
     }
     res<-list(newNodelist, edgelist,eadf,DegreeL)
   }
   res   
}

Degree.list<-function(nodelist, edgelist, Directed=F)
{
   Nn<-as.integer(length(nodelist[,1]))
   Ne<-as.integer(length(edgelist[,1]))
   NIDs<-as.integer(nodelist[,1])
   if (Directed)
   {
     InDegreeL<-vector(Nn, mode="integer")
     OutDegreeL<-vector(Nn, mode="integer")
   }
   else
   {
     DegreeL<-vector(Nn, mode="integer")
   }
   for (i in 1:Ne)
   {
     if (Directed)
     {
        OutDegreeL<-.C("addDegree", NIDs, Nn, as.integer(edgelist[i,2]), OutDegreeL)[[4]]
        InDegreeL<-.C("addDegree", NIDs, Nn, as.integer(edgelist[i,3]), InDegreeL)[[4]]
     }
     else
     {
       DegreeL<-.C("addDegree", NIDs, Nn, as.integer(edgelist[i,2]), DegreeL)[[4]]
       DegreeL<-.C("addDegree", NIDs, Nn, as.integer(edgelist[i,3]), DegreeL)[[4]]
     }
   }
   if (Directed)
   {
     res<-list(InDegreeL, OutDegreeL)
   }
   else
   {
     res<-list(DegreeL)
   }
   res   
}
# if edge is to be removed, update the degrees of both nodes, i.e. minus 1 respectively
Minus.DegreeL<-function(nodelist, edge, Directed=F, DegreeL=NULL, InDegreeL=NULL, OutDegreeL=NULL)
{
   NIDs<-as.integer(nodelist[,1])
   Nn<-as.integer(length(nodelist[,1]))
   if (Directed)
   {
     InDegreeL=.C("minusDegree", NIDs, Nn, as.integer(edge[3]), InDegreeL)[[4]]
     OutDegreeL<-.C("minusDegree", NIDs, Nn, as.integer(edge[2]), OutDegreeL)[[4]]
     res<-list(InDegreeL, OutDegreeL)     
   }
   else
   {
     DegreeL<-.C("minusDegree", NIDs, Nn, as.integer(edge[2]), DegreeL)[[4]]
     DegreeL<-.C("minusDegree", NIDs, Nn, as.integer(edge[3]), DegreeL)[[4]]
     res<-list(DegreeL)
   }
   res
}
#Sum(): To redefine an edge with the summary of edge weights;---typ=1
#Min(): To choose or redefine an edge with the minimum edge weight from the multiple-edges;---typ=2 
#Max(): To choose or redefine an edge with the maximum edge weight;---typ=3
#Mean(): To redefine an edge with the mean of edge weights, and weighted mean could be also considered.---typ=4, w is for weighted Mean(), default is (1/n,...,1/n)
Redef.functions<-function(v, typ=1)
{
  if (typ<1||typ>4||typ%%1!=0)
  stop("Please specify the type of function using 1,2,3 or 4")
  if (length(v)==0)
  stop ("Vector v can't be empty")
  if (typ==1)
  {
    res<-sum(v)
  }
  else if(typ==2)
  {
    res<-min(v)
  }
  else if(typ==3)
  {
    res<-max(v) 
  }
  else
  {
    res<-sum(v)/length(v)
  }
  res 
}

is.exception<-function(idlist,exception)
{
  if(length(exception)==0)
  {
     if (length(idlist)==0)
     res<-NULL
     else
     res<-rep(F, length(idlist))
  }
  else
  {
     if (length(idlist)==0)
     res<-NULL
     else
     {
       res<-c()
       for (i in 1:length(idlist))
       {
         idx<-which(exception==idlist)
         if(length(idx)==0)
         res<-c(res, F)
         else
         res<-c(res, T)         
       }
     }
  }
  res
}