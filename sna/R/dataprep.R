######################################################################
#
# dataprep.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/27/13
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains various routines for preparing/preprocessing data
# for use with the sna package.
#
# Contents:
#
#   add.isolates
#   as.edgelist.sna
#   as.sociomatrix.sna
#   diag.remove
#   ego.extract
#   event2dichot
#   gt
#   gvectorize
#   interval.graph
#   is.edgelist.sna
#   lower.tri.remove
#   make.stochastic
#   nties
#   sr2css
#   stackcount
#   symmetrize
#   upper.tri.remove
#
######################################################################


#add.isolates - Add isolates to one or more graphs
add.isolates<-function(dat,n,return.as.edgelist=FALSE){
  if(!return.as.edgelist){
    #Pre-process the raw input
    dat<-as.sociomatrix.sna(dat)
    if(is.list(dat))
      return(lapply(dat,add.isolates,n=n,return.as.edgelist=return.as.edgelist))
    #End pre-processing
    if(length(dim(dat))>2){
      d<-array(dim=c(dim(dat)[1],dim(dat)[2]+n,dim(dat)[3]+n))
      d[,,]<-0
      for(i in 1:dim(dat)[1])
         d[i,1:dim(dat)[2],1:dim(dat)[2]]<-dat[i,,]
    }else{
      d<-matrix(nrow=dim(dat)[1]+n,ncol=dim(dat)[2]+n)
      d[,]<-0
      d[1:dim(dat)[2],1:dim(dat)[2]]<-dat
    }   
    d
  }else{
    #Pre-process the raw input
    dat<-as.edgelist.sna(dat)
    if(is.list(dat))
      return(lapply(dat,add.isolates,n=n,return.as.edgelist=return.as.edgelist))
    #End pre-processing
    attr(dat,"n")<-attr(dat,"n")+n
    dat
  }
}


#Force the input into edgelist form.  Network size, directedness, and vertex
#names are stored as attributes, since they cannot otherwise be included
as.edgelist.sna<-function(x, attrname=NULL, as.digraph=TRUE, suppress.diag=FALSE, force.bipartite=FALSE){
  #In case of lists, process independently
  if(is.list(x)&&(!(class(x)%in%c("network","matrix.csr","matrix.csc", "matrix.ssr","matrix.ssc", "matrix.hb","data.frame"))))
    return(lapply(x,as.edgelist.sna, attrname=attrname,  as.digraph=as.digraph, suppress.diag=suppress.diag, force.bipartite=force.bipartite))
  #Begin with network objects
  if(class(x)=="network"){
    require("network")  #Must have network library to process network objects
    out<-as.matrix.network.edgelist(x,attrname=attrname,as.sna.edgelist=TRUE)
    #This should be fine unless we have an old version of network (<1.7);
    #here, we perform triage for old style objects.
    if(!("as.sna.edgelist"%in%names(formals(as.matrix.network.edgelist)))){
      if(NCOL(out)==2)                        #If needed, add edge values
        out<-cbind(out,rep(1,NROW(out)))
      if(suppress.diag&&has.loops(x))
        out<-out[!(out[,1]==out[,2]),]
      if((!is.directed(x))&&as.digraph){
        if(has.loops(x)){
          temp<-out[,1]==out[,2]
          if(any(temp)){
            temp2<-out[temp,]
            out<-out[!temp,]
            out<-rbind(out,out[,c(2,1,3)])
            out<-rbind(out,temp2)
          }else
            out<-rbind(out,out[,c(2,1,3)])
        }else
          out<-rbind(out,out[,c(2,1,3)])
      }
      attr(out,"n")<-network.size(x)
      attr(out,"vnames")<-network.vertex.names(x)
    }
    if(is.bipartite(x)) #Unneeded for new objects, but does no harm
      attr(out,"bipartite")<-get.network.attribute(x,"bipartite")
    else if(force.bipartite)
      out<-as.edgelist.sna(out,attrname=attrname,as.digraph=as.digraph, suppress.diag=suppress.diag,force.bipartite=force.bipartite)
  } else
  #Not a network -- is this a sparse matrix (from SparseM)?
  if(class(x)%in%c("matrix.csr","matrix.csc","matrix.ssr","matrix.ssc", "matrix.hb")){
    require("SparseM")   #Need SparseM for this
    if(force.bipartite||(!is.null(attr(x,"bipartite")))|| (x@dimension[1]!=x@dimension[2])){
      nr<-x@dimension[1]
      nc<-x@dimension[2]
      val<-x@ra
      if((!suppress.diag)&&(class(x)%in%c("matrix.ssr","matrix.ssc"))){
        snd<-rep(1:nr,each=diff(x@ia))
        rec<-nr+x@ja
        out<-cbind(snd,rec,val)
        out<-rbind(out,out[,c(2,1,3)])
      }else{
        snd<-switch(class(x),
          matrix.csr=rep(1:nr,each=diff(x@ia)),
          matrix.csc=x@ja,
          matrix.ssr=c(rep(1:nr,each=diff(x@ia)),x@ja),
          matrix.ssc=c(x@ja,rep(1:nr,each=diff(x@ia)))
        )
        rec<-switch(class(x),
          matrix.csr=nr+x@ja,
          matrix.csc=rep(nr+(1:nc),each=diff(x@ia)),
          matrix.ssr=c(nr+x@ja,rep(1:n,each=diff(x@ia))),
          matrix.ssc=c(rep(nr+(1:nc),each=diff(x@ia)),x@ja)
        )
        out<-cbind(snd,rec,val)
        out<-rbind(out,out[,c(2,1,3)])
      }
      attr(out,"n")<-nr+nc
      attr(out,"vnames")<-NULL    #No dimnames for these objects
      attr(out,"bipartite")<-nr
    }else{
      n<-x@dimension[1]
      val<-x@ra
      if((!suppress.diag)&&(class(x)%in%c("matrix.ssr","matrix.ssc"))){
        snd<-rep(1:n,times=diff(x@ia))
        rec<-x@ja
        temp<-snd==rec
        out<-cbind(snd,rec,val)
        temp2<-out[temp,]
        out<-out[!temp,]
        out<-rbind(out,out[,c(2,1,3)])
        out<-rbind(out,temp2)
      }else{
        snd<-switch(class(x),
          matrix.csr=rep(1:n,times=diff(x@ia)),
          matrix.csc=x@ja,
          matrix.ssr=c(rep(1:n,times=diff(x@ia)),x@ja),
          matrix.ssc=c(x@ja,rep(1:n,times=diff(x@ia)))
        )
        rec<-switch(class(x),
          matrix.csr=x@ja,
          matrix.csc=rep(1:n,times=diff(x@ia)),
          matrix.ssr=c(x@ja,rep(1:n,times=diff(x@ia))),
          matrix.ssc=c(rep(1:n,times=diff(x@ia)),x@ja)
        )
        out<-cbind(snd,rec,val)
        if(suppress.diag)
          out<-out[!(out[,1]==out[,2]),]
      }
      attr(out,"n")<-n
      attr(out,"vnames")<-NULL    #No dimnames for these objects
    }
    if(force.bipartite&&(is.null(attr(out,"bipartite"))))
      out<-as.edgelist.sna(out,attrname=attrname,as.digraph=as.digraph, suppress.diag=suppress.diag,force.bipartite=force.bipartite)
  } else
  #Matrix or data frame case
  if(is.matrix(x)||is.data.frame(x)){ 
    if((NCOL(x)==3)&&(!is.null(attr(x,"n")))){  #Is this already an edgelist?
      out<-x
      if(force.bipartite&&(is.null(attr(out,"bipartite")))){ #Treat as bipartite
        out[,2]<-out[,2]+attr(x,"n")
        out<-rbind(out,out[,c(2,1,3)])
        attr(out,"n")<-attr(x,"n")*2
        attr(out,"bipartite")<-attr(x,"n")
        if(!is.null(attr(x,"vnames")))
          attr(out,"vnames")<-c(attr(x,"vnames"),attr(x,"vnames"))
        else
          attr(out,"vnames")<-NULL
      }
    }else if((NCOL(x)==2)&&(!is.null(attr(x,"n")))){  #Is this an edgelist w/out vals?
      out<-cbind(x,rep(1,NROW(x)))
      attr(out,"n")<-attr(x,"n")
      attr(out,"bipartite")<-attr(x,"bipartite")
      attr(out,"vnames")<-attr(x,"vnames")
      if(force.bipartite&&(is.null(attr(out,"bipartite")))){ #Treat as bipartite
        out[,2]<-out[,2]+attr(x,"n")
        out<-rbind(out,out[,c(2,1,3)])
        attr(out,"n")<-attr(x,"n")*2
        attr(out,"bipartite")<-attr(x,"n")
        if(!is.null(attr(x,"vnames")))
          attr(out,"vnames")<-c(attr(x,"vnames"),attr(x,"vnames"))
        else
          attr(out,"vnames")<-NULL
      }
    }else if(force.bipartite||(!is.null(attr(x,"bipartite")))|| (NROW(x)!=NCOL(x))){  #Assume this is a bipartite graph
      mask<-is.na(x)|(x!=0)
      if(sum(mask)>0){
        snd<-row(x)[mask]
        rec<-NROW(x)+col(x)[mask]
        val<-x[mask]
      }else{
        snd<-vector()
        rec<-vector()
        val<-vector()
      }
      out<-cbind(snd,rec,val)
      out<-rbind(out,out[,c(2,1,3)])
      attr(out,"n")<-NROW(x)+NCOL(x)
      attr(out,"vnames")<-c(rownames(x),colnames(x))
      attr(out,"bipartite")<-NROW(x)
    }else{                                 #Assume this is an adjmat
      mask<-is.na(x)|(x!=0)
      snd<-row(x)[mask]
      rec<-col(x)[mask]
      val<-x[mask]
      out<-cbind(snd,rec,val)
      attr(out,"n")<-NROW(x)
      attr(out,"vnames")<-rownames(x)
    }
  }else
  #Array case 
  if(is.array(x)){
      dx<-dim(x)
      ldx<-length(dx)
      if(ldx==2){                                   #Two-dimensional array
        if((dx[2]==3)&&(!is.null(attr(x,"n")))){  #Is this already an edgelist?
          out<-as.matrix(x)
          attr(out,"n")<-attr(x,"n")
          attr(out,"bipartite")<-attr(x,"bipartite")
          attr(out,"vnames")<-attr(x,"vnames")
        }
        if((NCOL(x)==2)&&(!is.null(attr(x,"n")))){  #Is this an edgelist w/out vals?
          out<-cbind(as.matrix(x),rep(1,NROW(x)))
          attr(out,"n")<-attr(x,"n")
          attr(out,"bipartite")<-attr(x,"bipartite")
          attr(out,"vnames")<-attr(x,"vnames")
        }else if(force.bipartite||(!is.null(attr(x,"bipartite")))|| (NROW(x)!=NCOL(x))){  #Assume this is a bipartite graph
          mask<-is.na(x)|(x!=0)
          if(sum(mask)>0){
            snd<-row(x)[mask]
            rec<-NROW(x)+col(x)[mask]
            val<-x[mask]
          }else{
            sna<-vector()
            rec<-vector()
            val<-vector()
          }
          out<-cbind(snd,rec,val)
          out<-rbind(out,out[,c(2,1,3)])
          attr(out,"n")<-NROW(x)+NCOL(x)
          attr(out,"vnames")<-c(dimnames(x)[[1]],dimnames(x)[[2]])
          attr(out,"bipartite")<-NROW(x)
        }else{                                 #Assume this is an adjmat
          mask<-is.na(x)|(x!=0)
          snd<-row(x)[mask]
          rec<-col(x)[mask]
          val<-x[mask]
          out<-cbind(snd,rec,val)
          attr(out,"n")<-NROW(x)
          attr(out,"vnames")<-dimnames(x)[[1]]
        }
        if(force.bipartite&&(is.null(attr(out,"bipartite")))){ #Treat as bipartite
          out[,2]<-out[,2]+attr(x,"n")
          out<-rbind(out,out[,c(2,1,3)])
          attr(out,"n")<-attr(x,"n")*2
          attr(out,"bipartite")<-attr(x,"n")
          if(!is.null(attr(x,"vnames")))
            attr(out,"vnames")<-c(attr(x,"vnames"),attr(x,"vnames"))
          else
            attr(out,"vnames")<-NULL
        }
      }else if(ldx==3){                           #Three-dimensional array
        out<-unlist(apply(x,1,function(z){list(as.edgelist.sna(z, attrname=attrname,as.digraph=as.digraph,suppress.diag=suppress.diag,force.bipartite=force.bipartite))}),recursive=FALSE)
      }else
        stop("Array input to as.edgelist.sna must either be a proper edgelist, an adjacency matrix, or an adjacency array.\n")
  }else{
    stop("as.edgelist.sna input must be an adjacency matrix/array, edgelist matrix, network, or sparse matrix, or list thereof.\n")
  }
  #Return the result
  out
}


#Force the input into sociomatrix form.  This function includes an sna
#wrapper to the network function as.sociomatrix, for global happiness.
as.sociomatrix.sna<-function(x, attrname=NULL, simplify=TRUE, force.bipartite=FALSE){
  #If passed a list, operate on each element
  if(is.list(x)&&(!(class(x)%in%c("network","matrix.csr","matrix.csc", "matrix.ssr","matrix.ssc", "matrix.hb","data.frame")))){
    g<-lapply(x,as.sociomatrix.sna,attrname=attrname,simplify=simplify, force.bipartite=force.bipartite)
    #Otherwise, start with network
  }else if(class(x)=="network"){
    require("network")  #Must have network library to process network objects
    g<-as.sociomatrix(x, attrname=attrname, simplify=simplify)
  #Not a network -- is this a sparse matrix (from SparseM)?
  }else if(class(x)%in%c("matrix.csr","matrix.csc","matrix.ssr","matrix.ssc", "matrix.hb")){
    require("SparseM")   #Need SparseM for this
    bip<-attr(x,"bipartite")
    g<-as.matrix(x)      #Coerce to matrix form, and pass on
    attr(g,"bipartite")<-bip
  }else{
    #Coerce to adjacency matrix form -- by now, no other classes involved
    n<-attr(x,"n")                        #Grab attributes before they get lost
    bip<-attr(x,"bipartite")
    vnam<-attr(x,"vnames")
    if(is.array(x)&&(length(dim(x))==2))  #Quick diversion for 2-d arrays
      x<-as.matrix(x)
    if(is.data.frame(x))                  #Coerce data frames to matrices
      x<-as.matrix(x)
    if(is.matrix(x)){
      if((NCOL(x)%in%c(2,3))&&(!is.null(n))){     #sna edgelist
        if(NCOL(x)==2)
          x<-cbind(x,rep(1,NROW(x)))
        g<-matrix(0,n,n)
        if(NROW(x)>0)
          g[x[,1:2,drop=FALSE]]<-x[,3]
        rownames(g)<-vnam
        colnames(g)<-vnam
      }else if(force.bipartite||(!is.null(bip))||(NROW(x)!=NCOL(x))){    #Bipartite adjmat
        nr<-NROW(x)
        nc<-NCOL(x)
        g<-matrix(0,nr+nc,nr+nc)
        g[1:nr,(nr+1):(nr+nc)]<-x
        g[(nr+1):(nr+nc),1:nr]<-t(x)
        rownames(g)<-vnam
        colnames(g)<-vnam
      }else{                                             #Regular adjmat
        g<-x
      }  
    }else if(is.array(x)){                   #If an array, test for type
      if(length(dim(x))!=3)
        stop("as.sociomatrix.sna input must be an adjacency matrix/array, network, data frame, sparse matrix, or list.")
      if(force.bipartite||(!is.null(attr(x,"bipartite")))|| (dim(x)[2]!=dim(x)[3])){      #Bipartite stack
        dx<-dim(x)
        nr<-dx[2]
        nc<-dx[3]
        g<-array(0,dim=c(dx[1],nr+nc,nr+nc))
        for(i in 1:dx[1]){
          g[i,1:nr,(nr+1):(nr+nc)]<-x[i,,]
          g[i,(nr+1):(nr+nc),1:nr]<-t(x[i,,])
        }
      }else{                                          #Adjacency stack
        g<-x
      }
    }else{
      stop("as.sociomatrix.sna input must be an adjacency matrix/array, network, or list.")
    }
  }
  #Convert into the appropriate return format
  if(is.list(g)){   #Collapse if needed
    if(length(g)==1){
      g<-g[[1]]
      if((!simplify)&&(length(dim(g))==3)){  #Coerce to a list of matrices?
        out<-list()
        for(i in 1:dim(g)[1])
          out[[i]]<-g[i,,]
      }else{
        out<-g
      }
    }else{
      #Coerce to array form?
      if(simplify){
        dims<-sapply(g,dim)
        if(is.list(dims)){      #Dims must not be of equal length
          mats<-sapply(dims,length)
          mats[mats==1]<-0
          mats[mats==2]<-1
          mats[mats==3]<-sapply(dims[mats==3],"[[",1)
          mats<-cumsum(mats)
          dims<-sapply(dims,"[",2)
        }else{                  #Dims are of equal length
          if(NROW(dims)==3)      #Determine number of matrices per entry
            mats<-cumsum(dims[1,])
          else
            mats<-1:NCOL(dims)
          dims<-dims[2,]         #Get ncols
        }
        if((!any(is.null(dims)))&&(length(unique(dims))==1)&&(all(mats>0))){
          out<-array(dim=c(mats[length(mats)],dims[1],dims[1]))
          for(i in 1:length(mats))
            out[(c(0,mats)[i]+1):(mats[i]),,]<-g[[i]]
        }else
          out<-g
      }else
        out<-g
    }
  }else{
    if((!simplify)&&(length(dim(g))==3)){  #Coerce to a list of matrices?
      out<-list()
      for(i in 1:dim(g)[1])
        out[[i]]<-g[i,,]
    }else
      out<-g
  }
  #Return the result
  out
}



#diag.remove - NA the diagonals of adjacency matrices in a graph stack
diag.remove<-function(dat,remove.val=NA){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(lapply(dat,diag.remove,remove.val=remove.val))
   #End pre-processing
   if(length(dim(dat))>2){
      d<-dat
      for(i in 1:dim(dat)[1])
         diag(d[i,,])<-remove.val
   }
   else{
      d<-dat
      diag(d)<-remove.val
   }   
   d
}


#ego.extract - Extract ego nets from an input graph, returning them as a
#list of graphs.
ego.extract<-function(dat,ego=NULL,neighborhood=c("combined","in","out")){
  #Pre-process the raw input
  d<-as.sociomatrix.sna(dat)
  if(is.list(d))
    return(lapply(d,ego.extract,ego=ego,neighborhood=neighborhood))
  else if(length(dim(dat))==3)
    return(apply(d,1,ego.extract,ego=ego,neighborhood=neighborhood))
  #End pre-processing
  #Set input arguments
  if(is.null(ego))
    ego<-1:NROW(d)
  neighborhood<-match.arg(neighborhood)  
  #Extract the selected ego nets
  enet<-list()
  for(i in 1:length(ego)){             #Walk the ego list
    sel<-switch(neighborhood,               #Grab the alters
      "in"=(1:NROW(d))[d[,ego[i]]>0],
      "out"=(1:NROW(d))[d[ego[i],]>0],
      "combined"=(1:NROW(d))[(d[ego[i],]>0)|(d[,ego[i]]>0)]
    )
    if(length(sel)>0)
      sel<-c(ego[i],sel[sel!=ego[i]])  #Force ego to be first
    else
      sel<-ego[i]
    enet[[i]]<-d[sel,sel,drop=FALSE]   #Perform the extraction
  }
  #Return the result
  if(!is.null(rownames(d)))          #Try to name the egos....
    names(enet)<-rownames(d)[ego]
  else if(!is.null(colnames(d)))
    names(enet)<-colnames(d)[ego]
  else
    names(enet)<-ego
  enet
}


#event2dichot - Convert an observed event matrix to a dichotomous matrix.  
#Methods are quantile, mean, rquantile, rmean, cquantile, cmean, absolute, rank,
#rrank, and crank.  Thresh specifies the cutoff, in terms of whatever method is 
#used (if applicable).
event2dichot<-function(m,method="quantile",thresh=0.5,leq=FALSE){
   #Pre-process the raw input
   m<-as.sociomatrix.sna(m)
   if(is.list(m))
     return(lapply(m,event2dichot,method=method,thresh=thresh,leq=leq))
   #End pre-processing
   if(method=="quantile"){
      q<-quantile(m,thresh,na.rm=TRUE, names=FALSE)
      out<-as.numeric(m>q)
   } else if(method=="rquantile"){
      q<-quantile(m[1,],thresh,na.rm=TRUE, names=FALSE)
      out<-as.numeric(m[1,]>q)
      for(i in 2:dim(m)[1]){      
         q<-quantile(m[i,],thresh,na.rm=TRUE, names=FALSE)
         out<-rbind(out,as.numeric(m[i,]>q))
      }
   } else if(method=="cquantile"){
      q<-quantile(m[,1],thresh,na.rm=TRUE, names=FALSE)
      out<-as.numeric(m[,1]>q)
      for(i in 2:dim(m)[2]){      
         q<-quantile(m[,i],thresh,na.rm=TRUE, names=FALSE)
         out<-cbind(out,as.numeric(m[,i]>q))
      }
   } else if(method=="mean"){
      q<-mean(m)
      out<-as.numeric(m>q)
   } else if(method=="rmean"){
      q<-mean(m[1,])
      out<-as.numeric(m[1,]>q)
      for(i in 2:dim(m)[1]){      
         q<-mean(m[i,])
         out<-rbind(out,as.numeric(m[i,]>q))
      }
   } else if(method=="cmean"){
      q<-mean(m[,1])
      out<-as.numeric(m[,1]>q)
      for(i in 2:dim(m)[2]){      
         q<-mean(m[,i])
         out<-rbind(out,as.numeric(m[,i]>q))
      }
   } else if(method=="absolute"){
      out<-as.numeric(m>thresh)
   } else if(method=="rank"){
      o<-order(m)
      out<-as.numeric((max(o)-o+1)<thresh)
   } else if(method=="rrank"){
      o<-order(m[1,])
      out<-as.numeric((max(o)-o+1)<thresh)
      for(i in 2:dim(m)[1]){      
         o<-order(m[i,])
         out<-rbind(out,as.numeric((max(o)-o+1)<thresh))
      }
   } else if(method=="crank"){
      o<-order(m[,1])
      out<-as.numeric((max(o)-o+1)<thresh)
      for(i in 2:dim(m)[2]){      
         o<-order(m[,i])
         out<-cbind(out,as.numeric((max(o)-o+1)<thresh))
      }
   }
   if(leq==TRUE)
      out<-1-out
   if(is.null(dim(out))!=is.null(dim(m)))
      out<-array(out,dim=dim(m))
   else
      if(dim(out)!=dim(m))
         out<-array(out,dim=dim(m))
   out
}


#gt - "Graph transpose"; transposition of one or more networks
gt<-function(x, return.as.edgelist=FALSE){
  if(return.as.edgelist){
    #Pre-process the raw input
    x<-as.edgelist.sna(x)
    if(is.list(x))
      return(lapply(x,gt,return.as.edgelist=TRUE))
    #End pre-processing
    n<-attr(x,"n")
    vnames<-attr(x,"vnames")
    bipartite<-attr(x,"bipartite")
    x<-x[,c(2,1,3)]
    attr(x,"n")<-n
    attr(x,"vnames")<-vnames
    attr(x,"bipartite")<-bipartite
    x
  }else{
    #Pre-process the raw input
    x<-as.sociomatrix.sna(x)
    if(is.list(x))
      return(lapply(x,gt,return.as.edgelist=FALSE))
    #End pre-processing
    if(length(dim(x))==3){
      aperm(x,c(1,3,2))
    }else
      t(x)
  }
}


#gvectorize - Vectorization of adjacency matrices
gvectorize<-function(mats,mode="digraph",diag=FALSE,censor.as.na=TRUE){
   #Pre-process the raw input
   mats<-as.sociomatrix.sna(mats)
   if(is.list(mats))
     return(lapply(mats,gvectorize,mode=mode,diag=diag, censor.as.na=censor.as.na))
   #End pre-processing
   #Build the input data structures
   if(length(dim(mats))>2){
      m<-dim(mats)[1]
      n<-dim(mats)[2]
      n<-dim(mats)[3]
      d<-mats
   }else{
      m<-1
      n<-dim(mats)[1]
      o<-dim(mats)[2]
      d<-array(dim=c(1,n,o))
      d[1,,]<-mats
   }
   #If using NA censoring, turn unused parts of the matrices to NAs and vectorize
   if(censor.as.na){
      if(mode=="graph")
         d<-upper.tri.remove(d)
      if(!diag)
         d<-diag.remove(d)
      out<-apply(d,1,as.vector)
   }else{   #Otherwise, vectorize only the useful parts
      if(mode=="graph")
        mask<-apply(d,1,lower.tri,diag=diag)
      else{
        if(diag)
          mask<-matrix(TRUE,nrow=dim(d)[2]*dim(d)[3],ncol=dim(d)[1])
        else
          mask<-apply(d,1,function(z){diag(NROW(z))==0})
      }
      out<-apply(d,1,as.vector)
      if(m==1)
         out<-out[mask]
      else
         out<-matrix(out[mask],ncol=m)      
   }
   out
}


#interval.graph - Construct one or more interval graphs (and exchangeability 
#vectors) from a set of spells
interval.graph<-function(slist,type="simple",diag=FALSE){
   #Note that each slice of slist must have one spell per row, with col 1 containing the spell type,
   #col 2 containing the spell onset, and col 3 containing the spell termination.  If there are multiple
   #slices present, they must be indexed by the first dimension of the array.
   #First, the preliminaries
   o<-list()
   m<-stackcount(slist)          #Get the number of stacks
   if(m==1){
      d<-array(dim=c(m,dim(slist)[1],dim(slist)[2]))
      d[1,,]<-slist
   }else
      d<-slist
   ns<-dim(d)[2]                     #Get the number of spells
   o$exchange.list<-d[,,1]   #Exchange list is just the vector of spell types
   #Now, for the graph itself...
   o$graph<-array(dim=c(m,ns,ns))
   for(i in 1:ns)
      for(j in 1:ns)
         o$graph[,i,j]<-switch(type,
            simple=as.numeric((d[,i,2]<=d[,j,3])&(d[,i,3]>=d[,j,2])),  #"Start before the end, end after the beginning"
            overlap=pmax(pmin(d[,i,3],d[,j,3])-pmax(d[,i,2],d[,j,2]),0),
            fracxy=pmax(pmin(d[,i,3],d[,j,3])-pmax(d[,i,2],d[,j,2]),0)/(d[,i,3]-d[,i,2]),
            fracyx=pmax(pmin(d[,i,3],d[,j,3])-pmax(d[,i,2],d[,j,2]),0)/(d[,j,3]-d[,j,2]),
            jntfrac=2*pmax(pmin(d[,i,3],d[,j,3])-pmax(d[,i,2],d[,j,2]),0)/(d[,i,3]-d[,i,2]+d[,j,3]-d[,j,2])
         )
   #Patch up those loose ends.
   if(m==1)
      o$graph<-o$graph[1,,]
   if(!diag)
      o$graph<-diag.remove(o$graph,remove.val=0)
   #Return the data structure
   o
}


#is.edgelist.sna - check to see if a data object is an sna edgelist
is.edgelist.sna<-function(x){
  if(class(x)=="list")
    return(sapply(x,is.edgelist.sna))
  if(!(class(x)%in%c("matrix","array")))
    FALSE
  else if(length(dim(x))!=2)
    FALSE
  else if(dim(x)[2]!=3)
    FALSE
  else if(is.null(attr(x,"n")))
    FALSE
  else
    TRUE
}


#lower.tri.remove - NA the lower triangles of adjacency matrices in a graph 
#stack
lower.tri.remove<-function(dat,remove.val=NA){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(lapply(dat,lower.tri.remove,val=remove.val))
   #End pre-processing
   if(length(dim(dat))>2){
      d<-dat
      for(i in 1:dim(dat)[1]){
         temp<-d[i,,]
         temp[lower.tri(temp,diag=FALSE)]<-remove.val
         d[i,,]<-temp
      }
   }
   else{
      d<-dat
      d[lower.tri(d,diag=FALSE)]<-remove.val
   }   
   d
}


#make.stochastic - Make a graph stack row, column, or row-column stochastic
make.stochastic<-function(dat,mode="rowcol",tol=0.005,maxiter=prod(dim(dat))*100,anneal.decay=0.01,errpow=1){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(lapply(dat,make.stochastic,mode=mode,tol=tol,maxiter=maxiter, anneal.decay=anneal.decay,errpow=errpow))
   #End pre-processing
   #Organize the data
   m<-stackcount(dat)
   if(m==1){
      n<-dim(dat)[1]
      o<-dim(dat)[2]
      d<-array(dim=c(m,n,o))
      d[1,,]<-dat
   }else{
      n<-dim(dat)[2]
      o<-dim(dat)[3]
      d<-dat
   }
   #Stochasticize
   if(mode=="row"){
      for(i in 1:m)
         d[i,,]<-d[i,,]/t(sapply(apply(d[i,,],1,sum),rep,o))
   }else if(mode=="col"){
      for(i in 1:m)
         d[i,,]<-d[i,,]/sapply(apply(d[i,,],2,sum),rep,n)
   }else if(mode=="rowcol"){
      for(i in 1:m){
         f<-d[i,,]/t(sapply(apply(d[i,,],1,sum),rep,o))   #Seed with the row-stochastic form
         f<-f/sapply(apply(f,2,sum),rep,n)   #Col-stochasticize for good measure (sometimes this works)
         edgelist<-cbind(rep(1:n,o),rep(1:o,rep(n,o)))
         edgelist<-edgelist[d[i,,][edgelist]>0,]   #Skip edges which are forced to be zero-valued
         err<-sum(abs(apply(f,2,sum)-rep(1,o))^errpow,abs(apply(f,1,sum)-rep(1,n))^errpow)
         iter<-0
         while((err>(n+o)*tol)&(iter<maxiter)){  #Right now, use an annealer to find an approximate solution
            edge<-sample(1:dim(edgelist)[1],1)
            x<-edgelist[edge,1]
            y<-edgelist[edge,2]
            draw<-max(0,min(rnorm(1,f[x,y],d[i,x,y]/10),d[i,x,y]))
            nerr<-err-abs(sum(f[x,])-1)^errpow-abs(sum(f[,y])-1)^errpow+abs(sum(f[x,][-y])+draw-1)^errpow+abs(sum(f[,y][-x])+draw-1)^errpow
            if((nerr<err)|(runif(1,0,1)<exp(-anneal.decay*iter))){
               f[x,y]<-draw
               err<-nerr
            }
            iter<-iter+1
         }
         d[i,,]<-f
         if(err>(n+o)*tol)
            warning(paste("Annealer unable to reduce total error below apx",round(err,digits=7),"in matrix",i,". Hope that's OK....\n"))
      }
   }else if(mode=="total"){
         for(i in 1:m)
            d[i,,]<-d[i,,]/sum(d[i,,])
   }
   #Patch NaN values
   d[is.nan(d)]<-0
   #Return the output
   if(m==1)
      d[1,,]
   else
      d
}


#nties - Find the number of ties in a given graph or stack
nties<- function(dat,mode="digraph",diag=FALSE){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(lapply(dat,nties,mode=mode,diag=diag))
   #End pre-processing
   #Did someone send us a stack?
   if(length(dim(dat))>2)
      shiftit<-1
   else
      shiftit<-0
   #Get size params
   n<-dim(dat)[1+shiftit]
   m<-dim(dat)[2+shiftit]
   #Sanity check for hypergraphs
   if(mode=="hgraph")
      diag<-TRUE
   #Get the initial count
   count<-switch(mode,
      digraph = n*n,
      graph = (n*n-n)/2+n,
      hgraph = n*m
   )
   #Modify for diag, if needed
   if(!diag)
      count<-count-n
   #Return the needed info
   if(shiftit)
      rep(count,dim(dat)[1])
   else
      count                
}


#sr2css - Convert a row-wise self-report matrix to a CSS matrix with missing 
#observations.
sr2css<-function(net){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(net)
   if(is.list(net))
     return(lapply(net))
   #End pre-processing
   n<-dim(net)[1]
   css<-array(dim=c(n,n,n))
   for(i in 1:n){
      css[i,,]<-NA
      css[i,i,]<-net[i,]
   }
   css
}


#stackcount -How many matrices in a given stack?
stackcount<-function(d){
   #Pre-process the raw input
   d<-as.edgelist.sna(d)
   #End pre-processing
   if(is.list(d))
     length(d)
   else
     1
}


#symmetrize - Convert a graph or graph stack to a symmetric form.  Current rules
#for symmetrizing include "upper" and "lower" diagonals, "weak" connectedness 
#rule, and a "strong" connectedness rule.  If return.as.edgelist=TRUE, the
#data is processed and returned in sna edgelist form.
symmetrize<-function(mats,rule="weak",return.as.edgelist=FALSE){
   if(!return.as.edgelist){                 #Adjacency matrix form
     #Pre-process the raw input
     mats<-as.sociomatrix.sna(mats)
     if(is.list(mats))
       return(lapply(mats,symmetrize,rule=rule, return.as.edgelist=return.as.edgelist))
     #End pre-processing
     #Build the input data structures
     if(length(dim(mats))>2){
        m<-dim(mats)[1]
        n<-dim(mats)[2]
        o<-dim(mats)[3]
        d<-mats
     }else{
        m<-1
        n<-dim(mats)[1]
        o<-dim(mats)[2]
        d<-array(dim=c(1,n,o))
        d[1,,]<-mats
     }
     #Apply the symmetry rule
     for(i in 1:m){
        if(rule=="upper"){
           d[i,,][lower.tri(d[i,,])]<-t(d[i,,])[lower.tri(d[i,,])]
        }else if(rule=="lower"){
           d[i,,][upper.tri(d[i,,])]<-t(d[i,,])[upper.tri(d[i,,])]
        }else if(rule=="weak"){
           d[i,,]<-matrix(as.numeric(d[i,,]|t(d[i,,])),nrow=n,ncol=o)
        }else if(rule=="strong"){
           d[i,,]<-matrix(as.numeric(d[i,,]&t(d[i,,])),nrow=n,ncol=o)
        }
     }
     #Return the symmetrized matrix
     if(m==1)
        out<-d[1,,]
     else
        out<-d
     out
   }else{                                  #Edgelist matrix form
     #Pre-process the raw input
     mats<-as.edgelist.sna(mats)
     if(is.list(mats))
       return(lapply(mats,symmetrize,rule=rule, return.as.edgelist=return.as.edgelist))
     #End pre-processing
     n<-attr(mats,"n")
     vn<-attr(mats,"vnames")
     bip<-attr(mats,"bipartite")
     if(!is.null(bip))
       return(mats)                   #Return unaltered if bipartite
     #Apply the symmetry rule
     if(rule=="upper"){
       loops<-mats[mats[,1]==mats[,2],,drop=FALSE]
       upedge<-mats[mats[,1]<mats[,2],,drop=FALSE]
       mats<-rbind(upedge,upedge[,c(2,1,3)],loops)
     }else if(rule=="lower"){
       loops<-mats[mats[,1]==mats[,2],,drop=FALSE]
       loedge<-mats[mats[,1]>mats[,2],,drop=FALSE]
       mats<-rbind(loedge,loedge[,c(2,1,3)],loops)
     }else if(rule=="weak"){
       isloop<-mats[,1]==mats[,2]
       loops<-mats[isloop,,drop=FALSE]
       mats<-mats[!isloop,,drop=FALSE]
       dc<-.C("dyadcode_R",as.double(mats),as.integer(n),as.integer(NROW(mats)), dc=as.double(rep(0,NROW(mats))),PACKAGE="sna",NAOK=TRUE)$dc
       isdup<-duplicated(dc)
       mats<-mats[!isdup,,drop=FALSE]
       mats<-rbind(mats,mats[,c(2,1,3)],loops)
     }else if(rule=="strong"){
       isloop<-mats[,1]==mats[,2]
       loops<-mats[isloop,,drop=FALSE]
       mats<-mats[!isloop,,drop=FALSE]
       dc<-.C("dyadcode_R",as.double(mats),as.integer(n),as.integer(NROW(mats)), dc=as.double(rep(0,NROW(mats))),PACKAGE="sna",NAOK=TRUE)$dc
       isdup<-duplicated(dc)
       mats<-mats[isdup,,drop=FALSE]
       mats<-rbind(mats,mats[,c(2,1,3)],loops)
     }
     #Patch up the attributes and return
     attr(mats,"n")<-n
     attr(mats,"vnames")<-vn
     mats
   }
}


#upper.tri.remove - NA the upper triangles of adjacency matrices in a graph 
#stack
upper.tri.remove<-function(dat,remove.val=NA){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(lapply(dat,upper.tri.remove,remove.val=remove.val))
   #End pre-processing
   if(length(dim(dat))>2){
      d<-dat
      for(i in 1:dim(dat)[1]){
         temp<-d[i,,]
         temp[upper.tri(temp,diag=FALSE)]<-remove.val
         d[i,,]<-temp
      }
   }
   else{
      d<-dat
      d[upper.tri(d,diag=FALSE)]<-remove.val
   }   
   d
}


