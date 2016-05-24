######################################################################
#
# roles.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 4/8/10
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains functions relating to role analysis.
#
# Contents:
#   blockmodel
#   blockmodel.expand
#   equiv.clust
#   plot.blockmodel
#   plot.equiv.clust
#   print.blockmodel
#   print.equiv.clust
#   print.summary.blockmodel
#   redist
#   sedist
#   summary.blockmodel
#
######################################################################


#blockmodel - Generate blockmodels based on partitions of network positions
blockmodel<-function(dat,ec,k=NULL,h=NULL,block.content="density",plabels=NULL,glabels=NULL,rlabels=NULL,mode="digraph",diag=FALSE){
   #First, extract the blocks
   if(class(ec)=="equiv.clust")
     b<-cutree(ec$cluster,k,h)
   else if(class(ec)=="hclust")
     b<-cutree(ec,k,h)
   else
     b<-ec
   #Prepare the data
   dat<-as.sociomatrix.sna(dat,simplify=TRUE)
   if(is.list(dat))
     stop("Blockmodel requires input graphs to be of identical order.")
   n<-dim(dat)[2]
   if(length(dim(dat))>2)
      d<-dat
   else{
      d<-array(dim=c(1,n,n))
      d[1,,]<-dat
   }
   if(!diag)
      d<-diag.remove(d)
   #Get labels
   if(is.null(plabels)){
     if(class(ec)=="equiv.clust")
       plabels<-ec$plabels
     else
       plabels<-1:length(b)
   }
   if(is.null(glabels)){
     if(class(ec)=="equiv.clust")
       glabels<-ec$glabels
     else
       glabels<-1:length(b)
   }
   #Now, construct a model
   rn<-max(b)
   rm<-dim(d)[1]
   if(is.null(rlabels))      #Add labels for roles if needed
      rlabels<-paste("Block",1:rn)
   bm<-array(dim=c(rm,rn,rn))
   for(i in 1:rm)
      for(j in 1:rn)
         for(k in 1:rn){
            if(block.content=="density")
               bm[i,j,k]<-mean(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
            else if(block.content=="meanrowsum"){
               bm[i,j,k]<-mean(apply(d[i,b==j,b==k,drop=FALSE],2,sum,na.rm=TRUE))
            }else if(block.content=="meancolsum"){
               bm[i,j,k]<-mean(apply(d[i,b==j,b==k,drop=FALSE],3,sum,na.rm=TRUE))
            }else if(block.content=="sum"){
               bm[i,j,k]<-sum(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
            }else if(block.content=="median"){
               bm[i,j,k]<-median(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
            }else if(block.content=="min"){
               bm[i,j,k]<-min(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
            }else if(block.content=="max"){
               bm[i,j,k]<-max(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
            }else if(block.content=="types"){
               temp<-mean(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
               if(is.nan(temp))    #Is this a nan block (due to having only one actor)?
                  bm[i,j,k]<-"NA"
               else if(temp==0)    #Is this a null block?
                  bm[i,j,k]<-"null"
               else if(temp==1)   #How about a complete block?
                  bm[i,j,k]<-"complete"
               else if(all(apply(d[i,b==j,b==k,drop=FALSE],2,sum,na.rm=TRUE)>0,apply(d[i,b==j,b==k,drop=FALSE],3,sum,na.rm=TRUE)>0))
                  bm[i,j,k]<-"1 covered"   #1 covered block
               else if(all(apply(d[i,b==j,b==k,drop=FALSE],2,sum,na.rm=TRUE)>0))
                  bm[i,j,k]<-"1 row-covered"   #1 row-covered block
               else if(all(apply(d[i,b==j,b==k,drop=FALSE],3,sum,na.rm=TRUE)>0))
                  bm[i,j,k]<-"1 col-covered"   #1 col-covered block
               else
                  bm[i,j,k]<-"other"   #other block
            }
         }
   #Prepare the output object
   if(class(ec)=="equiv.clust")
     pord<-ec$cluster$order
   else if(class(ec)=="hclust")
     pord<-ec$order
   else
     pord<-order(ec)
   o<-list()
   o$block.membership<-b[pord]
   o$order.vector<-pord
   o$block.content<-block.content
   if(length(dim(dat))>2){
      o$blocked.data<-dat[,pord,pord]
      dimnames(o$blocked.data)<-list(glabels,plabels[pord],plabels[pord])
   }else{
      o$blocked.data<-dat[pord,pord]      
      dimnames(o$blocked.data)<-list(plabels[pord],plabels[pord])
   }
   if(dim(bm)[1]==1){
      o$block.model<-bm[1,,]
      rownames(o$block.model)<-rlabels
      colnames(o$block.model)<-rlabels
   }else{
      o$block.model<-bm
      dimnames(o$block.model)<-list(glabels,rlabels,rlabels)
   }
   o$plabels<-plabels[pord]
   o$glabels<-glabels
   o$rlabels<-rlabels
   o$cluster.method<-switch(class(ec),
     equiv.clust=ec$cluster.method,
     hclust=ec$method,
     "Prespecified"
   )
   o$equiv.fun<-switch(class(ec),
     equiv.clust=ec$equiv.fun,
     "None"
   )
   o$equiv.metric<-switch(class(ec),
     equiv.clust=ec$metric,
     "None"
   )
   class(o)<-"blockmodel"
   o   
}


#blockmodel.expand - Generate a graph (or stack) from a given blockmodel using 
#particular expansion rules
blockmodel.expand<-function(b,ev,mode="digraph",diag=FALSE){
   #First, get some useful parameters and such
   en<-sum(ev)
   el<-length(ev)
   bn<-max(b$block.membership)
   bm<-stackcount(b$block.model)
   if(bm>1)
      block.model<-b$block.model
   else{
      block.model<-array(dim=c(1,bn,bn))
      block.model[1,,]<-b$block.model
   }
   #Now, perform the expansion)
   expanded<-array(dim=c(bm,en,en))
   for(i in 1:bm){
      if(b$block.content=="density"){
         tp<-matrix(nrow=en,ncol=en)
         for(j in 1:el)
            for(k in 1:el)
               tp[(cumsum(ev)[j]-ev[j]+1):(cumsum(ev)[j]),(cumsum(ev)[k]-ev[k]+1):(cumsum(ev)[k])]<-block.model[i,j,k]
         tp[is.na(tp)|is.nan(tp)]<-0   #Fill in any NA or NaN blocks with zero
         expanded[i,,]<-rgraph(en,1,tprob=tp,mode=mode,diag=diag)
      }else
         stop(paste("\nContent type",b$block.content,"not supported yet.\n"))
   }
   #Return the output data
   if(dim(expanded)[1]>1)
      expanded
   else
      expanded[1,,]
}


#equiv.clust - Find clusters of positions based on an equivalence relation
equiv.clust<-function(dat,g=NULL,equiv.dist=NULL,equiv.fun="sedist",method="hamming",mode="digraph",diag=FALSE,cluster.method="complete",glabels=NULL,plabels=NULL,...){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   #End pre-processing
   #First, find the equivalence distances using the appropriate function and method
   if(is.null(g)){             #Set g to all graphs, if needed
     if(is.list(dat))
       g<-1:length(dat)
     else if(is.array(dat))
       g<-1:dim(dat)[1]
     else
       g<-1
   }     
   if(is.null(equiv.dist)){
     equiv.dist.fun<-match.fun(equiv.fun)
     equiv.dist<-equiv.dist.fun(dat,g=g,method=method,joint.analysis=TRUE, mode=mode,diag=diag,code.diss=TRUE,...)
   }
   #Generate the output object
   o<-list()
   #Produce the hierarchical clustering
   o$cluster<-hclust(as.dist(equiv.dist),method=cluster.method)
   #Generate labels
   if(is.null(glabels)){
     if(is.list(dat))
       glabels<-names(dat)[g]
     else
       glabels<-dimnames(dat)[[1]][g]
   }
   if(is.null(plabels)){
     if(is.list(dat))
       plabels<-dimnames(dat[[g]])[[2]]
     else
       plabels<-dimnames(dat)[[2]]
   }
   #Set the output class and take care of other details
   o$metric<-method
   o$equiv.fun<-equiv.fun
   o$cluster.method<-cluster.method
   if((length(dim(dat))==1)&(length(glabels)>1))
      o$glabels<-glabels[1]
   else
      o$glabels<-glabels
   o$plabels<-plabels
   class(o)<-"equiv.clust"
   #Produce the output
   o
}


#plot.blockmodel - Plotting for blockmodel objects
plot.blockmodel<-function(x,...){
   #Save old settings
   oldpar<-par(no.readonly=TRUE)
   on.exit(par(oldpar))
   #Get new settings from data
   n<-dim(x$blocked.data)[2]
   m<-stackcount(x$blocked.data)
   if(!is.null(x$plabels))
      plab<-x$plabels
   else
      plab<-(1:n)[x$order.vector]
   if(!is.null(x$glabels))
      glab<-x$glabels
   else
      glab<-1:m
   #Now, plot the blocked data
   par(mfrow=c(floor(sqrt(m)),ceiling(m/floor(sqrt(m)))))
   if(m>1)
      for(i in 1:m){
         plot.sociomatrix(x$blocked.data[i,,],labels=list(plab,plab), main=paste("Relation - ",glab[i]),drawlines=FALSE)
         for(j in 2:n)
            if(x$block.membership[j]!=x$block.membership[j-1])
               abline(v=j-0.5,h=j-0.5,lty=3)
      }
   else{
      plot.sociomatrix(x$blocked.data,labels=list(plab,plab), main=paste("Relation - ",glab[1]),drawlines=FALSE)
      for(j in 2:n)
         if(x$block.membership[j]!=x$block.membership[j-1])
            abline(v=j-0.5,h=j-0.5,lty=3)
   }
}


#plot.equiv.clust - Plotting for equivalence clustering objects
plot.equiv.clust<-function(x,labels=NULL,...){
   if(is.null(labels))
      plot(x$cluster,labales=x$labels,...)
   else
      plot(x$cluster,labels=labels,...)
}


#print.blockmodel - Printing for blockmodel objects
print.blockmodel<-function(x,...){
   cat("\nNetwork Blockmodel:\n\n")
   cat("Block membership:\n\n")
   if(is.null(x$plabels))                    #Get position labels
      plab<-(1:length(x$block.membership))[x$order.vector]
   else
      plab<-x$plabels
   temp<-matrix(x$block.membership,nrow=1)
   dimnames(temp)<-list("",plab)
   print(temp[1,order(x$order.vector)])  #Print in original order
   cat("\nReduced form blockmodel:\n\n")
   if(length(dim(x$block.model))>2){
      for(i in 1:dim(x$block.model)[1]){
         temp<-x$block.model[i,,]
         dimnames(temp)<-list(x$rlabels,x$rlabels)
         cat("\t",x$glabels[i],"\n") 
         print(temp)
         cat("\n")
      }
   }else{
         temp<-x$block.model
         dimnames(temp)<-list(x$rlabels,x$rlabels)
         cat("\t",x$glabels,"\n") 
         print(temp)
   }
}


#print.equiv.clust - Printing for equiv.clust objects
print.equiv.clust<-function(x, ...){
  cat("Position Clustering:\n\n")
  cat("\tEquivalence function:",x$equiv.fun,"\n")
  cat("\tEquivalence metric:",x$metric,"\n")
  cat("\tCluster method:",x$cluster.method,"\n")
  cat("\tGraph order:",length(x$cluster$order),"\n\n")
}


#print.summary.blockmodel - Printing for blockmodel summary objects
print.summary.blockmodel<-function(x,...){
   cat("\nNetwork Blockmodel:\n\n")

   cat("\nGeneral information:\n\n")
   cat("\tEquivalence function: ",x$equiv.fun,"\n")
   cat("\tEquivalence metric: ",x$equiv.metric,"\n")
   cat("\tClustering method: ",x$cluster.method,"\n")
   cat("\tBlockmodel content: ",x$block.content,"\n")

   cat("\n\nBlock membership by actor:\n\n")
   if(is.null(x$plabels))                    #Get position labels
      plab<-(1:length(x$block.membership))[x$order.vector]
   else
      plab<-x$plabels
   temp<-matrix(x$block.membership,nrow=1)
   dimnames(temp)<-list("",plab)
   print(temp[1,order(x$order.vector)])  #Print in original order

   cat("\n\nBlock membership by block:\n\n")
   for(i in 1:max(x$block.membership))
      cat("\t",x$rlabels[i],":",plab[x$block.membership==i],"\n")
   
   cat("\n\nReduced form blockmodel:\n\n")
   if(length(dim(x$block.model))>2){
      for(i in 1:dim(x$block.model)[1]){
         temp<-x$block.model[i,,]
         dimnames(temp)<-list(x$rlabels,x$rlabels)
         cat("\t",x$glabels[i],"\n") 
         print(temp)
         cat("\n")
      }
   }else{
         temp<-x$block.model
         dimnames(temp)<-list(x$rlabels,x$rlabels)
         cat("\t",x$glabels,"\n") 
         print(temp)
   }

   cat("\n\nBlocked data:\n\n")
   if(length(dim(x$block.model))>2){
      for(i in 1:dim(x$block.model)[1]){
         temp<-x$blocked.data[i,,]
         dimnames(temp)<-list(plab,plab)
         cat("\t",x$glabels[i],"\n") 
         print(temp)
         cat("\n")
      }
   }else{
         temp<-x$blocked.data
         dimnames(temp)<-list(plab,plab)
         cat("\t",x$glabels,"\n") 
         print(temp)
   }

}


#redist - Find a matrix of distances between positions based on regular 
#equivalence
redist<-function(dat, g=NULL, method=c("catrege"), mode="digraph", diag=FALSE, seed.partition=NULL, code.diss=TRUE, ...){
  #Internal function to compute neighborhoods for CATREGE
  neighb<-function(){
    nmat<-array(0,dim=c(r,n,n))
    for(i in 1:n)
      for(j in 1:n)
        if(d[i,j]>0)
          nmat[d[i,j],i,part1[j]]<-TRUE
    nmat
  }
  #Prep the data
  dat<-as.sociomatrix.sna(dat,simplify=TRUE)
  if(is.list(dat))
    stop("redist requires input graphs to be of identical order.")
  if(is.null(g))
    g<-1:dim(dat)[1]
  if(length(dim(dat)) > 2) {
    n <- dim(dat)[2]
    m <- length(g)
    d <- dat[g, , ]
  }else{
    n <- dim(dat)[2]
    m <- 1
    d <- array(dim = c(m, n, n))
    d[1, , ] <- dat
  }
  if(mode == "graph")
    d <- symmetrize(d)
  if(m==1)
    d<-array(d,dim=c(1,n,n))
  if (!diag)
    d <- diag.remove(d,0) #Currently, treat as zeros
  #Build the categorical matrix
  da<-array(dim=c(2*m,n,n))   #First, perform symmetric interleaving
  for(i in 1:m){
    da[i*2-1,,]<-d[i,,]
    da[i*2,,]<-t(d[i,,])
  }
  d<-apply(da,c(2,3),paste,collapse=" ")  #Convert to strings
  vals<-apply(sapply((1:2^(2*m))-1,function(z){(z%/%2^((1:(2*m))-1))%%2}),2, paste,collapse=" ")  #Obtain all possible strings
  r<-length(vals)-1   #Non-null values
  d<-apply(d,c(1,2),match,vals)-1   #Replace with numeric values
#  print(d[1:15,1:15])
#  vals<-sort(unique(as.vector(d)))             #Obtain unique values
#  print(vals)
#  r<-length(vals-1)                        #Get number of unique values
#  vals0<-grep("NA",vals)                 #Fix zeros
#  print(vals0)
#  vals0<-c(vals0,((1:r)[-vals0])[as.numeric(gsub(" ","",vals[-vals0]))==0]) 
#  print(vals0)
#  d<-apply(d,c(1,2),match,vals)          #Replace vals with numerics
#  d[d%in%vals0]<-0                       #Set zeros
  #Compute the equivalence
  if(match.arg(method)=="catrege"){
    outpart<-vector()
    if(is.null(seed.partition))
      part1<-rep(1,n)   #Create initial partition memberships
    else
      part1<-seed.partition
    flag<-TRUE
    while(flag){
      nmat<-neighb()  #Compute neighborhoods, using current partition
      outpart<-rbind(outpart,part1)
      flag<-FALSE     #Set change flag
      part2<-1:n
      for(i in 2:n)
        for(j in 1:(i-1))
          if(part1[i]==part1[j]){
            if(all(nmat[,i,]==nmat[,j,]))
              part2[i]<-part2[j]
            else
              flag<-TRUE
          }
      part1<-part2
    }
    eq<-matrix(0,n,n)
    for(i in 1:n)
      for(j in 1:n)
        eq[i,j]<-max((1:NROW(outpart))[outpart[,i]==outpart[,j]])
  }
  #Transform and rescale to distance form if required
  if(!code.diss)
    eq
  else{
    if(max(eq)==min(eq))
      matrix(0,NROW(eq),NCOL(eq))
    else
      (max(eq)-eq)/(max(eq)-min(eq))
  }
}


#sedist - Find a matrix of distances between positions based on structural 
#equivalence
sedist<-function(dat,g=c(1:dim(dat)[1]),method="hamming",joint.analysis=FALSE,mode="digraph",diag=FALSE,code.diss=FALSE){
   #First, prepare the data
   dat<-as.sociomatrix.sna(dat,simplify=TRUE)
   if(is.list(dat))
     stop("sedist requires input graphs to be of identical order.")
   if(length(dim(dat))>2){
      n<-dim(dat)[2]
      m<-length(g)
      d<-dat[g,,,drop=FALSE]
   }else{
      n<-dim(dat)[2]
      m<-1
      d<-array(dim=c(m,n,n))
      d[1,,]<-dat
   }
   if(!diag)
      d<-diag.remove(d)
   #Are we conducting a joint analysis?
   if(joint.analysis){
      o<-array(dim=c(1,n,n))
      #Build the data matrix
      v<-vector()
      for(i in 1:n)
         v<-cbind(v,c(as.vector(d[,i,]),as.vector(d[,,i])))
      #Proceed by method
      if(method=="correlation"){
         o[1,,]<-cor(v,use="pairwise")
         #Reverse code?
         if(code.diss)
            o<--o
      }else if(method=="euclidean"){
         for(i in 1:n)
            for(j in 1:n)
               o[1,i,j]<-sqrt(sum((v[,i]-v[,j])^2,na.rm=TRUE))         
      }else if(method=="hamming"){
         for(i in 1:n)
            for(j in 1:n)
               o[1,i,j]<-sum(abs(v[,i]-v[,j]),na.rm=TRUE)
      }else if(method=="gamma"){
         for(i in 1:n)
            for(j in 1:n){
               concord<-sum(as.numeric(v[,i]==v[,j]),na.rm=TRUE)
               discord<-sum(as.numeric(v[,i]!=v[,j]),na.rm=TRUE)
               o[1,i,j]<-(concord-discord)/(concord+discord)
            }                  
         #Reverse code?
         if(code.diss)
            o<--o
      }else if(method=="exact"){
         for(i in 1:n)
            for(j in 1:n)
               o[1,i,j]<-as.numeric(any(v[!(is.na(v[,i])|is.na(v[,j])),i]!=v[!(is.na(v[,i])|is.na(v[,j])),j]))
      }
   }else{  #Analyze each graph seperately
      o<-array(dim=c(m,n,n))
      for(k in 1:m){
         #Build the data matrix
         v<-vector()
         for(i in 1:n)
            v<-cbind(v,c(as.vector(d[k,i,]),as.vector(d[k,,i])))
         #Proceed by method
         if(method=="correlation"){
            o[k,,]<-cor(v,use="pairwise")
            o[k,,][is.na(o[k,,])]<-0
            #Reverse code?
            if(code.diss)
               o[k,,]<--o[k,,]
         }else if(method=="euclidean"){
            for(i in 1:n)
               for(j in 1:n)
                  o[k,i,j]<-sqrt(sum((v[,i]-v[,j])^2,na.rm=TRUE))         
         }else if(method=="hamming"){
            for(i in 1:n)
               for(j in 1:n)
                  o[k,i,j]<-sum(abs(v[,i]-v[,j]),na.rm=TRUE)
         }else if(method=="gamma"){
            for(i in 1:n)
               for(j in 1:n){
                  concord<-sum(as.numeric(v[,i]==v[,j]),na.rm=TRUE)
                  discord<-sum(as.numeric(v[,i]!=v[,j]),na.rm=TRUE)
                  o[k,i,j]<-(concord-discord)/(concord+discord)
               }                  
            #Reverse code?
            if(code.diss)
               o[k,,]<--o[k,,]
         }else if(method=="exact"){
            for(i in 1:n)
               for(j in 1:n)
                  o[k,i,j]<-as.numeric(any(v[!(is.na(v[,i])|is.na(v[,j])),i]!=v[!(is.na(v[,i])|is.na(v[,j])),j]))
         }
      }
   }
   #Dump the output
   if(dim(o)[1]==1)
      as.matrix(o[1,,])
   else
      o
}


#summary.blockmodel - Detailed printing for blockmodel objects
summary.blockmodel<-function(object, ...){
   o<-object
   class(o)<-"summary.blockmodel"
   o
}
