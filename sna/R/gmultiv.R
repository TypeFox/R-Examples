######################################################################
#
# gmultiv.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 4/30/09
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines associated with multivariate analysis of
# graph sets.
#
# Contents:
#   centralgraph
#   gclust.boxstats
#   gclust.centralgraph
#   gcor
#   gcov
#   gdist.plotdiff
#   gdist.plotstats
#   gscor
#   gscov
#   hdist
#   sdmat
#   structdist
#
######################################################################


#centralgraph - Find the central graph of a graph stack
centralgraph<-function(dat,normalize=FALSE){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     stop("Identical graph orders required in centralgraph.")
   #End pre-processing
   #Check to see if someone foolishly called this with one graph
   if(length(dim(dat))==2)
      out<-dat
   else{
      if(normalize)
         out<-apply(dat,c(2,3),mean,na.rm=TRUE)
      else
         out<-matrix(data=as.numeric(apply(dat,c(2,3),mean,na.rm=TRUE)>=0.5), nrow=dim(dat)[2],ncol=dim(dat)[2])
   }
   out
}


#gclust.boxstats - Plot statistics associated with clusters
gclust.boxstats<-function(h,k,meas,...){
   #h must be an hclust object, k the number of groups, and meas the group measurement vector
   out<-matrix(nrow=length(meas),ncol=k)
   gmat<-matrix(nrow=length(meas),ncol=2)
   gmat[,1]<-c(1:length(meas))
   gmat[,2]<-cutree(h,k=k)
   for(i in 1:k){
      out[1:length(meas[gmat[gmat[,2]==i,1]]),i]<-meas[gmat[gmat[,2]==i,1]]
   }
   boxplot(data.frame(out),...)
}


#gclust.centralgraph - Get central graphs associated with clusters
gclust.centralgraph<-function(h,k,dat,...){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     stop("Identical graph orders required in gclust.centralgraph.")
   #End pre-processing
   #h must be an hclust object, k the number of groups, and dat a collection
   #of graphs (with identical order)
   out<-array(dim=c(k,dim(dat)[2],dim(dat)[2]))
   gmat<-matrix(nrow=dim(dat)[1],ncol=2)
   gmat[,1]<-c(1:dim(dat)[1])
   gmat[,2]<-cutree(h,k=k)
   for(i in 1:k)
      out[i,,]<-centralgraph(dat[gmat[gmat[,2]==i,1],,],...)
   out
}


#gcor - Correlation between two or more graphs.
gcor<-function(dat,dat2=NULL,g1=NULL,g2=NULL,diag=FALSE,mode="digraph"){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     stop("Identical graph orders required in gcor.")
   if(!is.null(dat2)){
     dat2<-as.sociomatrix.sna(dat2)
     if(is.list(dat2))
       stop("Identical graph orders required in gcor.")
   }
   #End pre-processing
   #Collect data and various parameters
   if(is.null(g1))
     g1<-1:dim(dat)[1]
   if(is.null(g2))
     g2<-1:dim(dat)[1]
   if(!is.null(dat2)){
      if(length(dim(dat))>2)
         temp1<-dat
      else{
         temp1<-array(dim=c(1,dim(dat)[2],dim(dat)[2]))
         temp1[1,,]<-dat
      }
      if(length(dim(dat2))>2)
         temp2<-dat2
      else{
         temp2<-array(dim=c(1,dim(dat2)[2],dim(dat2)[2]))
         temp2[1,,]<-dat2
      }
      if(dim(temp1)[2]>dim(temp2)[2])
         temp2<-add.isolates(temp2,dim(temp1)[2]-dim(temp2)[2])
      if(dim(temp2)[2]>dim(temp1)[2])
         temp1<-add.isolates(temp1,dim(temp2)[2]-dim(temp1)[2])
      n<-dim(temp1)[2]
      gn<-dim(temp1)[1]+dim(temp2)[1]
      gn1<-dim(temp1)[1]
      gn2<-dim(temp2)[1]
      d<-array(dim=c(gn,n,n))
      d[1:gn1,,]<-temp1
      d[(gn1+1):(gn2+gn1),,]<-temp2
      g1<-1:gn1
      g2<-(gn1+1):(gn1+gn2)
   }else{
      d<-dat
      n<-dim(dat)[2]
      gn<-dim(dat)[1]
      gn1<-length(g1)
      gn2<-length(g2)
   }
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #Compute the graph correlation matrix
   gd<-matrix(nrow=gn1,ncol=gn2)
   rownames(gd)<-g1
   colnames(gd)<-g2
   for(i in 1:gn1)
      for(j in 1:gn2)
         gd[i,j]<-cor(as.vector(d[g1[i],,]),as.vector(d[g2[j],,]), use="complete.obs")
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      gd[1,1]
   else
      gd
}


#gcov - Covariance between two or more graphs.
gcov<-function(dat,dat2=NULL,g1=NULL,g2=NULL,diag=FALSE,mode="digraph"){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     stop("Identical graph orders required in gcov.")
   if(!is.null(dat2)){
     dat2<-as.sociomatrix.sna(dat2)
     if(is.list(dat2))
       stop("Identical graph orders required in gcov.")
   }
   #End pre-processing
   #Collect data and various parameters
   if(is.null(g1))
     g1<-1:dim(dat)[1]
   if(is.null(g2))
     g2<-1:dim(dat)[1]
   if(!is.null(dat2)){
      if(length(dim(dat))>2)
         temp1<-dat
      else{
         temp1<-array(dim=c(1,dim(dat)[2],dim(dat)[2]))
         temp1[1,,]<-dat
      }
      if(length(dim(dat2))>2)
         temp2<-dat2
      else{
         temp2<-array(dim=c(1,dim(dat2)[2],dim(dat2)[2]))
         temp2[1,,]<-dat2
      }
      if(dim(temp1)[2]>dim(temp2)[2])
         temp2<-add.isolates(temp2,dim(temp1)[2]-dim(temp2)[2])
      if(dim(temp2)[2]>dim(temp1)[2])
         temp1<-add.isolates(temp1,dim(temp2)[2]-dim(temp1)[2])
      n<-dim(temp1)[2]
      gn<-dim(temp1)[1]+dim(temp2)[1]
      gn1<-dim(temp1)[1]
      gn2<-dim(temp2)[1]
      d<-array(dim=c(gn,n,n))
      d[1:gn1,,]<-temp1
      d[(gn1+1):(gn2+gn1),,]<-temp2
      g1<-1:gn1
      g2<-(gn1+1):(gn1+gn2)
   }else{
      d<-dat
      n<-dim(dat)[2]
      gn<-dim(dat)[1]
      gn1<-length(g1)
      gn2<-length(g2)
   }
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #Compute the graph covariance matrix
   gd<-matrix(nrow=gn1,ncol=gn2)
   rownames(gd)<-g1
   colnames(gd)<-g2
   for(i in 1:gn1)
      for(j in 1:gn2)
         gd[i,j]<-cov(as.vector(d[g1[i],,]),as.vector(d[g2[j],,]), use="complete.obs")
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      gd[1,1]
   else
      gd
}


#gdist.plotdiff - Plot differences in graph-level statistics against inter-graph distances
gdist.plotdiff<-function(d,meas,method="manhattan",jitter=TRUE,xlab="Inter-Graph Distance",ylab="Measure Distance",lm.line=FALSE,...){
   #Note that d must be a matrix of distances, and meas must be a matrix of measures
   #Create a matrix of differences in graph-level statistics
   md<-dist(meas,method=method)
   #Vectorize
   dv<-as.vector(as.dist(d))
   mdv<-as.vector(md)
   if(jitter){
      dv<-jitter(dv)
      mdv<-jitter(mdv)
   }
   #Plot the results
   plot(dv,mdv,xlab=xlab,ylab=ylab,...)   
   #If needed, add a line fit
   abline(lm(mdv~dv),col="red")
}


#gdist.plotstats - Plot statistics associated with graphs against (projected) inter-graph distances
gdist.plotstats<-function(d,meas,siz.lim=c(0,0.15),rescale="quantile",display.scale="radius",display.type="circleray",cex=0.5,pch=1,labels=NULL,pos=1,labels.cex=1,legend=NULL,legend.xy=NULL,legend.cex=1,...){
   #d must be a matrix of distances (e.g., from structdist or hdist) and meas a matrix or vector of measures
   #Perform an MDS on the distances
   xy<-cmdscale(as.dist(d))
   n<-dim(xy)[1]
   #Adjust and rescale the measure(s) prior to display
   if(is.null(dim(meas)))
      m<-matrix(meas,ncol=1)
   else
      m<-meas
   nm<-dim(m)[2]
   if(rescale=="quantile"){   #Rescale by (empirical) quantiles (ordinal rescaling)
      m<-apply(m,2,order)
      m<-sweep(m,2,apply(m,2,min))
      m<-sweep(m,2,apply(m,2,max),"/")
   }else if(rescale=="affine"){   #Rescale the measure(s) by affine transformation to the [0,1] interval
      m<-sweep(m,2,apply(m,2,min))
      m<-sweep(m,2,apply(m,2,max),"/")
   }else if(rescale=="normalize"){   #Rescale the measure(s) by their maximum values
      m<-sweep(m,2,apply(m,2,max),"/")
   }
   #Determine how large our drawn symbols are to be
   if(display.scale=="area")  #If we're using area scaling, we need to take square roots
      m<-sqrt(m)
   msize<-m*siz.lim[2]+siz.lim[1]   #Now, express the scaled measures as fractions of the plotting range
   pwid<-max(xy)-min(xy)  #Grab width of plotting range
   msize<-msize*pwid          #Adjust the msize matrix.
   #Plot the data
   plot(xy,xlab=expression(lambda[1]),ylab=expression(lambda[2]),cex=cex,pch=pch,xlim=c(min(xy),max(xy)),ylim=c(min(xy),max(xy)),...)   #Plot the graphs' MDS positions
   if(display.type=="poly"){  #Plot measures using polygons
      for(i in 1:nm){
         for(j in 1:n){
            x<-xy[j,1]+sin(2*pi*((0:nm)/nm))*msize[j,i]
            y<-xy[j,2]+cos(2*pi*((0:nm)/nm))*msize[j,i]
            lines(x,y,col=i)
         }
      }
   }else if(display.type=="circle"){  #Plot measures using circles
      for(i in 1:nm){
         for(j in 1:n){
            x<-xy[j,1]+sin(2*pi*((0:500)/500))*msize[j,i]
            y<-xy[j,2]+cos(2*pi*((0:500)/500))*msize[j,i]
            lines(x,y,col=i)
         }
      }
   }else if(display.type=="ray"){  #Plot measures using rays
      for(i in 1:nm){
         for(j in 1:n){
            lines(c(xy[j,1],xy[j,1]+sin(2*pi*((i-1)/nm))*msize[j,i]),c(xy[j,2],xy[j,2]+cos(2*pi*((i-1)/nm))*msize[j,i]),col=i)
         }
      }      
   }else if(display.type=="polyray"){  #Plot measures using polys and rays
      for(i in 1:nm){
         for(j in 1:n){
            x<-xy[j,1]+sin(2*pi*((0:nm)/nm))*msize[j,i]
            y<-xy[j,2]+cos(2*pi*((0:nm)/nm))*msize[j,i]
            lines(x,y,col=i)
            lines(c(xy[j,1],xy[j,1]+sin(2*pi*((i-1)/nm))*msize[j,i]),c(xy[j,2],xy[j,2]+cos(2*pi*((i-1)/nm))*msize[j,i]),col=i)
         }
      }      
   }else if(display.type=="circleray"){  #Plot measures using circles and rays
      for(i in 1:nm){
         for(j in 1:n){
            x<-xy[j,1]+sin(2*pi*((0:500)/500))*msize[j,i]
            y<-xy[j,2]+cos(2*pi*((0:500)/500))*msize[j,i]
            lines(x,y,col=i)
            lines(c(xy[j,1],xy[j,1]+sin(2*pi*((i-1)/nm))*msize[j,i]),c(xy[j,2],xy[j,2]+cos(2*pi*((i-1)/nm))*msize[j,i]),col=i)
         }
      }      
   }
   #Add labels, if needed
   if(!is.null(labels))
      text(xy[,1],xy[,2],labels,pos=pos,cex=labels.cex)
   #Add legend?
   if(!is.null(legend)){
      if(is.null(legend.xy))
         legend.xy<-c(min(xy),max(xy))
      legend(legend.xy[1],legend.xy[2],legend=legend,fill=1:nm,cex=legend.cex)
   }
}


#gscor - Structural correlation between two or more graphs.
gscor<-function(dat,dat2=NULL,g1=NULL,g2=NULL,diag=FALSE,mode="digraph",method="anneal",reps=1000,prob.init=0.9,prob.decay=0.85,freeze.time=25,full.neighborhood=TRUE,exchange.list=0){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     stop("Identical graph orders required in gscor.")
   if(!is.null(dat2)){
     dat2<-as.sociomatrix.sna(dat2)
     if(is.list(dat2))
       stop("Identical graph orders required in gscor.")
   }
   #End pre-processing
   #Collect data and various parameters
   if(is.null(g1))
     g1<-1:dim(dat)[1]
   if(is.null(g2))
     g2<-1:dim(dat)[1]
   if(!is.null(dat2)){
      if(length(dim(dat))>2)
         temp1<-dat
      else{
         temp1<-array(dim=c(1,dim(dat)[2],dim(dat)[2]))
         temp1[1,,]<-dat
      }
      if(length(dim(dat2))>2)
         temp2<-dat2
      else{
         temp2<-array(dim=c(1,dim(dat2)[2],dim(dat2)[2]))
         temp2[1,,]<-dat2
      }
      if(dim(temp1)[2]>dim(temp2)[2])
         temp2<-add.isolates(temp2,dim(temp1)[2]-dim(temp2)[2])
      if(dim(temp2)[2]>dim(temp1)[2])
         temp1<-add.isolates(temp1,dim(temp2)[2]-dim(temp1)[2])
      n<-dim(temp1)[2]
      gn<-dim(temp1)[1]+dim(temp2)[1]
      gn1<-dim(temp1)[1]
      gn2<-dim(temp2)[1]
      d<-array(dim=c(gn,n,n))
      d[1:gn1,,]<-temp1
      d[(gn1+1):(gn2+gn1),,]<-temp2
      g1<-1:gn1
      g2<-(gn1+1):(gn1+gn2)
   }else{
      d<-dat
      n<-dim(dat)[2]
      gn<-dim(dat)[1]
      gn1<-length(g1)
      gn2<-length(g2)
   }
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,gn*n),nrow=gn,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,gn)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Compute the structural correlation matrix
   gd<-matrix(nrow=gn1,ncol=gn2)
   rownames(gd)<-g1
   colnames(gd)<-g2
   if(method=="none"){
      for(i in 1:gn1)
         for(j in 1:gn2){
            d1<-d[g1[i],order(el[1,]),order(el[1,])]  #Reorder d1
            d2<-d[g2[j],order(el[2,]),order(el[2,])]  #Reorder d2
            if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
               stop("Illegal exchange list; lists must be comparable!\n")
            gd[i,j]<-cor(as.vector(d1),as.vector(d2),use="complete.obs")
         }
   }else if(method=="exhaustive"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.exhaustive(d[g1[i],,],d[g2[j],,],function(m1,m2){cor(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max")
   }else if(method=="anneal"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.anneal(d[g1[i],,],d[g2[j],,],function(m1,m2){cor(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max",prob.init=prob.init,prob.decay=prob.decay,freeze.time=freeze.time,full.neighborhood=full.neighborhood)
   }else if(method=="hillclimb"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.hillclimb(d[g1[i],,],d[g2[j],,],function(m1,m2){cor(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max")
   }else if(method=="mc"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.mc(d[g1[i],,],d[g2[j],,],function(m1,m2){cor(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max",draws=reps)
   }
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      gd[1,1]
   else
      gd
}


#gscov - Structural covariance between two or more graphs.
gscov<-function(dat,dat2=NULL,g1=NULL,g2=NULL,diag=FALSE,mode="digraph",method="anneal",reps=1000,prob.init=0.9,prob.decay=0.85,freeze.time=25,full.neighborhood=TRUE,exchange.list=0){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     stop("Identical graph orders required in gscov.")
   if(!is.null(dat2)){
     dat2<-as.sociomatrix.sna(dat2)
     if(is.list(dat2))
       stop("Identical graph orders required in gscov.")
   }
   #End pre-processing
   #Collect data and various parameters
   if(is.null(g1))
     g1<-1:dim(dat)[1]
   if(is.null(g2))
     g2<-1:dim(dat)[1]
   if(!is.null(dat2)){
      if(length(dim(dat))>2)
         temp1<-dat
      else{
         temp1<-array(dim=c(1,dim(dat)[2],dim(dat)[2]))
         temp1[1,,]<-dat
      }
      if(length(dim(dat2))>2)
         temp2<-dat2
      else{
         temp2<-array(dim=c(1,dim(dat2)[2],dim(dat2)[2]))
         temp2[1,,]<-dat2
      }
      if(dim(temp1)[2]>dim(temp2)[2])
         temp2<-add.isolates(temp2,dim(temp1)[2]-dim(temp2)[2])
      if(dim(temp2)[2]>dim(temp1)[2])
         temp1<-add.isolates(temp1,dim(temp2)[2]-dim(temp1)[2])
      n<-dim(temp1)[2]
      gn<-dim(temp1)[1]+dim(temp2)[1]
      gn1<-dim(temp1)[1]
      gn2<-dim(temp2)[1]
      d<-array(dim=c(gn,n,n))
      d[1:gn1,,]<-temp1
      d[(gn1+1):(gn2+gn1),,]<-temp2
      g1<-1:gn1
      g2<-(gn1+1):(gn1+gn2)
   }else{
      d<-dat
      n<-dim(dat)[2]
      gn<-dim(dat)[1]
      gn1<-length(g1)
      gn2<-length(g2)
   }
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,gn*n),nrow=gn,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,gn)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Compute the structural covariance matrix
   gd<-matrix(nrow=gn1,ncol=gn2)
   rownames(gd)<-g1
   colnames(gd)<-g2
   if(method=="none"){
      for(i in 1:gn1)
         for(j in 1:gn2){
            d1<-d[g1[i],order(el[1,]),order(el[1,])]  #Reorder d1
            d2<-d[g2[j],order(el[2,]),order(el[2,])]  #Reorder d2
            if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
               stop("Illegal exchange list; lists must be comparable!\n")
            gd[i,j]<-cov(as.vector(d1),as.vector(d2),use="complete.obs")
         }
   }else if(method=="exhaustive"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.exhaustive(d[g1[i],,],d[g2[j],,],function(m1,m2){cov(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max")
   }else if(method=="anneal"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.anneal(d[g1[i],,],d[g2[j],,],function(m1,m2){cov(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max",prob.init=prob.init,prob.decay=prob.decay,freeze.time=freeze.time,full.neighborhood=full.neighborhood)
   }else if(method=="hillclimb"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.hillclimb(d[g1[i],,],d[g2[j],,],function(m1,m2){cov(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max")
   }else if(method=="mc"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.mc(d[g1[i],,],d[g2[j],,],function(m1,m2){cov(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max",draws=reps)
   }
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      gd[1,1]
   else
      gd
}


#hdist - Find the Hamming distances between two or more graphs.
hdist<-function(dat,dat2=NULL,g1=NULL,g2=NULL,normalize=FALSE,diag=FALSE,mode="digraph"){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     stop("Identical graph orders required in hdist.")
   if(!is.null(dat2)){
     dat2<-as.sociomatrix.sna(dat2)
     if(is.list(dat2))
       stop("Identical graph orders required in hdist.")
   }
   #End pre-processing
   #Collect data and various parameters
   if(is.null(g1))
     g1<-1:dim(dat)[1]
   if(is.null(g2))
     g2<-1:dim(dat)[1]
   if(!is.null(dat2)){
      if(length(dim(dat))>2)
         temp1<-dat
      else{
         temp1<-array(dim=c(1,dim(dat)[2],dim(dat)[2]))
         temp1[1,,]<-dat
      }
      if(length(dim(dat2))>2)
         temp2<-dat2
      else{
         temp2<-array(dim=c(1,dim(dat2)[2],dim(dat2)[2]))
         temp2[1,,]<-dat2
      }
      if(dim(temp1)[2]>dim(temp2)[2])
         temp2<-add.isolates(temp2,dim(temp1)[2]-dim(temp2)[2])
      if(dim(temp2)[2]>dim(temp1)[2])
         temp1<-add.isolates(temp1,dim(temp2)[2]-dim(temp1)[2])
      n<-dim(temp1)[2]
      gn<-dim(temp1)[1]+dim(temp2)[1]
      gn1<-dim(temp1)[1]
      gn2<-dim(temp2)[1]
      d<-array(dim=c(gn,n,n))
      d[1:gn1,,]<-temp1
      d[(gn1+1):(gn2+gn1),,]<-temp2
      g1<-1:gn1
      g2<-(gn1+1):(gn1+gn2)
   }else{
      d<-dat
      n<-dim(dat)[2]
      gn<-dim(dat)[1]
      gn1<-length(g1)
      gn2<-length(g2)
   }
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #Compute the raw distance matrix
   hd<-matrix(nrow=gn1,ncol=gn2)
   rownames(hd)<-g1
   colnames(hd)<-g2
   for(i in 1:gn1)
      for(j in 1:gn2)
         hd[i,j]<-sum(abs(d[g1[i],,]-d[g2[j],,]),na.rm=TRUE)
   #Normalize if need be
   if(normalize)
      hd<-hd/nties(dat[1,,],mode=mode,diag=diag)
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      hd[1,1]
   else
      hd
}


#sdmat - Estimate the matrix of structural distances among a set of unlabeled 
#graphs.
#NOTE: This is redundant, as the included functionality is now included within 
#hdist and structdist, but the function is left for reasons of compatibility as 
#well as speed (currently, the distance functions don't check for duplicate 
#calculations when building distance matrices).
sdmat<-function(dat,normalize=FALSE,diag=FALSE,mode="digraph",output="matrix",method="mc",exchange.list=NULL,...){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     stop("Identical graph orders required in sdmat.")
   #End pre-processing
   if(is.null(exchange.list))
     exchange.list<-rep(0,dim(dat)[2])
   m<-dim(dat)[1]
   sdm<-matrix(nrow=m,ncol=m)
   for(i in 1:m)
      if(i<m)
         for(j in i:m)
         sdm[i,j]<-structdist(dat,g1=i,g2=j,normalize=normalize,diag=diag, mode=mode,method=method,exchange.list=exchange.list,...)
   diag(sdm)<-0
   for(i in 1:m)
      sdm[i,c(1:i)[-i]]<-sdm[c(1:i)[-i],i]
   if(output=="dist"){
      sdm<-as.dist(sdm)
   }
   sdm
}


#structdist - Estimate the structural distance between two or more unlabeled 
#graphs.
structdist<-function(dat,g1=NULL,g2=NULL,normalize=FALSE,diag=FALSE,mode="digraph",method="anneal",reps=1000,prob.init=0.9,prob.decay=0.85,freeze.time=25,full.neighborhood=TRUE,mut=0.05,pop=20,trials=5,exchange.list=NULL){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     stop("Identical graph orders required in structdist.")
   #End pre-processing
   #Collect data and various parameters
   d<-dat
   n<-dim(dat)[2]
   if(is.null(g1))
     g1<-1:dim(dat)[1]
   if(is.null(g2))
     g2<-1:dim(dat)[1]
   if(is.null(exchange.list))
     exchange.list<-rep(0,dim(dat)[2])
   gn<-dim(dat)[1]
   gn1<-length(g1)
   gn2<-length(g2)
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){             #Single number case
         el<-matrix(rep(exchange.list,gn*n),nrow=gn,ncol=n)
      }else{                                    #Vector case
         el<-sapply(exchange.list,rep,gn)
      }  
   }else                #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Compute the distance matrix
   hd<-matrix(nrow=gn1,ncol=gn2)
   rownames(hd)<-g1
   colnames(hd)<-g2
   if(method=="exhaustive"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            hd[i,j]<-lab.optimize.exhaustive(d[g1[i],,],d[g2[j],,],function(m1,m2){sum(abs(m1-m2),na.rm=TRUE)},exchange.list=el[c(g1[i],g2[j]),],seek="min")
   }else if(method=="anneal"){
      for(i in 1:gn1)
         for(j in 1:gn2){
            hd[i,j]<-lab.optimize.anneal(d[g1[i],,],d[g2[j],,],function(m1,m2){sum(abs(m1-m2),na.rm=TRUE)},exchange.list=el[c(g1[i],g2[j]),],seek="min",prob.init=prob.init,prob.decay=prob.decay,freeze.time=freeze.time,full.neighborhood=full.neighborhood)
         }
   }else if(method=="hillclimb"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            hd[i,j]<-lab.optimize.hillclimb(d[g1[i],,],d[g2[j],,],function(m1,m2){sum(abs(m1-m2),na.rm=TRUE)},exchange.list=el[c(g1[i],g2[j]),],seek="min")
   }else if(method=="mc"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            hd[i,j]<-lab.optimize.mc(d[g1[i],,],d[g2[j],,],function(m1,m2){sum(abs(m1-m2),na.rm=TRUE)},exchange.list=el[c(g1[i],g2[j]),],seek="min",draws=reps)
   }else if(method=="ga"){
      #This is broken right now - exit with a warning
      stop("Sorry, GA mode is not currently supported.\n")
   }else if(method=="none"){
      for(i in 1:gn1)
         for(j in 1:gn2){
            d1<-d[g1[i],order(el[1,]),order(el[1,])]  #Reorder d1
            d2<-d[g2[j],order(el[2,]),order(el[2,])]  #Reorder d2
            if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
               stop("Illegal exchange list; lists must be comparable!\n")
            hd[i,j]<-sum(abs(d1-d2),na.rm=TRUE)
         }
   }else{
      cat("Method",method,"not implemented yet.\n")
   }
   #Normalize if need be
   if(normalize)
      hd<-hd/nties(dat[1,,],mode=mode,diag=diag)
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      hd[1,1]
   else
      hd
}

