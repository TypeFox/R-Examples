######################################################################
#
# nli.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/27/13
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines for calculating node-level indices.
# 
# Contents:
#   betweenness
#   bonpow
#   closeness
#   degree
#   evcent
#   flowbet
#   graphcent
#   infocent
#   stresscent
#
######################################################################


#betweenness - Find the betweenness centralities of network positions
betweenness<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE,ignore.eval=TRUE){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],betweenness,g=1,nodes=nodes,gmode=gmode, diag=diag,tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp, rescale=rescale,ignore.eval=ignore.eval))
   #End pre-processing
   n<-attr(dat,"n")
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      #Note that I'm currently kludging some of these cases...could be iffy.
      if(!(cmode%in%c("directed","undirected"))){
        star<-rbind(rep(1,n-1),2:n,rep(1,n-1))
        if(gmode=="graph")
          star<-cbind(star,star[,c(2,1,3)])
        attr(star,"n")<-n
        bet<-betweenness(star,g=1,nodes=1:n,gmode=gmode,diag=diag,tmaxdev=FALSE, cmode=cmode,geodist.precomp=NULL,rescale=FALSE,ignore.eval=ignore.eval)
        bet<-sum(max(bet)-bet)
      }
      if(gmode=="graph")
        cmode<-"undirected"
      bet<-switch(cmode,
         directed = (n-1)^2*(n-2),
         undirected = (n-1)^2*(n-2)/2,
         bet
      )
   }else{
      #First, set things up
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(cmode=="undirected")   #Symmetrize if need be
        dat<-symmetrize(dat,rule="weak",return.as.edgelist=TRUE)
      meas<-switch(cmode,
        undirected=0,
        directed=0,
        endpoints=1,
        proximalsrc=2,
        proximaltar=3,
        proximalsum=4,
        lengthscaled=5,
        linearscaled=6
      )
      if(!is.null(geodist.precomp)){
        if(is.null(geodist.precomp$gdist) || is.null(geodist.precomp$counts) || is.null(geodist.precomp$predecessors)){
          warning("Precomputed geodist output must include distance, count, and predecessor information (at least one of which was missing in geodist.precomp).  Re-computing on the fly.\n")
          precomp<-FALSE
        }else
          precomp<-TRUE
      }else{
        precomp<-FALSE
      }
      #Do the computation
      bet<-.Call("betweenness_R",dat,n,NROW(dat),meas,precomp,ignore.eval, geodist.precomp$gdist,geodist.precomp$counts,geodist.precomp$predecessors,NAOK=TRUE,PACKAGE="sna")
      if((cmode=="undirected")||(gmode=="graph"))
         bet<-bet/2
      #Return the results
      if(rescale)
         bet<-bet/sum(bet)
      bet<-bet[nodes]
   }
   bet
}


#bonpow - Find the Bonacich power centrality scores of network positions
bonpow<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,exponent=1,rescale=FALSE,tol=1e-7){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],bonpow,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,exponent=exponent,rescale=rescale,tol=tol))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,bonpow,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,exponent=exponent,rescale=rescale,tol=tol))
   #End pre-processing
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      if(gmode=="graph")
         ev<-(dim(dat)[2]-2)*sqrt(dim(dat)[2]/2)
      else
         ev<-sqrt(dim(dat)[2])*(dim(dat)[2]-1)
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         d<-dat[g,,]
      else
         d<-dat
      n<-dim(d)[1]
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(!diag)
         diag(d)<-0
      #Make an identity matrix
      id<-matrix(rep(0,n*n),nrow=n)
      diag(id)<-1
      #Do the computation
      ev<-apply(solve(id-exponent*d,tol=tol)%*%d,1,sum)  #This works, when it works.
      #Apply the Bonacich scaling, by default (sum of squared ev=n)
      ev<-ev*sqrt(n/sum((ev)^2))
      if(rescale)
         ev<-ev/sum(ev)
      ev[nodes]
   }
   ev
}


#closeness - Find the closeness centralities of network positions
closeness<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE,ignore.eval=TRUE){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],closeness,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp,rescale=rescale,ignore.eval=ignore.eval))
   #End pre-processing
   n<-attr(dat,"n")
   if(gmode=="graph"){
     cmode<-switch(cmode,
       directed = "undirected",
       undirected = "undirected",
       suminvdir = "siminvundir",
       suminvundir = "suminvundir",
     )
   }
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      clo<-switch(cmode,
         directed = (n-1)*(1-1/n),    #Depends on n subst for max distance
         undirected = (n-2)*(n-1)/(2*n-3),
         suminvdir = (n-1)*(n-1),
         suminvundir = (n-2-(n-2)/2)*(n-1),
      )
   }else{
      #First, prepare the data
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(cmode%in%c("undirected","suminvundir"))   #Symmetrize if need be
        dat<-symmetrize(dat,rule="weak",return.as.edgelist=TRUE)
      #Do the computation
      if(is.null(geodist.precomp))
         gd<-geodist(dat,count.paths=FALSE,predecessors=FALSE, ignore.eval=ignore.eval)
      else
         gd<-geodist.precomp
      diag(gd$gdist)<-NA
      clo<-switch(cmode,
        directed = (n-1)/rowSums(gd$gdist,na.rm=TRUE),
        undirected = (n-1)/rowSums(gd$gdist,na.rm=TRUE),
        suminvdir = rowSums(1/gd$gdist,na.rm=TRUE)/(n-1),
        suminvundir = rowSums(1/gd$gdist,na.rm=TRUE)/(n-1)
      )
      if(rescale)
         clo<-clo/sum(clo)
      clo<-clo[nodes]
   }
   #Return the results
   clo
}


#degree - Find the degree centralities of network positions
degree<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="freeman",rescale=FALSE,ignore.eval=FALSE){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],degree,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,rescale=rescale))
   #End pre-processing
   n<-attr(dat,"n")
   if(gmode=="graph")
     cmode<-"indegree"
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      if(gmode=="digraph")
        deg<-switch(cmode,
           indegree = (n-1)*(n-1+diag),
           outdegree = (n-1)*(n-1+diag),
           freeman = (n-1)*(2*(n-1)-2+diag)
        )
      else
        deg<-switch(cmode,
           indegree = (n-1)*(n-2+diag),
           outdegree = (n-1)*(n-2+diag),
           freeman = (n-1)*(2*(n-1)-2+diag)
        )
   }else{
      #Set things up
      m<-NROW(dat)
      cm<-switch(cmode,
        indegree = 0,
        outdegree = 1,
        freeman = 2
      )
      if(!(cmode%in%c("indegree","outdegree","freeman")))
        stop("Unknown cmode in degree.\n")
      #Calculate the scores
      deg<-.C("degree_R",as.double(dat),as.integer(m),as.integer(cm), as.integer(diag),as.integer(ignore.eval),deg=as.double(rep(0,n)), PACKAGE="sna",NAOK=TRUE)$deg
      if(rescale)
         deg<-deg/sum(deg)
      if(!is.null(nodes))
        deg<-deg[nodes]
   }
   deg
}


#evcent - Find the eigenvector centralities of network positions
evcent<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,rescale=FALSE,ignore.eval=FALSE,tol=1e-10,maxiter=1e5,use.eigen=FALSE){
   #Pre-process the raw input
   if(!use.eigen){
     dat<-as.edgelist.sna(dat)
     if(is.list(dat))
       return(sapply(dat[g],evcent,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,rescale=rescale,maxiter=maxiter,use.eigen=use.eigen))
   }else{
     dat<-as.sociomatrix.sna(dat,simplify=FALSE)
     if(is.list(dat))
       return(sapply(dat[g],evcent,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,rescale=rescale,maxiter=maxiter,use.eigen=use.eigen))
   }
   #End pre-processing
   if(use.eigen)
     n<-NROW(dat)
   else
     n<-attr(dat,"n")
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      if(gmode=="graph"){
         ev<-sqrt(2)/2*(n-2)
      }else
         ev<-n-1
   }else{
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(!diag){
        if(use.eigen)
          diag(dat)<-0
        else
          dat[dat[,1]==dat[,2],3]<-0
      }
      #Do the computation
      if(use.eigen){
        ev<-eigen(dat)$vectors[,1]
      }else{
        ev<-.C("evcent_R",as.double(dat),as.integer(n),as.integer(NROW(dat)), ev=as.double(rep(1,n)),as.double(tol),as.integer(maxiter),as.integer(1),as.integer(ignore.eval),NAOK=TRUE,PACKAGE="sna")$ev
      }
      if(rescale)
         ev<-ev/sum(ev)
      ev<-ev[nodes]
   }
   ev
}


#flowbet - Find the flow betweenness of network positions
#Note: this routine is fully non-optimized
flowbet<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="rawflow",rescale=FALSE,ignore.eval=FALSE){
  #Pre-process the raw input
  dat<-as.sociomatrix.sna(dat)
  if(is.list(dat))
    return(sapply(dat[g],flowbet,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,rescale=rescale,ignore.eval=ignore.eval))
  else if((length(g)>1)&&(length(dim(dat))>2))
    return(apply(dat[g,,],1,flowbet,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,rescale=rescale,ignore.eval=ignore.eval))
  #End pre-processing
  n<-NROW(dat)
  if(ignore.eval)
    dat<-dat>0
  if(tmaxdev){
    #We got off easy: just return the theoretical maximum deviation for the centralization routine
    flo<-switch(cmode,   #This only works if we assume unit capacities!
      rawflow=(n-1)^2*(n-2)/(1+(gmode=="graph")),
      normflow=n-1,
      fracflow=(n-1)^2*(n-2)/(1+(gmode=="graph"))
    )
  }else{
    #Wrapper for the Edmonds-Karp max-flow algorithm
    mflow<-function(x,src,snk){
      .C("maxflow_EK_R",as.double(x),as.integer(NROW(x)),as.integer(src-1),
        as.integer(snk-1),flow=as.double(0),NAOK=TRUE,PACKAGE="sna")$flow
    }
    #Start by obtaining all-pairs max-solutions
    maxflo<-matrix(Inf,n,n)
    if(gmode=="digraph"){
      for(i in 1:n)
        for(j in 1:n)
          if(i!=j)
            maxflo[i,j]<-mflow(dat,i,j)
    }else{
      for(i in 1:n)
        for(j in (i:n)[-1])
          maxflo[i,j]<-mflow(dat,i,j)
      maxflo[lower.tri(maxflo)]<-t(maxflo)[lower.tri(maxflo)]
    }
    if(cmode=="normflow"){
      flo<-maxflo
      diag(flo)<-0
      maxoflo<-rep(0,n)
      for(i in 1:n)
        maxoflo[i]<-sum(flo[-i,-i])
    }
    #Compute the flow betweenness scores
    flo<-rep(0,n)
    for(i in 1:n){
      for(j in 1:n)
        for(k in 1:n)
          if((i!=j)&&(i!=k)&&(j!=k)&&((gmode=="digraph")||j<k) &&(maxflo[j,k]>0)){
            redflow<-mflow(dat[-i,-i],j-(j>i),k-(k>i))
            flo[i]<-switch(cmode,
              rawflow=flo[i]+maxflo[j,k]-redflow,
              normflow=flo[i]+maxflo[j,k]-redflow,
              fracflow=flo[i]+(maxflo[j,k]-redflow)/maxflo[j,k]
            )
          }
    }
    if(cmode=="normflow")
      flo<-flo/maxoflo*(1+(gmode=="graph"))
    if(rescale)
      flo<-flo/sum(flo)
    if(is.null(nodes))
      nodes<-1:n
    flo<-flo[nodes]
  }
  flo 
}


#graphcent - Find the graph centralities of network positions
graphcent<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE,ignore.eval=TRUE){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],graphcent,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp,rescale=rescale,ignore.eval=ignore.eval))
   #End pre-processing
   if(gmode=="graph")   #If the data is symmetric, treat it as such
      cmode<-"undirected"
    n<-attr(dat,"n")
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      gc<-switch(cmode,
         directed = (n-1)*(1-1/n),  #Depends on n subst for infinite distance
         undirected = (n-1)/2
      )
   }else{
      #First, prepare the data
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(cmode=="undirected")   #Symmetrize if need be
        dat<-symmetrize(dat,rule="weak",return.as.edgelist=TRUE)
      #Do the computation
      if(is.null(geodist.precomp))
         gd<-geodist(dat,count.paths=FALSE,predecessors=FALSE, ignore.eval=ignore.eval)
      else
         gd<-geodist.precomp
      gc<-apply(gd$gdist,1,max)
      gc<-1/gc
      if(rescale)
         gc<-gc/sum(gc)
      gc<-gc[nodes]
   }
   #Return the results
   gc
}


# infocent - Find actor information centrality scores
# Wasserman & Faust pp. 192-197; based on code generously submitted by David
# Barron (thanks!) and tweaked by myself to enable compatibility with the
# centralization() routine.
infocent <- function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,cmode="weak",tmaxdev=FALSE,rescale=FALSE,tol=1e-20){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],infocent,g=1,nodes=nodes,gmode=gmode,diag=diag, cmode=cmode,tmaxdev=tmaxdev,rescale=rescale,tol=tol))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,infocent,g=1,nodes=nodes,gmode=gmode,diag=diag, cmode=cmode,tmaxdev=tmaxdev,rescale=rescale,tol=tol))
   #End pre-processing
   if(tmaxdev){  #If necessary, return the theoretical maximum deviation
      #We don't know the real maximum value...return the lone dyad instead
      m<-matrix(0,nrow=dim(dat)[2],ncol=dim(dat)[2])
      m[1,2]<-1
      m[2,1]<-1
      IC<-infocent(m,1,rescale=rescale)  #Get ICs for dyad
      cent<-sum(max(IC)-IC,na.rm=TRUE)    #Return the theoretical max deviation 
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         m<-dat[g,,]
      else
         m<-dat
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:dim(dat)[2]
      if(sum(m != t(m),na.rm=TRUE) > 0)   #test to see if directed
         m <- symmetrize(m,rule=cmode)    #if not, we have to symmetrize...
      n <- dim(m)[1]
      if(!diag) 
         diag(m)<-NA   # if diag=F set diagonal to NA
      iso <- is.isolate(m,1:n,diag=diag) # check for isolates
      ix <- which(!iso)
      m <- m[ix,ix]           # remove any isolates (can't invert A otherwise)
      A<-1-m
      A[m==0] <- 1
      diag(A) <- 1 + apply(m, 1, sum, na.rm=TRUE)
      Cn <- solve(A,tol=tol)
      Tr <- sum(diag(Cn))
      R <- apply(Cn, 1, sum)
      IC <- 1/(diag(Cn) + (Tr - 2*R)/n)   # Actor information centrality
      #Add back the isolates
      cent<-rep(0,n)
      cent[ix]<-IC
      #Rescale if needed
      if(rescale)
         cent<-cent/sum(cent)
      #Subset as requested
      cent<-cent[nodes]
   }
   #Return the result
   cent
}


loadcent<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE,ignore.eval=TRUE){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],loadcent,g=1,nodes=nodes,gmode=gmode, diag=diag,tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp, rescale=rescale,ignore.eval=ignore.eval))
   #End pre-processing
   n<-attr(dat,"n")
   if(gmode=="graph")   #If the data is symmetric, treat it as such
      cmode<-"undirected"
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      lc<-switch(cmode,
         directed = (n-1)^2*(n-2),
         undirected = (n-1)^2*(n-2)/2
      )
   }else{
      #First, set things up
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(cmode=="undirected")   #Symmetrize if need be
        dat<-symmetrize(dat,rule="weak",return.as.edgelist=TRUE)
      else
        dat<-gt(dat,return.as.edgelist=TRUE)  #Transpose the input digraph
      if(!is.null(geodist.precomp)){
        if(is.null(geodist.precomp$gdist) || is.null(geodist.precomp$counts) || is.null(geodist.precomp$predecessors)){
          warning("Precomputed geodist output must include distance, count, and predecessor information (at least one of which was missing in geodist.precomp).  Re-computing on the fly.\n")
          precomp<-FALSE
        }else
          precomp<-TRUE
      }else{
        precomp<-FALSE
      }
      #Do the computation (we use the betweenness routine, oddly)
      lc<-.Call("betweenness_R",dat,n,NROW(dat),8,precomp,ignore.eval, geodist.precomp$gdist,geodist.precomp$counts,geodist.precomp$predecessors,NAOK=TRUE,PACKAGE="sna")
      #Return the results
      if(rescale)
         lc<-lc/sum(lc)
      lc<-lc[nodes]
   }
   lc
}


#prestige - Find actor prestige scores from one of several measures
prestige<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,cmode="indegree",tmaxdev=FALSE,rescale=FALSE,tol=1e-7){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],prestige,g=1,nodes=nodes,gmode=gmode,diag=diag, cmode=cmode,tmaxdev=tmaxdev,rescale=rescale,tol=tol))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,prestige,g=1,nodes=nodes,gmode=gmode,diag=diag, cmode=cmode,tmaxdev=tmaxdev,rescale=rescale,tol=tol))
   #End pre-processing
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      n<-dim(dat)[2]
      if(cmode=="indegree")
         p<-degree(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="indegree.rownorm")
         p<-degree(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="indegree.rowcolnorm")
         p<-degree(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="eigenvector")
         p<-evcent(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag)
      else if(cmode=="eigenvector.rownorm")
         p<-evcent(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag)
      else if(cmode=="eigenvector.colnorm")
         p<-evcent(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag)
      else if(cmode=="eigenvector.rowcolnorm")
         p<-evcent(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag)
      else if(cmode=="domain"){
         p<-(n-1)^2
      }else if(cmode=="domain.proximity"){
         p<-(n-1)^2
      }else
         stop(paste("Cmode",cmode,"unknown.\n"))      
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         d<-dat[g,,]
      else
         d<-dat
      n<-dim(d)[1]
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(!diag)
         diag(d)<-0
      #Now, perform the computation
      if(cmode=="indegree")
         p<-degree(dat=dat,g=g,nodes=nodes,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="indegree.rownorm")
         p<-degree(dat=make.stochastic(d,mode="row"),g=1,nodes=nodes,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="indegree.rowcolnorm")
         p<-degree(dat=make.stochastic(d,mode="rowcol"),g=1,nodes=nodes,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="eigenvector")
         p<-eigen(t(d))$vector[,1]      
      else if(cmode=="eigenvector.rownorm")
         p<-eigen(t(make.stochastic(d,mode="row")))$vector[,1]      
      else if(cmode=="eigenvector.colnorm")
         p<-eigen(t(make.stochastic(d,mode="col")))$vector[,1]      
      else if(cmode=="eigenvector.rowcolnorm")
         p<-eigen(t(make.stochastic(d,mode="rowcol")))$vector[,1]
      else if(cmode=="domain"){
         r<-reachability(d)
         p<-apply(r,2,sum)-1
      }else if(cmode=="domain.proximity"){
         g<-geodist(d)
         p<-(apply(g$counts>0,2,sum)-1)^2/(apply((g$counts>0)*(g$gdist),2,sum)*(n-1))
         p[is.nan(p)]<-0
      }else
         stop(paste("Cmode",cmode,"unknown.\n"))      
      if(rescale)
         p<-p/sum(p)
      p<-p[nodes]
   }  
   p
}


#stresscent - Find the stress centralities of network positions
stresscent<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE,ignore.eval=TRUE){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],stresscent,g=1,nodes=nodes,gmode=gmode, diag=diag,tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp, rescale=rescale,ignore.eval=ignore.eval))
   #End pre-processing
   n<-attr(dat,"n")
   if(gmode=="graph")   #If the data is symmetric, treat it as such
      cmode<-"undirected"
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      str<-switch(cmode,
         directed = (n-1)^2*(n-2),
         undirected = (n-1)^2*(n-2)/2
      )
   }else{
      #First, set things up
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(cmode=="undirected")   #Symmetrize if need be
        dat<-symmetrize(dat,rule="weak",return.as.edgelist=TRUE)
      if(!is.null(geodist.precomp)){
        if(is.null(geodist.precomp$gdist) || is.null(geodist.precomp$counts) || is.null(geodist.precomp$predecessors)){
          warning("Precomputed geodist output must include distance, count, and predecessor information (at least one of which was missing in geodist.precomp).  Re-computing on the fly.\n")
          precomp<-FALSE
        }else
          precomp<-TRUE
      }else{
        precomp<-FALSE
      }
      #Do the computation (we use the betweenness routine, oddly)
      str<-.Call("betweenness_R",dat,n,NROW(dat),7,precomp,ignore.eval, geodist.precomp$gdist,geodist.precomp$counts,geodist.precomp$predecessors,NAOK=TRUE,PACKAGE="sna")
      if(cmode=="undirected")
         str<-str/2
      #Return the results
      if(rescale)
         str<-str/sum(str)
      str<-str[nodes]
   }
   str
}
