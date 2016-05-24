######################################################################
#
# gli.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/27/13
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains various routines for the calculation of 
# graph-level indices.
#
# Contents:
#   centralization
#   connectedness
#   dyad.census
#   efficiency
#   gden
#   grecip
#   gtrans
#   hierarchy
#   lubness
#   mutuality
#   triad.census
#   triad.classify
#
######################################################################


#centralization - Find the centralization of a graph (for some arbitrary 
#centrality measure)
centralization<-function(dat,FUN,g=NULL,mode="digraph",diag=FALSE,normalize=TRUE,...){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(mapply(centralization,dat[g],MoreArgs=list(FUN=FUN,g=1,mode=mode, diag=diag, normalize=normalize,...)))
   }
   #End pre-processing
   #Find the centrality function
   fun<-match.fun(FUN)
   #Grab the vector of centralities
   cv<-fun(dat,g=g,gmode=mode,diag=diag,...)
   #Find the empirical maximum
   cmax<-max(cv)
   #Now, for the absolute deviations....
   cent<-sum(cmax-cv)
   #If we're normalizing, we'll need to get the theoretical max from our centrality function
   if(normalize)
      cent<-cent/fun(dat,g=g,gmode=mode,diag=diag,tmaxdev=TRUE,...)
   #Return the centralization
   cent
}


#connectedness - Find the Krackhardt connectedness of a graph or graph stack
connectedness<-function(dat,g=NULL){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],connectedness))
   }
   #End pre-processing
   g<-symmetrize(dat,rule="weak",return.as.edgelist=TRUE)
   n<-attr(g,"n")
   m<-NROW(g)
   if(n<=1)
     con<-1
   else
     con<-.C("connectedness_R",as.double(g),as.integer(n),as.integer(m), con=as.double(0.0),PACKAGE="sna",NAOK=TRUE)$con
   #Return the result
   con
}


#dyad.census - Return the Holland and Leinhardt MAN dyad census for a given 
#graph or graph stack
dyad.census<-function(dat,g=NULL){
   #Define an internal function to get the dyad census for a single mat
   intcalc<-function(m){
     n<-attr(m,"n")
     m<-m[m[,1]!=m[,2],,drop=FALSE]          #Kill loops, if any
     if(NROW(m)>0){
       dc<-.C("dyadcode_R",as.double(m),as.integer(n),as.integer(NROW(m)), dc=as.double(rep(0,NROW(m))),PACKAGE="sna",NAOK=TRUE)$dc
       mis<-is.na(m[,3])            #Count/remove missing dyads
       if(any(mis)){
         mis[dc%in%c(dc[mis])]<-TRUE
         dcm<-dc[mis]
         dmut<-sum(duplicated(dcm))
         dasym<-length(dcm)-2*dmut
         mc<-dmut+dasym
         dc<-dc[!mis]
       }else
         mc<-0
       mut<-sum(duplicated(dc))            #Find non-missing counts
       asym<-length(dc)-2*mut
       c(mut,asym,choose(n,2)-mut-asym-mc)
     }else
       c(0,0,choose(n,2))
   }
   #Organize the data
   dat<-as.edgelist.sna(dat)
   #Perform the census
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     man<-t(sapply(dat[g],intcalc))
   }else{
     man<-intcalc(dat)
   }
   if(length(man)==3)
     man<-matrix(man,nrow=1)
   colnames(man)<-c("Mut","Asym","Null")
   #Return the result
   man
}


#efficiency - Find the Krackhardt efficiency of a graph or graph stack
efficiency<-function(dat,g=NULL,diag=FALSE){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],efficiency,diag=diag))
   }
   #End pre-processing
   #Define an internal function, for convenience
   inteff<-function(g,diag){
      comsz<-component.dist(g,connected="weak")$csize
      reqedge<-sum(comsz-1)     #Get required edges
      maxv<-sum(comsz*(comsz-(!diag))-(comsz-1)) #Get maximum violations
      if(!diag)
         g<-diag.remove(g)
      edgec<-sum(g,na.rm=TRUE)  #Get count of actual edges
      1-(edgec-reqedge)/maxv
   }
   #Perform the actual calculation
   if(length(dim(dat))>2){
      if(is.null(g))
        g<-1:dim(dat)[1]
      eff<-apply(dat[g,,,drop=FALSE],1,inteff,diag=diag)
   }else
      eff<-inteff(dat,diag=diag)
   #Return the result
   eff
}


#gden - Compute the density of an input graph or graph stack.
gden<-function(dat,g=NULL,diag=FALSE,mode="digraph",ignore.eval=FALSE){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],gden,diag=diag,mode=mode,ignore.eval=ignore.eval))
   }
   #End pre-processing
   n<-attr(dat,"n")
   bip<-attr(dat,"bipartite")
   #If needed, remove loops and missing edges
   if((!diag)&&(!(mode%in%c("hgraph","twomode"))))
     dat<-dat[dat[,1]!=dat[,2],,drop=FALSE]
   nmis<-sum(is.na(dat[,3]))
   dat<-dat[!is.na(dat[,3]),,drop=FALSE]
   #Find number/value of ties, and counts
   if(n==0){
     den<-NaN
   }else if(n==1){
     if(!diag)
       den<-NaN
     else{
       if(ignore.eval)
         den<-(NROW(dat)>0)/(1-nmis)
       else
         den<-sum(dat[,3],na.rm=TRUE)/(1-nmis)
     }
   }else{
     if(ignore.eval)
       count<-NROW(dat)
     else
       count<-sum(dat[,3])
     nt<-switch(mode,
       digraph=n*(n-1)-nmis+diag*n,
       graph=n*(n-1)-nmis+diag*n,
       hgraph=bip*(n-bip)-nmis,
       twomode=bip*(n-bip)-nmis
     )
     den<-count/nt
   }
   #Return the result
   den
}


#grecip - Compute the reciprocity of an input graph or graph stack.
grecip<-function(dat,g=NULL,measure=c("dyadic","dyadic.nonnull","edgewise","edgewise.lrr","correlation")){
  #Pre-process the raw input
  dat<-as.edgelist.sna(dat)
  if(is.list(dat)){
    if(!is.null(g))
      dat<-dat[g]
  }
  #End pre-processing
  if(match.arg(measure)=="correlation"){       #Correlation measure
    if(!is.list(dat))                            #For simplicity, coerce to list
      dat<-list(dat)
    recip<-sapply(dat,function(z){             #Compute the measure
      n<-attr(z,"n")                           #Get counts
      nd<-choose(n,2)
      ne<-nd*2
      z<-z[z[,1]!=z[,2],,drop=FALSE]             #Remove loops
      #Handle a zillion special cases
      if(n<2){                                     #Only defined if n>1
        return(NA)
      }else if(n==2){                              #No var for n=2, make 0/1
        if(NROW(z)==0)
          return(1)
        else if(any(is.na(z[,3])))
          return(NA)
        else if(NROW(z)==1){
          if(z[1,3]!=0)
            return(0)
          else
            return(1)
        }else
          return((z[1,3]==z[2,3])+0)
      }
      #More general case
      if(NROW(z)>0){                             #Process edges
        emiss<-sum(is.na(z[,3]))                   #Missing edge count
        gm<-sum(z[,3],na.rm=TRUE)/(ne-emiss)       #Get graph mean and var
        gv<-(sum((z[,3]-gm)^2,na.rm=TRUE)+(ne-NROW(z))*gm^2)/(ne-emiss-1)
        if(gv==0)                                  #If var 0, treat as corr 1
          return(1)
        dc<-.C("dyadcode_R",as.double(z),as.integer(n),as.integer(NROW(z)), dc=as.double(rep(0,NROW(z))),PACKAGE="sna",NAOK=TRUE)$dc
        odc<-order(dc)
        zv<-z[odc,3]-gm  #Order by dyad ID, subtract graph mean
        dc<-dc[odc]
        dsum<-0          #Compute the product-moment correlation
        dmiss<-0
        dcount<-0
        i<-1
        while(i<=length(zv)){
          if((i<length(zv))&&(dc[i]==dc[i+1])){    #Two obs per dyad
            if(is.na(zv[i])||is.na(zv[i+1]))
              dmiss<-dmiss+1
            else{
              dsum<-dsum+zv[i]*zv[i+1]
              dcount<-dcount+1
            }
            i<-i+2
          }else{                 #One obs per dyad (other is 0, by defn)
            if(is.na(zv[i]))
              dmiss<-dmiss+1
            else{
              dsum<-dsum-gm*zv[i]
              dcount<-dcount+1
            }
            i<-i+1
          }
        }
        return(2*(dsum+gm^2*(nd-dcount-dmiss))/((2*nd-2*dmiss-1)*gv))
      }else{
        return(1)                                #All zeros - treat as corr 1
      }
    })
  }else{                                       #All other measures
    dc<-dyad.census(dat)   #Obtain the dyad census for all specified graphs
    #Compute the appropriate measure
    recip<-switch(match.arg(measure),
      dyadic=(dc[,1]+dc[,3])/(dc[,1]+dc[,2]+dc[,3]),
      dyadic.nonnull=dc[,1]/(dc[,1]+dc[,2]),
      edgewise=2*dc[,1]/(2*dc[,1]+dc[,2]),
      edgewise.lrr=log(dc[,1]*(dc[,1]+dc[,2]+dc[,3])/(dc[,1]+dc[,2]/2)^2)
    )
  }
  #Return the result
  recip
}


#gtrans - Compute the transitivity of an input graph or graph stack.
gtrans<-function(dat,g=NULL,diag=FALSE,mode="digraph",measure=c("weak","strong","weakcensus","strongcensus","rank","correlation"),use.adjacency=TRUE){
  #Perform some very crude triage - if we know the data format, and if
  #adjacency processing is obviously a bad idea, default back to edgelists
  if(use.adjacency&&(!(match.arg(measure)%in%c("correlation","rank")))){
    adjisok<-function(z){       #Is it OK to use adjacency (as far as we know?)
      if(is.edgelist.sna(z)){
        if(attr(z,"n")>40000)
          FALSE
        else if((attr(z,"n")>1000)&&(NROW(z)/attr(z,"n")^2<0.5))
          FALSE
        else
          TRUE
      }else if(class(z)=="matrix"){
        if(NCOL(z)>1000)
          FALSE
        else
          TRUE
      }else if(class(z)=="array"){
        if(dim(z)[2]>1000)
          FALSE
        else
          TRUE
      }else if(class(z)=="network"){
        require(network)
        if(network.size(z)>40000)
          FALSE
        else if((network.size(z)>1000)&& (network.edgecount(z)/network.size(z)^2<0.5))
          FALSE
        else
          TRUE
      }else
        TRUE
    }
    if(class(dat)=="list")
      adjcheck<-sapply(dat,adjisok)
    else
      adjcheck<-adjisok(dat)
    if(any(!adjcheck)){
      use.adjacency<-FALSE
      warning("gtrans called with use.adjacency=TRUE, but your data looks too large for that to work well.  Overriding to edgelist method.")
    }
  }
  if(use.adjacency&&(match.arg(measure)=="rank")){ #Only edgelist for rank
    use.adjacency<-FALSE
  }
  #End crude triage
  if((!use.adjacency)&&(match.arg(measure)=="correlation")){
    warning("Currently, non-adjacency computation for the correlation measure is not supported.  Defaulting to use.adjacency==TRUE in gtrans.\n")
    use.adjacency<-TRUE
  }
  if(use.adjacency){  #Use adjacency matrix - much faster for n<1000 or dense
    #Pre-process the raw input
    dat<-as.sociomatrix.sna(dat)
    if(is.list(dat)){
      if(is.null(g))
        g<-1:length(dat)
      return(sapply(dat[g],gtrans,diag=diag,mode=mode,measure=measure, use.adjacency=use.adjacency))
    }
    #End pre-processing
    n<-dim(dat)[2]
    if(length(dim(dat))>2){     #Is this a stack?
       if(!is.null(g)){                 #Were individual graphs selected?
          gn<-length(g)
          d<-dat[g,,]
       }else{
          d<-dat
          gn<-dim(dat)[1]
       }
    }else{
       d<-dat
       gn<-1
    }
    if(gn==1){     #Only one graph - convert to stack format
       temp<-array(dim=c(1,n,n))
       temp[1,,]<-d
       d<-temp
    }
    if(!diag)           #If not using the diagonal, remove it
       d<-diag.remove(d,remove.val=0)
    #Compute the appropriate transitivity indices
    t<-vector()
    for(i in 1:gn){
       #Prepare the transitivity test matrices
       if(match.arg(measure)!="correlation"){
         dt<-d[i,,]!=0
       }else{
         dt<-d[i,,]
       }
       dsqt<-(dt%*%dt)
       #NA the diagonal, if needed
       if(!diag){
          diag(dt)<-NA
          diag(dsqt)<-NA
        }
       #Compute the transitivity
       t[i]<-switch(match.arg(measure),
          strong=sum(dt*dsqt+(!dt)*(NCOL(d[i,,])-2-dsqt),na.rm=TRUE) / (choose(NCOL(d[i,,]),3)*6),
          strongcensus=sum(dt*dsqt+(!dt)*(NCOL(d[i,,])-2-dsqt),na.rm=TRUE),
          weak=sum(dt*dsqt,na.rm=TRUE)/sum(dsqt,na.rm=TRUE),
          weakcensus=sum(dt*dsqt,na.rm=TRUE),
          correlation=(function(x,y){
            tv<-var(x,use="pairwise.complete.obs")* var(y,use="pairwise.complete.obs")
            if(is.na(tv))
              NA
            else if(tv==0)
              all(x==y,na.rm=TRUE)+0
            else
              cor(x,y,use="pairwise.complete.obs")
          })(x=as.vector(dt),y=as.vector(dsqt))
       )
       if(is.nan(t[i]))  #By convention, map undefined case to 1
         t[i]<-1
    }
    #Return the result
    t
  }else{  #Use edgelist - much faster for large, sparse graphs
    #Pre-process the raw input
    dat<-as.edgelist.sna(dat)
    if(is.list(dat)){
      if(is.null(g))
        g<-1:length(dat)
      return(sapply(dat[g],gtrans,diag=diag,mode=mode,measure=measure, use.adjacency=use.adjacency))
    }
    #End pre-processing
    if(attr(dat,"n")<3)          #Consider vacuously transitive if n<3
      return(1)
    meas<-switch(match.arg(measure),
      "strong"=0,
      "strongcensus"=0,
      "weak"=1,
      "weakcensus"=1,
      "rank"=2,
      "correlation"=3
    )
    gt<-.C("transitivity_R",as.double(dat),as.integer(attr(dat,"n")), as.integer(NROW(dat)),gt=as.double(c(0,0)),as.integer(meas),as.integer(1),NAOK=TRUE,PACKAGE="sna")$gt
    if(match.arg(measure)%in%c("weak","strong","rank")){
      if(gt[2]==0)              #By convention, return 1 if no preconditions
        1
      else
        gt[1]/gt[2]
    }else
      gt[1]
  }
}


#hierarchy - Find the hierarchy score of a graph or graph stack
hierarchy<-function(dat,g=NULL,measure=c("reciprocity","krackhardt")){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],hierarchy,measure=measure))
   }
   #End pre-processing
   if(is.null(g))
     g<-1:stackcount(dat)
   if(match.arg(measure)=="reciprocity")  #Use reciprocity scores
         h<-1-grecip(dat,g)
   else if(match.arg(measure)=="krackhardt"){ #Calculate the Krackhardt reciprocity
      d<-array(dim=c(length(g),dim(dat)[2],dim(dat)[2]))
      if(length(dim(dat))>2)
         d<-dat[g,,,drop=FALSE]
      else
         d[1,,]<-dat
      h<-1-apply(d,1,function(x){r<-reachability(x); grecip(r,measure="dyadic.nonnull")})
   }
   #Return the result
   h
}


#lubness - Find Krackhardt's Least Upper Boundedness of a graph or graph stack
lubness<-function(dat,g=NULL){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],lubness))
   }
   #End pre-processing
   #Define an internal function, for convenience
   intlub<-function(g){
      r<-reachability(g)    #Get reachability (in directed paths) of g
      cd<-component.dist(g,connected="weak")  #Get weak components of g
      nolub<-0
      maxnolub<-0
      for(i in 1:max(cd$membership)){   #Walk through the components
         vi<-(1:dim(g)[1])[cd$membership==i]  #Get the vertices of component i
         if(length(vi)>2){  #Components must be of size 3
           #Accumulate violations
           viol<-as.double(0)
           viol<-.C("lubness_con_R",as.double(g[vi,vi]), as.double(length(vi)),as.integer(r[vi,vi]),viol=viol,PACKAGE="sna")$viol
           nolub<-nolub+viol
           #Also accumulate maximum violations
           maxnolub<-maxnolub+(length(vi)-1)*(length(vi)-2)/2 
         }
      }
      #Return 1-violations/max(violations)
      1-nolub/maxnolub
   }
   #Perform the actual calculation
   if(length(dim(dat))>2){
      if(!is.null(g))
        dat<-dat[g,,,drop=FALSE]
      lub<-apply(dat,1,intlub)
   }
   else
      lub<-intlub(dat)
   #Return the result
   lub
}


#mutuality - Find the number of mutual (i.e., reciprocated) edges in a graph
mutuality<-function(dat,g=NULL){
  #Pre-process the raw input
  dat<-as.edgelist.sna(dat)
  if(is.list(dat)){
    if(!is.null(g))
      dat<-dat[g]
  }
  #End pre-processing
  dc<-dyad.census(dat)   #Obtain the dyad census for all specified graphs
  dc[,1]                 #Return the mutual count
}


#triad.census - Conduct a Davis and Leinhardt triad census for a graph or graph stack
triad.census<-function(dat,g=NULL,mode=c("digraph","graph")){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   #End pre-processing
   #First, define the triad class vector
   tc<-switch(match.arg(mode),
     graph=0:3,
     digraph=c("003","012","102","021D","021U","021C","111D","111U","030T", "030C","201","120D","120U","120C","210","300")
   )
   #Obtain triad census scores
   if(!is.list(dat))
     dat<-list(dat)
   rnam<-names(dat)
   gm<-as.integer(switch(match.arg(mode),graph=0,digraph=1))
   tcm<-matrix(nrow=length(dat),ncol=length(tc))
   for(i in 1:length(dat)){
     n<-as.integer(attr(dat[[i]],"n"))
     m<-as.integer(NROW(dat[[i]]))
     tcv<-as.double(rep(0,length(tc)))
     if(n>2)
       tcm[i,]<-.C("triad_census_R",as.double(dat[[i]]),n,m,tcv=tcv,gm, as.integer(1),PACKAGE="sna", NAOK=TRUE)$tcv
   }
   colnames(tcm)<-tc
   rownames(tcm)<-rnam
   #Return the result
   tcm
}


#triad.classify - Return the Davis and Leinhardt classification of a given triad
triad.classify<-function(dat,g=1,tri=c(1,2,3),mode=c("digraph","graph")){
   #Zeroth step: extract the triad
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
      d<-dat[[g]][tri,tri]
   else if(length(dim(dat))==2)
      d<-dat[tri,tri]
   else
      d<-dat[g,tri,tri]
   #First, classify as NA if any entries are missing
   if(any(is.na(d[upper.tri(d)|lower.tri(d)])))
      return(NA)
   #Next, define the triad class vector
   tc<-switch(match.arg(mode),
     graph=0:3,
     digraph=c("003","012","102","021D","021U","021C","111D","111U","030T", "030C","201","120D","120U","120C","210","300")
   )
   #Classify the triad
   tt<-as.integer(0)
   gm<-as.integer(switch(match.arg(mode),graph=0,digraph=1))
   tt<-.C("triad_classify_R",as.integer(d),tt=tt,gm,PACKAGE="sna")$tt
   tc[tt+1]
}
