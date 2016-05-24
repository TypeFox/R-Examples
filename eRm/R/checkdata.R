# uses
#   component.dist
#   reachability
#   geodist
#   symmetrize
#   components.c
#   geodist.c
#
# from R package sna


# function to check for ill-conditioned data in the RM
#   requires package sna

##checkdata<-function(x)
##{
##    k<-ncol(x)
##    adj<-matrix(0,nc=k,nr=k)
##    for (i in 1:k) for(j in 1:k) {
##        adj[i,j]<- 1*any(x[,i]>x[,j],na.rm=TRUE)
##    }
##
##    #library(sna)
##    #adj <- diag.remove(adj)
##    # %print(adj)  # adjacency marix
##    cd <- component.dist(adj, connected = "strong")
##    cm <- cd$membership
##    cmp <- max(cm)
##
##
##    if(cmp>1) {
##         cat("Data:",deparse(substitute(x)),"are ill-conditioned\n")
##         cat("Number of strong components",cmp,"\n")
##         cat("Component membership of items: ",cm,"\n")
##    } else
##         cat("Data:",deparse(substitute(x)),"are well-conditioned\n")
##}
##
######################################################
component.dist<-
function (dat, connected = c("strong", "weak", "unilateral",
    "recursive"))
{
#   dat <- as.sociomatrix.sna(dat)
#   if (is.list(dat))
#       return(lapply(dat, component.dist, connected = connected))
#   else if (length(dim(dat)) > 2)
#       return(apply(dat, 1, component.dist, connected = connected))
    n <- dim(dat)[2]
    if (any(dat != t(dat)))
        dat <- switch(match.arg(connected), weak = symmetrize(dat,
            rule = "weak"), unilateral = reachability(dat), strong = symmetrize(reachability(dat),
            rule = "strong"), recursive = symmetrize(dat, rule = "strong"))
#   if (match.arg(connected) == "unilateral")
#       if (any(dat != t(dat)))
#           warning("Nonunique unilateral component partition detected in component.dist.  Problem vertices will be arbitrarily assigned to one of their components.\n")
    membership <- rep(0, n)
    membership <- .C("component_dist_R", as.double(dat), as.double(n),
        membership = as.double(membership), PACKAGE="eRm")$membership
    o <- list()
    o$membership <- membership
    o$csize <- vector()
    for (i in 1:max(membership)) o$csize[i] <- length(membership[membership ==
        i])
    o$cdist <- vector()
    for (i in 1:n) o$cdist[i] <- length(o$csize[o$csize == i])
    o
}

#reachability - Find the reachability matrix of a graph.
reachability<-function(dat,geodist.precomp=NULL){
   #Pre-process the raw input
#   dat<-as.sociomatrix.sna(dat)
#   if(is.list(dat))
#     return(lapply(dat,reachability,geodist.precomp=geodist.precomp))
#   else if(length(dim(dat))>2)
#     return(apply(dat,1,reachability,geodist.precomp=geodist.precomp))
#     return(unlist(apply(dat,1,function(x,geodist.precomp){list(reachability(x, geodist.precomp=geodist.precomp))},geodist.precomp=geodist.precomp),recursive=FALSE))
   #End pre-processing
   #Get the counts matrix
   if(is.null(geodist.precomp))
      cnt<-geodist(dat)$counts
   else
      cnt<-geodist.precomp$counts
   #Dichotomize and return
   apply(cnt>0,c(1,2),as.numeric)
}

#geodist - Find the numbers and lengths of geodesics among nodes in a graph
#using a BFS, a la Brandes (2000).  (Thanks, Ulrik!)
geodist<-function(dat,inf.replace=Inf){
   #Pre-process the raw input
#   dat<-as.sociomatrix.sna(dat)
#   if(is.list(dat))
#     return(lapply(dat,geodist,inf.replace=inf.replace))
#   else if(length(dim(dat))>2)
#     return(apply(dat,1,geodist,inf.replace=inf.replace))
   #End pre-processing
   n<-dim(dat)[2]
   #Initialize the matrices
   sigma<-matrix(0,nrow=n,ncol=n)
   gd<-matrix(Inf,nrow=n,ncol=n)
   #Perform the calculation
   geo<-.C("geodist_R",as.double(dat),as.double(n),gd=as.double(gd), sigma=as.double(sigma),NAOK=TRUE,PACKAGE="eRm")
   #Return the results
   o<-list()
   o$counts<-matrix(geo$sigma,n,n)
   o$gdist<-matrix(geo$gd,n,n)
   o$gdist[o$gdist==Inf]<-inf.replace  #Patch Infs, if desired
   o
}

#symmetrize - Convert a graph or graph stack to a symmetric form.  Current rules
#for symmetrizing include "upper" and "lower" diagonals, "weak" connectedness
#rule, and a "strong" connectedness rule.
symmetrize<-function(mats,rule="weak"){
   #Pre-process the raw input
#   mats<-as.sociomatrix.sna(mats)
#   if(is.list(mats))
#     return(lapply(mats,symmetrize,rule=rule))
   #End pre-processing
   #Build the input data structures
#   if(length(dim(mats))>2){
#      m<-dim(mats)[1]
#      n<-dim(mats)[2]
#      o<-dim(mats)[3]
#      d<-mats
#   }else{
      m<-1
      n<-dim(mats)[1]
      o<-dim(mats)[2]
      d<-array(dim=c(1,n,o))
      d[1,,]<-mats
#   }
   #Apply the symmetry rule
   for(i in 1:m){
      if(rule=="upper"){
#         temp<-d[i,,]
#         for(j in 1:n)
#            temp[j:n,j]<-temp[j,j:n]
#         d[i,,]<-temp
#      }else if(rule=="lower"){
#         temp<-d[i,,]
#         for(j in 1:n)
#            temp[j,j:n]<-temp[j:n,j]
#         d[i,,]<-temp
#      }else if(rule=="weak"){
#         d[i,,]<-matrix(as.numeric(d[i,,]|t(d[i,,])),nrow=n,ncol=o)
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
}
