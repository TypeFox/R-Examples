######################################################################
#
# misc.R
#
# Written by Carter T. Butts <buttsc@uci.edu>; portions contributed by
# David Hunter <dhunter@stat.psu.edu> and Mark S. Handcock
# <handcock@u.washington.edu>.
#
# Last Modified 02/26/13
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/network package
#
# This file contains various network routines which don't fit anywhere
# else (generally, utilities and the like).
#
# Contents:
#
#   as.color
#   mixingmatrix
#   network.density
#   is.color
#   is.discrete
#   is.discrete.character
#   is.discrete.numeric
#   print.mixingmatrix
#   which.matrix.type
#
######################################################################


#Given a vector of non-colors, try to coerce them into some reasonable
#color format.  This may not work well, but what the hell....
as.color<-function(x,opacity=1.0){
  if(opacity > 1 | opacity < 0){
    stop('opacity parameter must be a numeric value in the range 0 to 1')
  }
  colors<-x
  #Numeric rule: if integer leave as-is, otherwise convert to grayscale
  if(is.numeric(x)){
    if(any(x!=round(x),na.rm=TRUE)){
      colors<-gray((x-min(x))/(max(x)-min(x)))
    }else
      colors<-x
  }
  #Factor rule: categorical colorings
  if(is.factor(x)){
    colors<-match(levels(x)[x],levels(x))
  }
  #Character rule: if colors, retain as colors; else categorical
  if(is.character(x)){
    if(all(is.color(x)))
      colors<-x
    else{
      colors<-match(x,sort(unique(x)))
    }
  }
  # add transparency if not 1
  if(opacity < 1){
    colors<-grDevices::adjustcolor(colors,alpha.f=opacity)
  }
  return(colors)
}


#Return the mixing matrix for a network object, on a given attribute.  This
#is a relocated function from the ergm package; it probably belongs elsewhere,
#but is needed for the summary.network method (and in that sense is basic
#enough to include.
mixingmatrix <- function(nw, attrname) {
  if(!is.network(nw)){
    stop("mixingmatrix() requires a network object")
  }
  if(missing(attrname)){
    stop("attrname argument is missing. mixingmatrix() requires an an attribute name")
  }
  if(network.size(nw)==0){
    warning("mixing matrices not well-defined for graphs with no vertices.")
    type<-"directed"
    if(is.bipartite(nw))
      type<-"bipartite"
    tabu<-matrix(nrow=0,ncol=0)
    ans<-list(type=type,matrix=tabu)
    class(ans)<-"mixingmatrix"
    return(ans)
  }
  nodecov <- unlist(get.vertex.attribute(nw, attrname))
  u<-sort(unique(nodecov))
  # nodecovnum <- match(nodecov, u)
  el <- as.matrix.network.edgelist(nw)
  type <- "directed"
  if (is.bipartite(nw)) { # must have heads < tails now
    if (is.directed(nw)) 
      cat("Warning:  Bipartite networks are currently\n",
          "automatically treated as undirected\n")
    type <- "bipartite"
    rowswitch <- apply(el, 1, function(x) x[1]>x[2])
    el[rowswitch, 1:2] <- el[rowswitch, 2:1]
    nb1 <- get.network.attribute(nw,"bipartite")
    u<-sort(unique(nodecov[0:nb1]))
    From <- c(u, nodecov[el[,1]])
    u<-sort(unique(nodecov[(nb1+1):network.size(nw)]))
    To <- c(u, nodecov[el[,2]])
  }else{
    From <- c(u, nodecov[el[,1]])
    To <- c(u, nodecov[el[,2]])
  }
  tabu <- table(From, To)  # Add u,u diagonal to ensure each 
  # value is represented, then subtract it later
  diag(tabu) <- diag(tabu) - 1
  if(!is.directed(nw) && !is.bipartite(nw)){
    type <- "undirected"
    tabu <- tabu + t(tabu)
    diag(tabu) <- diag(tabu)/2
  }
  ans <- list(type=type, matrix=tabu)
  class(ans) <- "mixingmatrix"
  ans
}


# Return the density of the given network.  (This probably won't stay in
# this package....
#
network.density<-function(x,na.omit=TRUE,discount.bipartite=FALSE){
  if(!is.network(x))
    stop("network.density requires a network object.")
  if(network.size(x)==0){
    warning("Density is not well-defined for networks of order 0.")
    return(NaN)
  }
  if(is.multiplex(x))
    warning("Network is multiplex - no general way to define density.  Returning value for a non-multiplex network (hope that's what you wanted).\n")
  ec<-network.edgecount(x,na.omit=na.omit)
  n<-network.size(x)
  bip<-x%n%"bipartite"
  if(is.hyper(x)){
    if((bip>=0)&&(discount.bipartite)){
      pe<-choose(bip,1:bip)*choose(n-bip,1:(n-bip))*(1+is.directed(x))
    }else{
      if(has.loops(x))
        pe<-sum(choose(n,1:n))^(1+is.directed(x))
      else
        pe<-sum(choose(n,1:n))/(1+!is.directed(x))
    }
  }else{
    if((bip>=0)&&(discount.bipartite)){
      pe<-bip*(n-bip)*(1+is.directed(x))
    }else{
      pe<-n*(n-1)/(1+!is.directed(x))+(has.loops(x)*network.size(x))
    }
  }  
  ec/pe
}

# has.edges  checks if any of the specified vertex ids have edges (are not isolates)
has.edges<-function(net,v=seq_len(network.size(net))){
  if(network.size(net)==0){
    return(logical(0))
  }
  if(any(v < 1) | any(v > network.size(net))){
    stop("'v' argument must be a valid vertex id in is.isolate")
  }
  ins<-sapply(net$iel[v],length)
  outs<-sapply(net$oel[v],length)
  return(ins+outs != 0)
}


#Returns TRUE if x is a character in a known color format
is.color<-function(x){
  xic<-rep(FALSE,length(x))         #Assume not a color by default
  
  xc<-sapply(x,is.character)        #Must be a character string
  #For characters, must be a named color or a #RRGGBB/#RRGGBBAA sequence
  xic[xc]<-(x[xc]%in%colors())| ((nchar(x[xc])%in%c(7,9))&(substr(x[xc],1,1)=="#"))
  xic[is.na(x)]<-NA                 #Missing counts as missing
  #Return the result
  xic
}


is.discrete.numeric<-function(x){
 (is.numeric(x)|is.logical(x)) && mean(duplicated(x)) > 0.8
}


is.discrete.character<-function(x){
 (is.character(x)|is.logical(x)) && mean(duplicated(x)) > 0.8
}


is.discrete<-function(x){
 (is.numeric(x)|is.logical(x)|is.character(x)) && mean(duplicated(x)) > 0.8
}


#Print method for mixingmatrix objects
print.mixingmatrix <- function(x, ...) {
  m <- x$mat
  rn <- rownames(m)
  cn <- colnames(m)  
  if (x$type == "undirected") {
    dimnames(m) <- list(rn, cn)
    cat("Note:  Marginal totals can be misleading\n",
        "for undirected mixing matrices.\n")
  } else {
    total <- apply(m,1,sum)
    m <- cbind(m,total)
    total <- apply(m,2,sum)
    m <- rbind(m,total)
    rn <- c(rn, "Total")
    cn <- c(cn, "Total")
    if (x$type == "bipartite")
      dimnames(m) <- list(B1 = rn,B2 = cn)
    else
      dimnames(m) <- list(From = rn,To = cn)
  }
  print(m)
}


which.matrix.type<-function(x)
{
    if (!is.network(x)) {
        if (is.character(x<-as.matrix(x))){ 
          if (diff(dim(x))==0)
            out<-"adjacency"
          else if (dim(x)[2]==2)
            out<-"edgelist"
          else
            out<-"bipartite"
        }else if (!is.numeric(x))  
            out<-NA
        else if (diff(dim(x))==0)  
            out<-"adjacency"
        else if (max(abs(x),na.rm=TRUE)==1 && max(abs(x-as.integer(x)),na.rm=TRUE)==0)
            out<-"bipartite"
        else if (max(abs(x-as.integer(x))[,1:2],na.rm=TRUE)==0 && min(x[,1:2],na.rm=TRUE)>0)
            out<-"edgelist"
        else
            out<-NA
    }
    else {  # Very ad-hoc criteria for choosing; choice can be overridden.
	if (is.hyper(x))
            out<-"incidence"
        else if ((n<-x$gal$n)<14 || x$gal$mnext>n*n/2)
            out<-"adjacency"
        else
            out<-"edgelist"
    }
    out
}
