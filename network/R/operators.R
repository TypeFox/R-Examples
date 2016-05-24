######################################################################
#
# operators.R
#
# Written by Carter T. Butts <buttsc@uci.edu>; portions contributed by
# David Hunter <dhunter@stat.psu.edu> and Mark S. Handcock
# <handcock@u.washington.edu>.
#
# Last Modified 01/28/11
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/network package
#
# This file contains various operators which take networks as inputs.
#
# Contents:
#
# "$<-.network"
# "[.network"
# "[<-.network"
# "%e%"
# "%e%<-"
# "%eattr%"
# "%eattr%<-"
# "%n%"
# "%n%<-"
# "%nattr%"
# "%nattr%<-"
# "%s%"
# "%v%"
# "%v%<-"
# "%vattr%"
# "%vattr%<-"
# "+"
# "+.default"
# "+.network"
# "-"
# "-.default"
# "-.network"
# "*"
# "*.default"
# "*.network"
# "!.network"
# "|.network"
# "&.network"
# "%*%.network"
# "%c%"
# "%c%.network"
# networkOperatorSetup
# prod.network
# sum.network
# 
######################################################################

# removed this function because it appears that '<-' is no longer a generic in R, so it was never getting called and the copy was not being made. See ticket #550
"<-.network"<-function(x,value){
  .Deprecated("network.copy or '<-' works just fine",msg="The network assignment S3 method '<-.network' has been deprecated because the operator '<-' is no longer an S3 generic in R so the .network version does not appear to be called. If you see this warning, please contact the maintainers to let us know you use this function")
  x<-network.copy(value)
  return(x)
}



# removed so that will dispatch to internal primitive method #642
#"$<-.network"<-function(x,i,value){
#  cl<-oldClass(x)
#  class(x)<-NULL
#  x[[i]]<-value
#  class(x)<-cl
#  return(x)
#}


"[.network"<-function(x,i,j,na.omit=FALSE){
  narg<-nargs()+missing(na.omit)
  n<-network.size(x)
  xnames <- network.vertex.names(x)
  if(missing(i))              #If missing, use 1:n
    i<-1:n
  if((narg>3)&&missing(j))
    j<-1:n
  if(is.matrix(i)&&(NCOL(i)==1))  #Vectorize if degenerate matrix
    i<-as.vector(i)
  if(is.matrix(i)){    #Still a matrix?
    if(is.logical(i)){                    #Subset w/T/F?
      j<-col(i)[i]
      i<-row(i)[i]
      out<-is.adjacent(x,i,j,na.omit=na.omit)
    }else{                                #Were we passed a pair list?
      if(is.character(i))
        i<-apply(i,c(1,2),match,xnames)
      out<-is.adjacent(x,i[,1],i[,2], na.omit=na.omit)
    }
  }else if((narg<3)&&missing(j)){   #Here, assume a list of cell numbers
    ir<-1+((i-1)%%n)
    ic<-1+((i-1)%/%n)
    out<-is.adjacent(x,ir,ic,na.omit=na.omit)
  }else{                      #Otherwise, assume a vector or submatrix
    if(is.character(i))
      i<-match(i,xnames)
    if(is.character(j))
      j<-match(j,xnames)
    i<-(1:n)[i]                 #Piggyback on R's internal tricks
    j<-(1:n)[j]
    if(length(i)==1){
      out<-is.adjacent(x,i,j,na.omit=na.omit)
    }else{
      if(length(j)==1){
        out<-is.adjacent(x,i,j,na.omit=na.omit)
      }else{
        jrep<-rep(j,rep.int(length(i),length(j)))
        if(length(i)>0)
          irep<-rep(i,times=ceiling(length(jrep)/length(i)))
        out<-matrix(is.adjacent(x,irep,jrep,na.omit=na.omit), length(i),length(j))
      }
    }
    if((!is.null(xnames))&&is.matrix(out))
      dimnames(out) <- list(xnames[i],xnames[j])
  }
  out+0                       #Coerce to numeric
}


"[<-.network"<-function(x,i,j,names.eval=NULL,add.edges=FALSE,value){
  #Check for hypergraphicity
  if(is.hyper(x))
    stop("Assignment operator overloading does not currently support hypergraphic networks.");
  #Set up the edge list to change
  narg<-nargs()+missing(names.eval)+missing(add.edges)
  n<-network.size(x)
  xnames <- network.vertex.names(x)
  if(missing(i))              #If missing, use 1:n
    i<-1:n
  if((narg>5)&&missing(j))
    j<-1:n
  if(is.matrix(i)&&(NCOL(i)==1))  #Vectorize if degenerate matrix
    i<-as.vector(i)
  if(is.matrix(i)){    #Still a matrix?
    if(is.logical(i)){                    #Subset w/T/F?
      j<-col(i)[i]
      i<-row(i)[i]
      el<-cbind(i,j)
    }else{                                #Were we passed a pair list?
      if(is.character(i))
        i<-apply(i,c(1,2),match,xnames)
      el<-i
    }
  }else if((narg<6)&&missing(j)){  #Here, assume a list of cell numbers
    el<-1+cbind((i-1)%%n,(i-1)%/%n)
  }else{                      #Otherwise, assume a vector or submatrix
    if(is.character(i))
      i<-match(i,xnames)
    if(is.character(j))
      j<-match(j,xnames)
    i<-(1:n)[i]                 #Piggyback on R's internal tricks
    j<-(1:n)[j]
    if(length(i)==1){
      el<-cbind(rep(i,length(j)),j)
    }else{
      if(length(j)==1)
        el<-cbind(i,rep(j,length(i)))
      else{
        jrep<-rep(j,rep.int(length(i),length(j)))
        if(length(i)>0)
          irep<-rep(i,times=ceiling(length(jrep)/length(i)))
        el<-cbind(irep,jrep)
      }
    }
  }
  #Set up values
  if(is.matrix(value))
    val<-value[cbind(match(el[,1],sort(unique(el[,1]))), match(el[,2],sort(unique(el[,2]))))]
  else
    val<-rep(as.vector(value),length=NROW(el))
  #Perform the changes
  if(is.null(names.eval)){  #If no names given, don't store values
    for (k in seq_along(val)) {
      eid <- get.edgeIDs(x,el[k,1],el[k,2],neighborhood="out", na.omit=FALSE)
      if (!is.na(val[k]) & val[k] == 0) {
        # delete edge
        if (length(eid) > 0) x<-delete.edges(x,eid)
      } else {
        if (length(eid) == 0 & (has.loops(x)|(el[k,1]!=el[k,2]))) {
          # add edge if needed
          x<-add.edges(x,as.list(el[k,1]),as.list(el[k,2]))
          eid <- get.edgeIDs(x,el[k,1],el[k,2],neighborhood="out", na.omit=FALSE)
        }
        if (is.na(val[k])) {
          set.edge.attribute(x,"na",TRUE,eid)   # set to NA
        } else if (val[k] == 1) {
          set.edge.attribute(x,"na",FALSE,eid)   # set to 1
        }
      }
    }
  }else{                   #An attribute name was given, so store values
    epresent<-vector()
    eid<-vector()
    valsl<-list()
    for(k in 1:NROW(el)){
      if(is.adjacent(x,el[k,1],el[k,2],na.omit=FALSE)){  #Collect extant edges
        loceid<-get.edgeIDs(x,el[k,1],el[k,2],neighborhood="out",na.omit=FALSE)
        if(add.edges){    #Need to know if we're adding/removing edges
          if(val[k]==0){     #If 0 and adding/removing, eliminate present edges
            x<-delete.edges(x,loceid)
          }else{             #Otherwise, add as normal
            valsl<-c(valsl,as.list(rep(val[k],length(loceid))))
            eid<-c(eid,loceid)
          }
        }else{
          valsl<-c(valsl,as.list(rep(val[k],length(loceid))))
          eid<-c(eid,loceid)
        }
        epresent[k]<-TRUE
      }else
        epresent[k]<-(val[k]==0)   #If zero, skip it; otherwise, add
    }
    if(sum(epresent)>0)               #Adjust attributes for extant edges
      x<-set.edge.attribute(x,names.eval,valsl,eid)
    if(add.edges&&(sum(!epresent)>0))           #Add new edges, if needed
      x<-add.edges(x,as.list(el[!epresent,1]),as.list(el[!epresent,2]), names.eval=as.list(rep(names.eval,sum(!epresent))),vals.eval=as.list(val[!epresent]))
  }
  #Return the modified graph
  x
}


"%e%"<-function(x,attrname){
  get.edge.value(x,attrname=attrname)
}


"%e%<-"<-function(x,attrname,value){
  set.edge.value(x,attrname=attrname,value=value)
}


"%eattr%"<-function(x,attrname){
  x %e% attrname
}


"%eattr%<-"<-function(x,attrname,value){
  x %e% attrname <- value
}


"%n%"<-function(x,attrname){
  get.network.attribute(x,attrname=attrname)
}


"%n%<-"<-function(x,attrname,value){
  set.network.attribute(x,attrname=attrname,value=value)
}


"%nattr%"<-function(x,attrname){
  x %n% attrname
}


"%nattr%<-"<-function(x,attrname,value){
  x %n% attrname <- value
}


"%s%"<-function(x,v){
  if(is.list(v))
    get.inducedSubgraph(x,v=v[[1]],alters=v[[2]])
  else
    get.inducedSubgraph(x,v=v)
}


"%v%"<-function(x,attrname){
  get.vertex.attribute(x,attrname=attrname)
}


"%v%<-"<-function(x,attrname,value){
  set.vertex.attribute(x,attrname=attrname,value=value)
}


"%vattr%"<-function(x,attrname){
  x %v% attrname
}


"%vattr%<-"<-function(x,attrname,value){
  x %v% attrname <- value
}

#"+"<-function(e1, e2, ...) UseMethod("+")
#
#"+.default"<-function(e1,e2,...) { (base::"+")(e1,e2) }
#
#"+.network"<-function(e1,e2,attrname=NULL,...){
#  e1<-as.sociomatrix(e1,attrname=attrname)
#  e2<-as.sociomatrix(e2,attrname=attrname)
#  network(e1+e2,ignore.eval=is.null(attrname),names.eval=attrname)
#}
"+.network"<-function(e1,e2){
  #Set things up
  outinf<-networkOperatorSetup(x=e1,y=e2)
  #Select edges to add; semantics are "adding" edges, which is like union
  #in the non-multigraph case, but actually results in accumulating edge copies
  #in for multiplex graphs.
  out<-outinf$net
  if(is.hyper(out)){            #Hypergraph; for now, return an apology
    stop("Elementwise operations on hypergraphs not yet supported.") 
  }else{                        #Dyadic network
    out<-outinf$net
    #For boolean addition, take the union of edge sets
    el<-rbind(outinf$elx,outinf$ely)
    elna<-rbind(outinf$elnax,outinf$elnay)
    if(!is.multiplex(out)){           #If not multiplex, remove duplicates
      el<-unique(el)
      elna<-unique(elna)
      if(NROW(el)*NROW(elna)>0){
        n<-network.size(out)
        elnum<-(el[,1]-1)+n*(el[,2]-1)
        elnanum<-(elna[,1]-1)+n*(elna[,2]-1)
        elna<-elna[!(elnanum%in%elnum),,drop=FALSE] #For union, NA loses
      }
    }
    if(NROW(el)>0)                             #Add non-missing edges
      add.edges(out,tail=el[,1],head=el[,2])
    if(NROW(elna)>0)                           #Add missing edges
      add.edges(out,tail=elna[,1],head=elna[,2], names.eval=replicate(NROW(elna),list("na")), vals.eval=replicate(NROW(elna),list(list(na=TRUE))))
  }
  #Return the resulting network
  out
}


#"-"<-function(e1, e2, ...) UseMethod("-")
#
#"-.default"<-function(e1,e2,...) { (base::"-")(e1,e2) }
#
"-.network"<-function(e1,e2){
  #Set things up
  outinf<-networkOperatorSetup(x=e1,y=e2)
  #Semantics here are "edge subtraction"; this is like "and not" for the
  #non-multiplex case, but in the latter we can think of it as subtracting
  #copies of edges (so if there were 5 copies of (i,j) in e1 and 2 copies in 
  #e2, we would be left with 3 copies).  Note that this means that NAs are
  #asymmetric: an edge in e2 will first cancel a "sure" edge, and then an
  #NA edge when the sure ones are exhausted.  NA edges in e2 don't cancel
  #sure edges in e1, but they render them unsure (i.e., NA).  NAs in e2
  #have no effect on remaining NAs in e1 (unsure vs unsure), nor on 0s.
  out<-outinf$net
  if(is.hyper(out)){            #Hypergraph; for now, return an apology
    stop("Elementwise operations on hypergraphs not yet supported.") 
  }else{                        #Dyadic network
    out<-outinf$net
    #For boolean subtraction, want edges in e1 that are not in e2
    el<-outinf$elx
    elna<-outinf$elnax
    if(!is.multiplex(out)){         #If not multiplex, cancellation is absolute
      n<-network.size(out)
      elnum<-(el[,1]-1)+n*(el[,2]-1)
      elnanum<-(elna[,1]-1)+n*(elna[,2]-1)
      elynum<-(outinf$ely[,1]-1)+n*(outinf$ely[,2]-1)
      elynanum<-(outinf$elnay[,1]-1)+n*(outinf$elnay[,2]-1)
      #For every edge or NA edge in x, kill it if in ely
      sel<-!(elnum%in%elynum)
      el<-el[sel,,drop=FALSE]
      elnum<-elnum[sel]
      sel<-!(elnanum%in%elynum)
      elna<-elna[sel,,drop=FALSE]
      elnanum<-elnanum[sel]
      #Now, for the remaining edges from x, set to NA if in elyna
      sel<-!(elnum%in%elynanum)
      elna<-rbind(elna,el[!sel,,drop=FALSE])
      el<-el[sel,,drop=FALSE]
      #Clean up any non-uniqueness (recall that el, elna started unique)
      elna<-unique(elna)
    }else{                          #If multiplex, cancellation is 1:1
      n<-network.size(out)
      elnum<-(el[,1]-1)+n*(el[,2]-1)
      elnanum<-(elna[,1]-1)+n*(elna[,2]-1)
      elynum<-(outinf$ely[,1]-1)+n*(outinf$ely[,2]-1)
      elynanum<-(outinf$elnay[,1]-1)+n*(outinf$elnay[,2]-1)
      #Every edge in ely kills one copy of the corresponding edge in el
      i<-1
      while((NROW(el)>0)&&(i<=length(elynum))){
        j<-match(elynum[i],elnum)
        if(is.na(j)){                #No match; increment i
          i<-i+1
        }else{                       #Match!  Cancel both and don't increment
          el<-el[-j,,drop=FALSE]
          elnum<-elnum[-j]
          elynum<-elynum[-i]
        }
      }
      #Every remaining ely kills one copy of the corresponding edge in elna
      i<-1
      while((NROW(elna)>0)&&(i<=length(elynum))){
        j<-match(elynum[i],elnanum)
        if(is.na(j)){                #No match; increment i
          i<-i+1
        }else{                       #Match!  Cancel both and don't increment
          elna<-elna[-j,,drop=FALSE]
          elnanum<-elnanum[-j]
          elynum<-elynum[-i]
        }
      }
      #Every elnay converts one corresponding el to elna
      i<-1
      while((NROW(el)>0)&&(i<=length(elynanum))){
        j<-match(elynanum[i],elnum)
        if(is.na(j)){                #No match; increment i
          i<-i+1
        }else{                       #Match!  Cancel both and don't increment
          elna<-rbind(elna,el[j,,drop=FALSE])
          el<-el[-j,,drop=FALSE]
          elnum<-elnum[-j]
          elynanum<-elynanum[-i]
        }
      }
     }
    if(NROW(el)>0)                             #Add non-missing edges
      add.edges(out,tail=el[,1],head=el[,2])
    if(NROW(elna)>0)                           #Add missing edges
      add.edges(out,tail=elna[,1],head=elna[,2], names.eval=replicate(NROW(elna),list("na")), vals.eval=replicate(NROW(elna),list(list(na=TRUE))))
  }
  #Return the resulting network
  out
}


#"*"<-function(e1, e2, ...) UseMethod("*")
#
#"*.default"<-function(e1,e2,...) { (base::"*")(e1,e2) }
#
"*.network"<-function(e1,e2){
  #Set things up
  outinf<-networkOperatorSetup(x=e1,y=e2)
  #Multiplication semantics here are like "and" in the non-multiplex case,
  #but in the multiplex case we assume that the number of edges is itself
  #multplied.  Multiplication is treated by pairing, so the number of sure
  #edges is sure(e1)*sure(e2), and the number of NA edges is
  #sure(e1)*NA(e2) + NA(e1)*sure(e2) + NA(e1)*NA(e2), where sure and NA are
  #here counts of the (i,j) edge that are non-missing or missing
  #(respectively).
  out<-outinf$net
  if(is.hyper(out)){            #Hypergraph; for now, return an apology
    stop("Elementwise operations on hypergraphs not yet supported.") 
  }else{                        #Dyadic network
    out<-outinf$net
    n<-network.size(out)
    el<-matrix(nrow=0,ncol=2)
    elna<-matrix(nrow=0,ncol=2)
    if(is.multiplex(out)){        #Multiplex case: add edge for every pair
      allpairs<-unique(rbind(outinf$elx,outinf$elnax,outinf$ely,outinf$elnay))
      allnum<-(allpairs[,1]-1)+n*(allpairs[,2]-1)
      elxnum<-(outinf$elx[,1]-1)+n*(outinf$elx[,2]-1)
      elxnanum<-(outinf$elnax[,1]-1)+n*(outinf$elnax[,2]-1)
      elynum<-(outinf$ely[,1]-1)+n*(outinf$ely[,2]-1)
      elynanum<-(outinf$elnay[,1]-1)+n*(outinf$elnay[,2]-1)
      allxcnt<-sapply(allnum,function(z,w){sum(z==w)},w=elxnum)
      allxnacnt<-sapply(allnum,function(z,w){sum(z==w)},w=elxnanum)
      allycnt<-sapply(allnum,function(z,w){sum(z==w)},w=elynum)
      allynacnt<-sapply(allnum,function(z,w){sum(z==w)},w=elynanum)
      el<-allpairs[rep(1:length(allnum),times=allxcnt*allycnt),,drop=FALSE]
      elna<-allpairs[rep(1:length(allnum),times=allxcnt*allynacnt+ allxnacnt*allycnt+allxnacnt*allynacnt),,drop=FALSE]
    }else{                       #Non-multiplex case: "and"
      elx<-unique(outinf$elx)
      elnax<-unique(outinf$elnax)
      ely<-unique(outinf$ely)
      elnay<-unique(outinf$elnay)
      elxnum<-(elx[,1]-1)+n*(elx[,2]-1)
      elxnanum<-(elnax[,1]-1)+n*(elnax[,2]-1)
      sel<-elxnanum%in%elxnum                    #Override NA with edges w/in x
      if(sum(sel)>0){
        elnax<-elnax[!sel,,drop=FALSE]
        elxnanum<-elxnanum[!sel,,drop=FALSE]
      }
      elynum<-(ely[,1]-1)+n*(ely[,2]-1)
      elynanum<-(elnay[,1]-1)+n*(elnay[,2]-1)
      sel<-elynanum%in%elynum                    #Override NA with edges w/in y
      if(sum(sel)>0){
        elnay<-elnay[!sel,,drop=FALSE]
        elynanum<-elynanum[!sel,,drop=FALSE]
      }
      #Check for matches across the "sure" edges
      ematch<-match(elxnum,elynum)
      el<-rbind(el,elx[!is.na(ematch),,drop=FALSE])
      elx<-elx[is.na(ematch),,drop=FALSE]              #Remove the matched cases
      elxnum<-elxnum[is.na(ematch)]
      if(length(ematch[!is.na(ematch)])>0){
        ely<-ely[-ematch[!is.na(ematch)],,drop=FALSE]
        elynum<-elynum[-ematch[!is.na(ematch)]]
      }
      #Match sure xs with unsure ys
      if(length(elxnum)*length(elynanum)>0){
        ematch<-match(elxnum,elynanum)
        elna<-rbind(elna,elx[!is.na(ematch),,drop=FALSE])
        elx<-elx[is.na(ematch),,drop=FALSE]            #Remove the matched cases
        elxnum<-elxnum[is.na(ematch)]
        if(length(ematch[!is.na(ematch)])>0){
          elnay<-elnay[-ematch[!is.na(ematch)],,drop=FALSE]
          elynanum<-elynanum[-ematch[!is.na(ematch)]]
        }
      }
      #Match sure ys with unsure xs
      if(length(elynum)*length(elxnanum)>0){
        ematch<-match(elynum,elxnanum)
        elna<-rbind(elna,ely[!is.na(ematch),,drop=FALSE])
        ely<-ely[is.na(ematch),,drop=FALSE]            #Remove the matched cases
        elynum<-elynum[is.na(ematch)]
        if(length(ematch[!is.na(ematch)])>0){
          elnax<-elnax[-ematch[!is.na(ematch)],,drop=FALSE]
          elxnanum<-elxnanum[-ematch[!is.na(ematch)]]
        }
      }
      #Match unsure xs with unsure ys
      if(length(elxnanum)*length(elynanum)>0){
        ematch<-match(elxnanum,elynanum)
        elna<-rbind(elna,elnax[!is.na(ematch),,drop=FALSE])
      }
    }
    if(NROW(el)>0)                             #Add non-missing edges
      add.edges(out,tail=el[,1],head=el[,2])
    if(NROW(elna)>0)                           #Add missing edges
      add.edges(out,tail=elna[,1],head=elna[,2], names.eval=replicate(NROW(elna),list("na")), vals.eval=replicate(NROW(elna),list(list(na=TRUE))))
  }
  #Return the resulting network
  out
}


"!.network"<-function(e1){
  #Set things up
  outinf<-networkOperatorSetup(x=e1)
  #Select edges to add; semantics are "not" which means that one takes the
  #non-multiplex complement of edges.  Any sure edge implies 0, an NA edge
  #without a sure edge implies NA, no sure or NA edge implies 1.
  out<-outinf$net
  if(is.hyper(out)){            #Hypergraph; for now, return an apology
    stop("Elementwise operations on hypergraphs not yet supported.") 
  }else{                        #Dyadic network
    out<-outinf$net
    n<-network.size(out)
    #Start with the complete graph, and cut things away
    el<-cbind(rep(1:n,each=n),rep(1:n,n))
    if(!is.directed(out))         #Needs to match order in networkOperatorSetup
      el<-el[el[,1]<=el[,2],]
    if(!has.loops(out))
      el<-el[el[,1]!=el[,2],]
    elnum<-(el[,1]-1)+n*(el[,2]-1)
    elna<-matrix(nrow=0,ncol=2)
    #Remove all sure edges
    elx<-unique(outinf$elx)
    elxnum<-(elx[,1]-1)+n*(elx[,2]-1)
    ematch<-match(elxnum,elnum)
    if(length(ematch[!is.na(ematch)])>0){
      el<-el[-ematch[!is.na(ematch)],,drop=FALSE]
      elnum<-elnum[-ematch[!is.na(ematch)]]
    }
    #Convert all unsure edges to NAs
    elnax<-unique(outinf$elnax)
    elxnanum<-(elnax[,1]-1)+n*(elnax[,2]-1)
    ematch<-match(elxnanum,elnum)
    if(length(ematch[!is.na(ematch)])>0){
      elna<-el[ematch[!is.na(ematch)],,drop=FALSE]
      el<-el[-ematch[!is.na(ematch)],,drop=FALSE]
    }
    if(NROW(el)>0)                             #Add non-missing edges
      add.edges(out,tail=el[,1],head=el[,2])
    if(NROW(elna)>0)                           #Add missing edges
      add.edges(out,tail=elna[,1],head=elna[,2], names.eval=replicate(NROW(elna),list("na")), vals.eval=replicate(NROW(elna),list(list(na=TRUE))))
  }
  #Return the resulting network
  out
}


"|.network"<-function(e1,e2){
  #Set things up
  outinf<-networkOperatorSetup(x=e1,y=e2)
  #Select edges to add; semantics are "or," which means that one takes the
  #non-multiplex union of edges (like the non-multiplex case of the +
  #operator).  Here, a sure edge in either input graph will override an NA,
  #and an NA will override a zero.
  out<-outinf$net
  if(is.hyper(out)){            #Hypergraph; for now, return an apology
    stop("Elementwise operations on hypergraphs not yet supported.") 
  }else{                        #Dyadic network
    out<-outinf$net
    #For boolean addition, take the union of edge sets
    el<-rbind(outinf$elx,outinf$ely)
    elna<-rbind(outinf$elnax,outinf$elnay)
    el<-unique(el)
    elna<-unique(elna)
    if(NROW(el)*NROW(elna)>0){
      n<-network.size(out)
      elnum<-(el[,1]-1)+n*(el[,2]-1)
      elnanum<-(elna[,1]-1)+n*(elna[,2]-1)
      elna<-elna[!(elnanum%in%elnum),,drop=FALSE] #For union, NA loses
    }
    if(NROW(el)>0)                             #Add non-missing edges
      add.edges(out,tail=el[,1],head=el[,2])
    if(NROW(elna)>0)                           #Add missing edges
      add.edges(out,tail=elna[,1],head=elna[,2], names.eval=replicate(NROW(elna),list("na")), vals.eval=replicate(NROW(elna),list(list(na=TRUE))))
  }
  #Return the resulting network
  out
}


"&.network"<-function(e1,e2){
  #Set things up
  outinf<-networkOperatorSetup(x=e1,y=e2)
  #Select edges to add; semantics are "and," which means that one places an
  #(i,j) edge if there exists a sure (i,j) edge in both e1 and e2.  If there
  #is not a sure edge in each but there is at least an unsure edge in each,
  #then we place an NA in the (i,j) slot.  Otherwise, we leave it empty.  This
  #is just like boolean "and" for non-multiplex graphs, but is not quite the
  #same in the multiplex case.
  out<-outinf$net
  if(is.hyper(out)){            #Hypergraph; for now, return an apology
    stop("Elementwise operations on hypergraphs not yet supported.") 
  }else{                        #Dyadic network
    out<-outinf$net
    n<-network.size(out)
    el<-matrix(nrow=0,ncol=2)
    elna<-matrix(nrow=0,ncol=2)
    elx<-unique(outinf$elx)
    elnax<-unique(outinf$elnax)
    ely<-unique(outinf$ely)
    elnay<-unique(outinf$elnay)
    elxnum<-(elx[,1]-1)+n*(elx[,2]-1)
    elxnanum<-(elnax[,1]-1)+n*(elnax[,2]-1)
    sel<-elxnanum%in%elxnum                    #Override NA with edges w/in x
    if(sum(sel)>0){
      elnax<-elnax[!sel,,drop=FALSE]
      elxnanum<-elxnanum[!sel,,drop=FALSE]
    }
    elynum<-(ely[,1]-1)+n*(ely[,2]-1)
    elynanum<-(elnay[,1]-1)+n*(elnay[,2]-1)
    sel<-elynanum%in%elynum                    #Override NA with edges w/in y
    if(sum(sel)>0){
      elnay<-elnay[!sel,,drop=FALSE]
      elynanum<-elynanum[!sel,,drop=FALSE]
    }
    #Check for matches across the "sure" edges
    ematch<-match(elxnum,elynum)
    el<-rbind(el,elx[!is.na(ematch),,drop=FALSE])
    elx<-elx[is.na(ematch),,drop=FALSE]              #Remove the matched cases
    elxnum<-elxnum[is.na(ematch)]
    if(length(ematch[!is.na(ematch)])>0){
      ely<-ely[-ematch[!is.na(ematch)],,drop=FALSE]
      elynum<-elynum[-ematch[!is.na(ematch)]]
    }
    #Match sure xs with unsure ys
    if(length(elxnum)*length(elynanum)>0){
      ematch<-match(elxnum,elynanum)
      elna<-rbind(elna,elx[!is.na(ematch),,drop=FALSE])
      elx<-elx[is.na(ematch),,drop=FALSE]            #Remove the matched cases
      elxnum<-elxnum[is.na(ematch)]
      if(length(ematch[!is.na(ematch)])>0){
        elnay<-elnay[-ematch[!is.na(ematch)],,drop=FALSE]
        elynanum<-elynanum[-ematch[!is.na(ematch)]]
      }
    }
    #Match sure ys with unsure xs
    if(length(elynum)*length(elxnanum)>0){
      ematch<-match(elynum,elxnanum)
      elna<-rbind(elna,ely[!is.na(ematch),,drop=FALSE])
      ely<-ely[is.na(ematch),,drop=FALSE]            #Remove the matched cases
      elynum<-elynum[is.na(ematch)]
      if(length(ematch[!is.na(ematch)])>0){
        elnax<-elnax[-ematch[!is.na(ematch)],,drop=FALSE]
        elxnanum<-elxnanum[-ematch[!is.na(ematch)]]
      }
    }
    #Match unsure xs with unsure ys
    if(length(elxnanum)*length(elynanum)>0){
      ematch<-match(elxnanum,elynanum)
      elna<-rbind(elna,elnax[!is.na(ematch),,drop=FALSE])
    }
    if(NROW(el)>0)                             #Add non-missing edges
      add.edges(out,tail=el[,1],head=el[,2])
    if(NROW(elna)>0)                           #Add missing edges
      add.edges(out,tail=elna[,1],head=elna[,2], names.eval=replicate(NROW(elna),list("na")), vals.eval=replicate(NROW(elna),list(list(na=TRUE))))
  }
  #Return the resulting network
  out
}


# --------------------------- %c% -------------------------------
# conditionally create this method, as it may allready have 
# been created and loaded by sna package

if (!exists('%c%')){

"%c%"<-function(e1,e2){
  UseMethod("%c%",e1)
}

}


"%c%.network"<-function(e1,e2){
  #Set things up
  net1<-networkOperatorSetup(x=e1)
  net2<-networkOperatorSetup(x=e2)
  if(is.bipartite(net1$net)){          #Find in/out set sizes for e1
    insz1<-net1$net%n%"bipartite"
    outsz1<-net1$net%n%"n"-net1$net%n%"bipartite"
  }else{
    insz1<-net1$net%n%"n"
    outsz1<-net1$net%n%"n"
  }
  if(is.bipartite(net2$net)){          #Find in/out set sizes for e2
    insz2<-net2$net%n%"bipartite"
    outsz2<-net2$net%n%"n"-net2$net%n%"bipartite"
  }else{
    insz2<-net2$net%n%"n"
    outsz2<-net2$net%n%"n"
  }
  if(outsz1!=insz2)
    stop("Non-conformable relations in %c%.  Cannot compose.")
  if(is.hyper(net1$net)||is.hyper(net2$net))  #Hypergraph; for now, stop
    stop("Elementwise operations on hypergraphs not yet supported.")
  #Test for vertex name matching (governs whether we treat as bipartite)
  if(is.network(e1))
    vnam1<-network.vertex.names(e1)
  else if(!is.null(attr(e1,"vnames")))
    vnam1<-attr(e1,"vnames")
  else if(is.matrix(e1)||is.data.frame(e1)||is.array(e1))
    vnam1<-row.names(e1)
  else
    vnam1<-NULL
  if(is.network(e2))
    vnam2<-network.vertex.names(e2)
  else if(!is.null(attr(e2,"vnames")))
    vnam2<-attr(e2,"vnames")
  else if(is.matrix(e2)||is.data.frame(e2)||is.array(e2))
    vnam2<-row.names(e2)
  else
    vnam2<-NULL
  if((!is.null(vnam1))&&(!is.null(vnam2))&&(length(vnam1)==length(vnam2)) &&all(vnam1==vnam2))
    vnammatch<-TRUE
  else
    vnammatch<-FALSE
  #Decide on bipartite representation and create graph
  if((!is.bipartite(net1$net))&&(!is.bipartite(net2$net))&&vnammatch)
    out<-network.initialize(insz1, directed=is.directed(net1$net)|is.directed(net2$net), loops=TRUE,multiple=is.multiplex(net1$net)|is.multiplex(net2$net))
  else
    out<-network.initialize(insz1+outsz2,bipartite=insz1, directed=is.directed(net1$net)|is.directed(net2$net),multiple=is.multiplex(net1$net)|is.multiplex(net2$net))
  #Accumulate edges (yeah, could be made more efficient -- cope with it)
  el<-matrix(nrow=0,ncol=2)
  elna<-matrix(nrow=0,ncol=2)
  bip1<-net1$net%n%"bipartite"
  bip2<-net2$net%n%"bipartite"
  if(!is.directed(net1$net)){  #Double the edges if undirected
    net1$elx<-rbind(net1$elx,net1$elx[net1$elx[,1]!=net1$elx[,2],2:1])
    net1$elnax<-rbind(net1$elnax,net1$elnax[net1$elnax[,1]!=net1$elnax[,2],2:1])
  }
  if(!is.directed(net2$net)){  #Double the edges if undirected
    net2$elx<-rbind(net2$elx,net2$elx[net2$elx[,1]!=net2$elx[,2],2:1])
    net2$elnax<-rbind(net2$elnax,net2$elnax[net2$elnax[,1]!=net2$elnax[,2],2:1])
  }
  if(NROW(net1$elx)>0){
    for(i in 1:NROW(net1$elx)){
      sel<-net2$elx[net2$elx[,1]==(net1$elx[i,2]-bip1),2]-bip2
      if(length(sel)>0)
        el<-rbind(el,cbind(rep(net1$elx[i,1],length(sel)),sel+insz1))
    }
  }
  if(NROW(net1$elnax)>0){
    for(i in 1:NROW(net1$elnax)){
      sel<-net2$elnax[net2$elnax[,1]==(net1$elnax[i,2]-bip1),2]-bip2
      if(length(sel)>0)
        elna<-rbind(elna,cbind(rep(net1$elnax[i,1],length(sel)),sel+insz1))
    }
  }
  if(!is.bipartite(out)){     #If not bipartite, remove the insz1 offset
    if(NROW(el)>0)
      el[,2]<-el[,2]-insz1
    if(NROW(elna)>0)
      elna[,2]<-elna[,2]-insz1
  }
  if(!is.multiplex(out)){     #If necessary, consolidate edges
    if(NROW(el)>1)
      el<-unique(el)
    if(NROW(elna)>1){
      elna<-unique(elna)
    }
    if(NROW(elna)*NROW(el)>0){
      sel<-rep(TRUE,NROW(elna))
      for(i in 1:NROW(elna)){
        if(any((el[,1]==elna[i,1])&(el[,2]==elna[i,2])))
          sel[i]<-FALSE
      }
      elna<-elna[sel,]
    }
  }
  #Add the edges
  if(NROW(el)>0)                             #Add non-missing edges
    add.edges(out,tail=el[,1],head=el[,2])
  if(NROW(elna)>0)                           #Add missing edges
    add.edges(out,tail=elna[,1],head=elna[,2], names.eval=replicate(NROW(elna),list("na")), vals.eval=replicate(NROW(elna),list(list(na=TRUE))))
  #Return the resulting network
  out
}


#Given one or two input networks, return the information needed to generate
#output for binary or unary operations.  The return value for this function is
#a list with elements:
# net: the output network (empty, but with attributes set)
# elx: the edgelist for the first network (non-missing)
# elnax: the list of missing edges for the first network
# ely: in the binary case, the edgelist for the second network (non-missing)
# elnay: in the binary case, the list of missing edges for the second network
networkOperatorSetup<-function(x,y=NULL){
  #Determine what attributes the output should have
  if(is.network(x)){
    nx<-network.size(x)     #Get size, directedness, multiplexity, bipartition
    dx<-is.directed(x)
    mx<-is.multiplex(x)
    hx<-is.hyper(x)
    lx<-has.loops(x)
    bx<-x%n%"bipartite"
    if(is.null(bx))
      bx<-FALSE
  }else{                     #If not a network object, resort to adj form
    x<-as.sociomatrix(x)
    if(NROW(x)!=NCOL(x)){    #Bipartite matrix
      nx<-NROW(x)+NCOL(x)
      dx<-FALSE
      mx<-FALSE
      hx<-FALSE
      lx<-FALSE
      bx<-NROW(x)
    }else{
      nx<-NROW(x)
      dx<-TRUE
      mx<-FALSE
      hx<-FALSE
      lx<-any(diag(x)!=0,na.rm=TRUE)
      bx<-FALSE
    }
  }
  if(is.null(y)){                 #If y is null, setup for unary operator
    n<-nx
    d<-dx
    m<-mx
    h<-hx
    b<-bx
    l<-lx
    x<-x
  }else{                          #Binary case
    if(is.network(y)){
      ny<-network.size(y)     #Get size, directedness, multiplexity, bipartition
      dy<-is.directed(y)
      my<-is.multiplex(y)
      hy<-is.hyper(y)
      ly<-has.loops(y)
      by<-y%n%"bipartite"
      if(is.null(by))
        by<-FALSE
    }else{                     #If not a network object, resort to adj form
      y<-as.sociomatrix(y)
      if(NROW(y)!=NCOL(y)){    #Bipartite matrix
        ny<-NROW(y)+NCOL(y)
        dy<-FALSE
        my<-FALSE
        hy<-FALSE
        ly<-FALSE
        by<-NROW(y)
      }else{
        ny<-NROW(y)
        dy<-TRUE
        my<-FALSE
        hy<-FALSE
        ly<-any(diag(y)!=0,na.rm=TRUE)
        by<-FALSE
      }
    }
    if(nx!=ny)                     #Make sure that our networks are conformable
      stop("Non-conformable networks (must have same numbers of vertices for elementwise operations).")
    if(bx!=by)
      stop("Non-conformable networks (must have same bipartite status for elementwise operations).")
    n<-nx                         #Output size=input size
    b<-bx                         #Output bipartition=input bipartition
    d<-dx|dy                      #Output directed if either input directed
    l<-lx|ly                      #Output has loops if either input does
    h<-hx|hy                      #Output hypergraphic if either input is
    m<-mx|my                      #Output multiplex if either input is
  }
  #Create the empty network object that will ultimately receive the edges
  net<-network.initialize(n=n,directed=d,hyper=h,loops=l,multiple=m,bipartite=b)
  #Create the edge lists; what the operator does with 'em isn't our problem
  if(h){                                           #Hypergraph
    stop("Elementwise operations not yet supported on hypergraphs.")
  }else{                                           #Dyadic network
    #Get the raw edge information
    if(is.network(x)){
      elx<-as.matrix(x,matrix.type="edgelist")
      elnax<-as.matrix(is.na(x),matrix.type="edgelist")
      if(d&(!dx)){      #Need to add two-way edges; BTW, can't have (!d)&dx...
        elx<-rbind(elx,elx[elx[,2]!=elx[,1],2:1,drop=FALSE])
        elnax<-rbind(elnax,elnax[,2:1])
      } else if (!dx){ # need to enforce edge ordering i<j for comparison
        # replace all rows where i<j with j,i
        elx[elx[,1]>elx[,2],]<-elx[elx[,1]>elx[,2],c(2,1)]
      }
    }else{
      elx<-which(x!=0,arr.ind=TRUE)
      elnax<-which(is.na(x),arr.ind=TRUE)
      if(!d){   #Sociomatrix already has two-way edges, so might need to remove
        elx<-elx[elx[,1]>=elx[,2],,drop=FALSE]
        elnax<-elnax[elnax[,1]>=elnax[,2],,drop=FALSE]
      } 
    }
    if(!is.null(y)){
      if(is.network(y)){
        ely<-as.matrix(y,matrix.type="edgelist")
        elnay<-as.matrix(is.na(y),matrix.type="edgelist")
        if(d&(!dy)){      #Need to add two-way edges; BTW, can't have (!d)&dy...
          ely<-rbind(ely,ely[ely[,2]!=ely[,1],2:1,drop=FALSE])
          elnay<-rbind(elnay,elnay[,2:1])
        } else if (!dy){ # need to enforce edge ordering i<j for comparison
          # replace all rows where i<j with j,i
          ely[ely[,1]>ely[,2],]<-ely[ely[,1]>ely[,2],c(2,1)]
        }
      }else{
        ely<-which(y!=0,arr.ind=TRUE)
        elnay<-which(is.na(y),arr.ind=TRUE)
        if(!d){  #Sociomatrix already has two-way edges, so might need to remove
          ely<-ely[ely[,1]>=ely[,2],,drop=FALSE]
          elnay<-elnay[elnay[,1]>=elnay[,2],d,rop=FALSE]
        }
      }
    }
    if(!l){        #Pre-emptively remove loops, as needed
      elx<-elx[elx[,1]!=elx[,2],,drop=FALSE]
      elnax<-elnax[elnax[,1]!=elnax[,2],,drop=FALSE]
      if(!is.null(y)){
        ely<-ely[ely[,1]!=ely[,2],,drop=FALSE]
        elnay<-elnay[elnay[,1]!=elnay[,2],,drop=FALSE]
      }
    }
    if(!m){        #Pre-emptively remove multiplex edges, as needed
      elx<-unique(elx)
      elnax<-unique(elnax)
      if(!is.null(y)){
        ely<-unique(ely)
        elnay<-unique(elnay)
      }
    }
  }
  #Return everything
  if(is.null(y))
    list(net=net,elx=elx,elnax=elnax)
  else
    list(net=net,elx=elx,elnax=elnax,ely=ely,elnay=elnay)
}


prod.network<-function(..., attrname=NULL, na.rm=FALSE){
  inargs<-list(...)
  y<-inargs[[1]]
  for(i in (1:length(inargs))[-1]){
    x<-as.sociomatrix(inargs[[i]],attrname=attrname)
    if(na.rm)
      x[is.na(x)]<-0
    ym<-as.sociomatrix(y,attrname=attrname)
    if(na.rm)
      ym[is.na(ym)]<-0
    y[,,names.eval=attrname,add.edges=TRUE]<-x*ym
  }
  y
}


sum.network<-function(..., attrname=NULL, na.rm=FALSE){
  inargs<-list(...)
  y<-inargs[[1]]
  for(i in (1:length(inargs))[-1]){
    x<-as.sociomatrix(inargs[[i]],attrname=attrname)
    if(na.rm)
      x[is.na(x)]<-0
    ym<-as.sociomatrix(y,attrname=attrname)
    if(na.rm)
      ym[is.na(ym)]<-0
    y[,,names.eval=attrname,add.edges=TRUE]<-x+ym
  }
  y
}
