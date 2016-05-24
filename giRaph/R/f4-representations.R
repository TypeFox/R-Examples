## f4-representations.R --- 
## Author          : Jens Henrik Badsberg, Claus Dethlefsen, Luca La Rocca
## Created On      : Tue Nov 30 16:50:00 2004
## Last Modified By: Luca La Rocca
## Last Modified On: Fri Mar 10 18:36:00 2006
## Update Count    : 51
## Status          : Unknown, Use with caution!
######################################################

## The four representations:
##     G=(V,E): incidenceList
##     I: incidenceMatrix
##     A: adjacencyList
##     X: adjacencyMatrix
##
## Equivalence (~) of these are as follows
##
## anyGraph:     G
## generalGraph: G ~ I
## multiGraph:   G ~ I ~ A
## simpleGraph:  G ~ I ~ A ~ X
##
## where equivalence means that for the given family of graphs
## no information is lost when coercing from one representation to the other.
## Character vertex identifiers are kept as dim-attributes
## on the X and I matrices, and on the A-list.

## construction and visualization

# constructor method for class 'incidenceList'
setMethod("initialize","incidenceList",
                 function(.Object,V=character(0),E=list()){
                     .Object@V<-new("vertexSet",V)
                     .Object@E<-new("edgeList",
                        E[unlist(lapply(E,function(x){isEmpty(x)||maxId(x)<=length(.Object@V)}))])
                     return(.Object)
                 } # end of function
                ) # end of SetMethod
# a valid 'incidenceList' object is returned

# constructor method for class 'incidenceMatrix'
setMethod("initialize","incidenceMatrix",
                 function(.Object,I=matrix(0,nrow=0,ncol=0)){
                     .Object@.Data<-as(I,"matrix")
                     n<-ncol(.Object@.Data) # number of vertices
                     m<-nrow(.Object@.Data) # number of edges
                     if(n>0){ # make valid character vertex identifiers
                        Vnames<-colnames(.Object@.Data)
                        if(is(Vnames,"NULL")||any(duplicated(Vnames)))
                            colnames(.Object)<-make.names(seq(1:n)) # default vertex names
                        else
                            colnames(.Object)<-make.names(Vnames,unique=TRUE) # vertex names from input
                        if(m>0){ # make valid entries
                            for(i in 1:m){ # for all edges
                                nonzeros<-.Object@.Data[i,]!=0
                                    .Object@.Data[i,nonzeros]<-as.numeric(factor(rank(.Object@.Data[i,nonzeros],
                                                                                      ties.method="min")))
                            } # end of for
                            nonempty<-apply(.Object@.Data,1,function(x){sum(x)>0}) #Êfind non-empty edges
                            .Object@.Data<-matrix(.Object@.Data[nonempty,],sum(nonempty),n)
                        } # end of if
                     } # end of if
                     rownames(.Object)<-NULL # there should not be any edge names
                     return(.Object)
                 } # end of function
                ) # end of SetMethod
# a valid 'incidenceMatrix' object is returned

# constructor method for class 'adjacencyList'
setMethod("initialize","adjacencyList",
          function(.Object,id=character(0),pa=list(),ch=list(),ne=list()){
            Vnames<-names(new("vertexSet",id))
            n<-length(Vnames)
            .Object@.Data<-rep(list(list()),n)
            names(.Object@.Data)<-Vnames
            if(n>0){ # there are vertices
              if(length(pa)==n){ # parents are given
                for(i in 1:n){ # add parents
                  aux<-as(pa[[i]],"integer")
                  aux<-aux[(aux>0)&(aux<=n)]
                  if(length(aux)>0){ # there are parents
                    .Object@.Data[[i]]$pa<-c(.Object@.Data[[i]]$pa,aux)
                    for(j in aux) .Object@.Data[[j]]$ch<-c(.Object@.Data[[j]]$ch,i) # fix children
                  } # end of if
                } # end of for
              } # end of if
              if(length(ch)==n){ # children are given
                for(i in 1:n){ # add children
                  aux<-as(ch[[i]],"integer")
                  aux<-aux[(aux>0)&(aux<=n)]
                  if(length(aux)>0){ # there are children
                    .Object@.Data[[i]]$ch<-c(.Object@.Data[[i]]$ch,aux)
                    for(j in aux) .Object@.Data[[j]]$pa<-c(.Object@.Data[[j]]$pa,i) # fix parents
                  } # end of if
                } # end of for
              } # end of if
              if(length(ne)==n){ # neighbours are given
                for(i in 1:n){ # add neighbours
                  aux<-as(ne[[i]],"integer")
                  aux<-aux[(aux>0)&(aux<=n)]
                  if(length(aux)>0){ # there are neighbours
                    for(j in aux){ # for all neighbour occurrences
                     if(j>i){ # not yet considered
                       .Object@.Data[[i]]$ne<-c(.Object@.Data[[i]]$ne,j)
                       .Object@.Data[[j]]$ne<-c(.Object@.Data[[j]]$ne,i)
                     } else if(j<i){ # maybe already considered
                       if((sum(j==aux)-sum(j==.Object@.Data[[i]]$ne))>0){ # more to add
                         .Object@.Data[[i]]$ne<-c(.Object@.Data[[i]]$ne,j)
                         .Object@.Data[[j]]$ne<-c(.Object@.Data[[j]]$ne,i)
                       } # end of if
                     } else .Object@.Data[[i]]$ne<-c(.Object@.Data[[i]]$ne,i) # j==i, that is a loop
                    } # end of for 
                  } # end of if
                } # end of for
              } # end of if
            } # end of if
            return(.Object)
          } # end of function
         ) # end of SetMethod
# a valid 'adjacencyList' object is returned

# constructor method for class 'adjacencyMatrix'
setMethod("initialize","adjacencyMatrix",
                 function(.Object,X=matrix(0,nrow=0,ncol=0)){
                     aux<-as(X,"matrix")
                     n<-min(nrow(aux),ncol(aux))
                     Rnames<-rownames(aux)
                     Cnames<-colnames(aux)
                     if(n>0){ # make valid entries and character vertex identifiers
                        .Object@.Data<-matrix(0,n,n)
                        if(n==1){ # single vertex
                            if(is(Rnames,"NULL")||is(Cnames,"NULL")||(Rnames[1]!=Cnames[1])){
                                rownames(.Object@.Data)<-"X1" # default vertex name
                            }else{ # vertex name from input
                                rownames(.Object@.Data)<-make.names(Rnames[1])
                            } # end of if-else
                        }else{ # two or more vertices
                            for(i in 1:n) for (j in seq(1,n)[-i]) if(aux[i,j]!=0) .Object@.Data[i,j]<-1
                            if(is(Rnames,"NULL")||is(Cnames,"NULL")||
                               any(Rnames[1:n]!=Cnames[1:n])||any(duplicated(Rnames[1:n]))){
                                 rownames(.Object@.Data)<-make.names(seq(1:n)) # default vertex names
                            }else{ #Êvertex names from input
                                rownames(.Object@.Data)<-make.names(Rnames[1:n],unique=TRUE)
                            } # end of if-else
                        } # end of if-else
                        colnames(.Object@.Data)<-rownames(.Object@.Data) # make the same
                     } # end of if
                     return(.Object)
                 } # end of function
                ) # end of SetMethod
# a valid 'adjacencyMatrix' object is returned

# show method for class 'incidenceList'
setMethod("show","incidenceList",
                 function(object){
                     cat("An object of class \"incidenceList\"",fill=T)
                     cat("V=")
                     show(object@V)
                     cat("E=")
                     showRel(object@E,object@V)
                 } # end of function
                ) # end of setMethod
# a shorter representation than the default one

# keeping default show method for class 'incidenceMatrix'

# show method for class 'adjacencyList'
setMethod("show","adjacencyList",
                 function(object){
                     cat("An object of class \"adjacencyList\"",fill=T)
                     if(!isEmpty(object)){ # something to show
                         Vnames<-names(object)
                         for(i in 1:length(object)){
                             blank<-rep(" ",nchar(Vnames[i]))
                             cat(Vnames[i]," <- {",sep="")
                             cat(Vnames[object@.Data[[i]]$pa],sep=",")
                             cat("}",fill=T)
                             cat(blank," -- {",sep="")
                             cat(Vnames[object@.Data[[i]]$ne],sep=",")
                             cat("}",fill=T)
                             cat(blank," -> {",sep="")
                             cat(Vnames[object@.Data[[i]]$ch],sep=",")
                             cat("}",fill=T)
                         } # end of for
                     } else{ # nothing to show
                       cat("list()",fill=T)
                     } # end of if-else
                 } # end of function
                ) # end of setMethod
# a shorter representation than the default one

# keeping default show method for class 'adjacencyMatrix'

## getting and setting information

# 'names' method for class 'incidenceList'
setMethod("names", "incidenceList", function(x) names(x@V))
# take the names from the 'vertexSet'

# 'names<-' replacement method for class 'incidenceList'
setReplaceMethod("names","incidenceList",
                         function(x,value){
                         aux<-new("vertexSet",value)
                         if(card(aux)==card(x@V)) x@V<-aux
                         else warning("Sorry, wrong cardinality: names unchanged...")
                         x # returns the possibly modified object
                         } # end of function
                        ) # end of setMethod
# set the names, if possible (number of vertices cannot be changed)

# 'names' method for class 'incidenceMatrix'
setMethod("names", "incidenceMatrix", function(x) colnames(x))
# take the names of the columns

# 'names<-' replacement method for class 'incidenceMatrix'
setReplaceMethod("names","incidenceMatrix",
                         function(x,value){
                         aux<-new("vertexSet",value)
                         if(card(aux)==ncol(x)) colnames(x)<-aux@.Data
                         else warning("Sorry, wrong cardinality: names unchanged...")
                         x # returns the possibly modified object
                         } # end of function
                        ) # end of setMethod
# set the names, if possible (number of vertices cannot be changed)

# 'names' method for class 'adjacencyList' is inherited from class 'list'

# 'names<-' replacement method for class 'adjacencyList'
setReplaceMethod("names","adjacencyList",
                         function(x,value){
                         aux<-new("vertexSet",value)
                         if(card(aux)==length(x@.Data)) names(x@.Data)<-aux@.Data
                         else warning("Sorry, wrong cardinality: names unchanged...")
                         x # returns the possibly modified object
                         } # end of function
                        ) # end of setMethod
# set the names, if possible (number of vertices cannot be changed)

# 'names' method for class 'adjacencyMatrix'
setMethod("names", "adjacencyMatrix", function(x) colnames(x))
# take the names of the columns

# 'names<-' replacement method for class 'adjacencyMatrix'
setReplaceMethod("names","adjacencyMatrix",
                         function(x,value){
                         aux<-new("vertexSet",value)
                         if(card(aux)==ncol(x)){
                           colnames(x)<-aux@.Data
                           rownames(x)<-aux@.Data
                         } else warning("Sorry, wrong cardinality: names unchanged...")
                         x # returns the possibly modified object
                         } # end of function
                        ) # end of setMethod
# set the names, if possible (number of vertices cannot be changed)

# 'card' method for class 'incidenceList'
setMethod("card", "incidenceList", function(object,...) list(v=card(object@V),e=card(object@E)))
# returns the number of vertices and edges

# 'card' method for class 'incidenceMatrix'
setMethod("card", "incidenceMatrix", function(object,...) list(v=ncol(object),e=nrow(object)))
# returns the number of vertices and edges

# 'card' method for class 'adjacencyList'
setMethod("card", "adjacencyList",
                  function(object,...){
                    n<-length(object@.Data)
                    m<-0
                    if(n>0){ # there are vertices
                      for(i in 1:n){ # add outgoing edges
                        m<-m+length(object@.Data[[i]]$ch)+sum(object@.Data[[i]]$ne>=i)
                      } # end of for
                    } # end of if
                    return(list(v=n,e=m))
                  } # end of function
                ) # end of setMethod
# returns the number of vertices and edges

# 'card' method for class 'adjacencyMatrix'
setMethod("card", "adjacencyMatrix",
                  function(object,...){
                    n<-ncol(object@.Data)
                    m<-0
                    if(n>1) # maybe edges
                      for(i in seq(1,n-1)) # scan rows
                        for(j in seq(i+1,n)) # scan columns
                          if(object@.Data[i,j]) m<-m+1 # edge (directed or undirected)
                          else if (object@.Data[j,i]) m<-m+1 # directed
                    return(list(v=n,e=m))
                  } # end of function
                ) # end of setMethod
# returns the number of vertices and edges

## property checking

# 'isEmpty' method for class 'incidenceList'
setMethod("isEmpty","incidenceList", function(object,...) isEmpty(object@V))
# an 'incidenceList' is empty if such is its 'vertexSet'

# 'isEmpty' method for class 'incidenceMatrix'
setMethod("isEmpty","incidenceMatrix", function(object,...) ncol(object)==0)
# an 'incidenceMatrix' with no columns represents no vertices

# 'isEmpty' method for class 'adiacencyList' is inherited from 'vector'

# 'isEmpty' method for class 'adiacencyMatrix'
setMethod("isEmpty","adjacencyMatrix", function(object,...) nrow(object)==0)
# an 'adjacencyMatrix' is empty if it has no entries

# 'isPresent' method for 'edge' in 'incidenceList'
setMethod("isPresent",c(el="edge",ou="incidenceList"),function(el,ou)return(isPresent(el,ou@E)))
# a 'logical' value answering the question is returned

# 'isPresent' method for 'undirectedEdge' in 'incidenceMatrix'
setMethod("isPresent",c(el="undirectedEdge",ou="incidenceMatrix"),
                      function(el,ou){
                        if(nrow(ou)>0&&!isEmpty(el)&&maxId(el)<=ncol(ou)){ # maybe
                          rowedge<-rep(0,ncol(ou))
                          rowedge[el@.Data]<-1
                          return(any(apply(ou@.Data,1,function(x) all(x==rowedge))))
                        } # end of if
                        else return(FALSE) # surely not
                      } # end of function
                    ) # end of setMethod
# a 'logical' value answering the question is returned

# 'isPresent' method for 'directedEdge' in 'incidenceMatrix'
setMethod("isPresent",c(el="directedEdge",ou="incidenceMatrix"),
                      function(el,ou){
                        if(nrow(ou)>0&&length(el)>1&&maxId(el)<=ncol(ou)){ # maybe
                          rowedge<-rep(0,ncol(ou))
                          rowedge[unlist(el@.Data)]<-rep(1:length(el),unlist(lapply(el@.Data,length)))
                          return(any(apply(ou@.Data,1,function(x) all(x==rowedge))))
                        } # end of if
                        else return(FALSE) # surely not
                      } # end of function
                    ) # end of setMethod
# a 'logical' value answering the question is returned

# 'isPresent' method for 'undirectedEdge' in 'adjacencyList'
setMethod("isPresent",c(el="undirectedEdge",ou="adjacencyList"),
                      function(el,ou){
                        n<-card(ou)$v # number of vertices
                        h<-card(el) # edge cardinality
                        if( n>0 && h>=1 && h<=2 && maxId(el)<=n ){ # maybe
                          return(any(ou@.Data[[el@.Data[1]]]$ne==el@.Data[h]))
                        } else{ # surely not
                          return(FALSE)
                        } # end of if else
                      } # end of function
                    ) # end of setMethod
# a 'logical' value answering the question is returned

# 'isPresent' method for 'directedEdge' in 'adjacencyList'
setMethod("isPresent",c(el="directedEdge",ou="adjacencyList"),
                      function(el,ou){
                        n<-card(ou)$v # number of vertices
                        if( n>0 && length(el)>1 && card(el)<=2 && maxId(el)<=n ){ # maybe
                          return(any(ou@.Data[[el@.Data[[1]]]]$ch==el@.Data[[2]]))
                        } else{ # surely not
                          return(FALSE)
                        } # end of if else
                      } # end of function
                    ) # end of setMethod
# a 'logical' value answering the question is returned

# 'isPresent' method for 'undirectedEdge' in 'adjacencyMatrix'
setMethod("isPresent",c(el="undirectedEdge",ou="adjacencyMatrix"),
                      function(el,ou){
                        n<-card(ou)$v # number of vertices
                        h<-card(el) # edge cardinality
                        if( n>1 && h==2 && maxId(el)<=n ){ # maybe
                          return(ou@.Data[el@.Data[1],el@.Data[2]]&&ou@.Data[el@.Data[2],el@.Data[1]])
                        } else{ # surely not
                          return(FALSE)
                        } # end of if else
                      } # end of function
                    ) # end of setMethod
# a 'logical' value answering the question is returned

# 'isPresent' method for 'directedEdge' in 'adjacencyMatrix'
setMethod("isPresent",c(el="directedEdge",ou="adjacencyMatrix"),
                      function(el,ou){
                        n<-card(ou)$v # number of vertices
                        h<-card(el) # edge cardinality
                        if( n>1 && length(el)>=2 && h<=2 && maxId(el)<=n ){ # maybe
                          return(ou@.Data[el@.Data[[1]],el@.Data[[2]]]&&!ou@.Data[el@.Data[[2]],el@.Data[[1]]])
                        } else{ # surely not
                          return(FALSE)
                        } # end of if else
                      } # end of function
                    ) # end of setMethod
# a 'logical' value answering the question is returned

# comparison method for class 'incidenceList'
setMethod("areTheSame",c("incidenceList","incidenceList"),
                 function(x,y){
                     res<-(areTheSame(x@V,y@V))&&(length(x@E)==length(y@E))
                     if(res) # maybe
                        res<-areTheSame(x@E,recode(y@E,y@V,x@V))
                     return(res)
                 } #Êend of function
         ) # end of setMethod
# a 'logical' value answering the question is returned

# comparison method for class 'incidenceMatrix'
setMethod("areTheSame",c("incidenceMatrix","incidenceMatrix"),
                 function(x,y){
                    res<-(setequal(names(x),names(y)))&&(nrow(x)==nrow(y))
                    if(res&&nrow(x)>0){ # maybe and non-trivial
                        m<-nrow(x) # number of edges
                        y@.Data<-matrix(y@.Data[,names(x)],nrow=m,ncol=card(x)$v,dimnames=list(NULL,names(x)))
                        unmatched<-rep(TRUE,m) # refers to 'y'
                        for(i in 1:m){ # match 'x@.Data[i,]' in 'y'
                            found<-FALSE # refers to 'x@.Data[i,]'
                            for(j in which(unmatched)){ # try all unmatched 'y@.Data[j,]'
                                if(all(x@.Data[i,]==y@.Data[j,])){ # matched
                                    found<-T
                                    unmatched[j]<-FALSE
                                    break
                                 } # end of if
                            } # end of for (j)
                            if(!found){ # not matched
                                res<-FALSE
                                break
                            } # end of if
                        } # end of for (i)
                     res<-res&&!any(unmatched)
                    } # end of if
                    return(res)
                 } # end of function
         ) # end of setMethod
# a 'logical' value answering the question is returned

# comparison method for class 'adjacencyList'
setMethod("areTheSame",c("adjacencyList","adjacencyList"),
                 function(x,y){
                    res<-FALSE
                    if(setequal(names(x),names(y))){ # maybe
                        y@.Data<-y@.Data[match(names(x),names(y))]
                        res<-rep(NA,length(x))
                        for(i in 1:length(res)){
                            res[i]<-setequal(x@.Data[[i]]$ne,y@.Data[[i]]$ne)&&
                            setequal(x@.Data[[i]]$pa,y@.Data[[i]]$pa)&&
                            setequal(x@.Data[[i]]$ch,y@.Data[[i]]$ch)
                        } # end of for
                        res<-all(res)
                    } #Êend of if
                    return(res)
                 } #Êend of function
         ) # end of setMethod
# a 'logical' value answering the question is returned

# comparison method for class 'adjacencyMatrix'
setMethod("areTheSame",c("adjacencyMatrix","adjacencyMatrix"),
                 function(x,y){
                    Xnames<-names(x)
                    res<-setequal(Xnames,names(y))
                    if(res&&nrow(x)>0){ # maybe and non-trivial
                        res<-all(x@.Data==y@.Data[Xnames,Xnames])
                    } # end of if
                    return(res)
                 } # end of function
         ) # end of setMethod
# a 'logical' value answering the question is returned

## extraction

# multi extractor method for class 'incidenceList'
setMethod("[","incidenceList",
          function(x,i,j=NA,drop=NA){
            oldV<-x@V
            x@V<-x@V[i]
            if(length(x@E)>0){ # there are edges
              auxE<-recode(x@E,oldV,x@V)
              noshorter<-logical(0) # edges referring to dropped vertices have been shortened by 'recode'
              for(i in 1:length(auxE))
                noshorter[i]<-(card(auxE[[i]])==card(x@E[[i]]))
              x@E<-auxE[noshorter]
            } # end of if
            return(x)
          } # end of function
         ) # end of setMethod
# the subgraph induced by 'i' is extracted

# single extractor method for class 'incidenceList'
setMethod("[[","incidenceList",
          function(x,i,j=NA,drop=NA) return(x@V[[i]])
         ) # end of setMethod
# the name of the i-th vertex is extracted

# multi extractor method for class 'incidenceMatrix'
setMethod("[","incidenceMatrix",
          function(x,i,j=NA,drop=NA){
            Vnames<-names(x)[i]
            if(length(Vnames)==0){ # empty output
              return(new("incidenceMatrix"))
            }else{
              whichEdges<-apply(matrix(x@.Data[,-i],nrow=card(x)$e,ncol=card(x)$v-length(Vnames)),
                                1,function(y) sum(y)==0)
              return(new("incidenceMatrix",matrix(x@.Data[whichEdges,i],nrow=sum(whichEdges),
                                                  ncol=length(Vnames),dimnames=list(NULL,Vnames))))
            } # end of if-else
          } # end of function
         ) # end of setMethod
# the subgraph induced by 'i' is extracted

# single extractor method for class 'incidenceMatrix'
setMethod("[[","incidenceMatrix",
          function(x,i,j=NA,drop=NA) return(names(x)[i])
         ) # end of setMethod
# the name of the i-th vertex is extracted

# multi extractor method for class 'adjacencyList'
setMethod("[","adjacencyList",
          function(x,i,j=NA,drop=NA){
            if(isEmpty(x)) return(x)
            else{
              A<-new("adjacencyList",id=names(x)[i]) # new empty adjacency list
              n<-length(A) # number of vertices in output
              if(n>0){ # there are vertices in output
                from<-names(x)
                to<-names(A)
                for(i in 1:n){ # deal with each of them
                    j<-match(to[i],from) # find its original position
                    aux<-match(from[x@.Data[[j]]$ne],to)
                    A@.Data[[i]]$ne<-aux[!is.na(aux)]
                    aux<-match(from[x@.Data[[j]]$pa],to)
                    A@.Data[[i]]$pa<-aux[!is.na(aux)]
                    aux<-match(from[x@.Data[[j]]$ch],to)
                    A@.Data[[i]]$ch<-aux[!is.na(aux)]
                } # end of for
              } # end of if
            } # end of if else
            return(A)
          } # end of function
         ) # end of setMethod
# the subgraph induced by 'i' is extracted

# single extractor method for class 'adjacencyList'
setMethod("[[","adjacencyList",
          function(x,i,j=NA,drop=NA) return(names(x)[i])
         ) # end of setMethod
# the name of the i-th vertex is extracted

# multi extractor method for class 'adjacencyMatrix'
setMethod("[","adjacencyMatrix",
          function(x,i,j=NA,drop=NA){
            Vnames<-names(x)[i]
            if(length(Vnames)==0) # empty output
              return(new("adjacencyMatrix"))
            else
              return(new("adjacencyMatrix",matrix(x@.Data[i,i],nrow=length(Vnames),
                                                  ncol=length(Vnames),dimnames=list(Vnames,Vnames))))
          } # end of function
         ) # end of setMethod
# the subgraph induced by 'i' is extracted

# single extractor method for class 'adjacencyMatrix'
setMethod("[[","adjacencyMatrix",
          function(x,i,j=NA,drop=NA) return(names(x)[i])
         ) # end of setMethod
# the name of the i-th vertex is extracted

## typecasting

# see file 'f5-conversions.R' for coerce methods

## operators

# see file 'f7-operators.R' for '+/-/*' methods
