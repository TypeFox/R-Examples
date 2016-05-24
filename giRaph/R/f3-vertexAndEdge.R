## f3-vertexAndEdge.R --- 
## Author          : Jens Henrik Badsberg, Claus Dethlefsen, Luca La Rocca
## Created On      : Tue Nov 30 14:32:00 2004
## Last Modified By: Luca La Rocca
## Last Modified On: Fri Feb 15 13:18:00 2008
## Update Count    : 28
## Status          : Unknown, Use with caution!
######################################################

### ----------------------------------vertices---------------------------------------

## construction and visualization

# constructor method for class 'vertexSet'
setMethod("initialize","vertexSet",
                 function(.Object,...){
                     id<-unique(as(unlist(list(...)),"character")) # get unique 'character' values
                     if(length(id)>0)
                         .Object@.Data<-make.names(id,unique=TRUE) # make unique valid names
                     return(.Object)
                 } # end of function
                ) # end of SetMethod
# 'new("vertexSet",x,y)' returns a 'vertexSet' object consisting of a vector of
# unique syntactically valid names made from the unique elements of 'x' and 'y'

# wrapper for 'vertexSet' construction
v <- function(...) new("vertexSet",...)
# just call the right initialize method

# show method for class 'vertexSet'
setMethod("show","vertexSet",
                 function(object)
                     cat("{",paste(object@.Data,collapse=","),"}",sep="",fill=T)
                ) # end of setMethod
# for example '{X1,X2,X3}' denotes a three element 'vertexSet'

## getting and setting information

# 'names' method for class 'vertexSet'
setMethod("names","vertexSet",function(x) x@.Data)
# data part is returned

# 'card' method is inherited from class 'vector'

## property checking

# 'isEmpty' method is inherited from class 'vector'

# comparison method for class 'vertexSet'
setMethod("areTheSame", c("vertexSet", "vertexSet"), function(x,y) setequal(x@.Data,y@.Data))
# a 'logical' value answering the question is returned

## extraction

# multi extractor method for class 'vertexSet'
setMethod("[","vertexSet",function(x,i,j=NA,drop=NA){
                            if(isEmpty(x)) x
                            else new("vertexSet",x@.Data[i])
                          } # end of function
         ) # end of setMethod
# a subset of vertices is extracted

# single extractor method for class 'vertexSet' gives a 'character' object

## typecasting

# typecasting method from 'vector' to 'vertexSet'
setAs("vector","vertexSet",function(from,to) new(to,from))
# just call the constructor with 'vector' input

# on the other hand, typecasting from 'vertexSet' to 'vector'
# is automatic and gives the '.Data' slot of the 'vertexSet' object

## operators

# union method for 'vertexSet' and 'vertexSet'
setMethod("+",c("vertexSet","vertexSet"),
                 function(e1,e2) new("vertexSet",union(e1@.Data,e2@.Data))
                ) # end of setMethod
# the union of the two vertex sets is returned

# intersection method for 'vertexSet' and 'vertexSet'
setMethod("*",c("vertexSet","vertexSet"),
                 function(e1,e2) new("vertexSet",intersect(e1@.Data,e2@.Data))
                ) # end of setMethod
# the intersection of the two vertex sets is returned

# asymmentric difference method for 'vertexSet' and 'vertexSet'
setMethod("-",c("vertexSet","vertexSet"),
                 function(e1,e2) new("vertexSet",setdiff(e1@.Data,e2@.Data))
                ) # end of setMethod
# the asymmetric difference of the two vertex sets is returned

### ----------------------------------edges------------------------------------

## construction and visualization

# constructor method for class 'undirectedEdge'
setMethod("initialize","undirectedEdge",
                 function(.Object,...){
                     id<-unique(as(unlist(list(...)),"integer")) # get unique 'integer' values
                     if(length(id)>0)
                         .Object@.Data<-id[id>0&!is.na(id)] # keep strictly positive values only
                         return(.Object)
                 } # end of function
                ) # end of SetMethod
# 'new("undirectedEdge",x,y)' returns an 'undirectedEdge' object consisting of a vector of
# unique strictly positive 'integer' numbers made from the unique elements of 'x' and 'y'

# constructor method for class 'directedEdge'
setMethod("initialize","directedEdge",
                 function(.Object,...){
                     Args <- list(...) # get the arguments
                     if(length(Args)==1) # handle single argument case
                         Args<-Args[[1]]
                     if(length(Args)==0) # empty object
                         .Object@.Data<-list()
                     else{ # process the arguments
                         if(is(Args,"list")){ # unlist and store order information
                             Ords<-numeric(0)
                             for(i in 1:length(Args))
                                 Ords<-c(Ords,rep(i,length(unlist(Args[[i]]))))
                             Args<-unlist(Args)
                         } else # order information is trivial
                             Ords<-seq(1,length(Args))
                         # make 'integer'
                         Args<-as(Args,"integer")
                         # keep strictly positive
                         Args<-Args[Args>0&!is.na(Args)]
                         Ords<-Ords[Args>0&!is.na(Args)]
                         # remove duplicates
                         dups<-duplicated(Args)
                         Args<-Args[!dups]
                         Ords<-Ords[!dups]
                         # construct the data part
                         res<-list()
                         len<-length(Args)
                         if(len>0){ # otherwise nothing to be done
                             preOrd<-0 # previous order
                             curVer<-0 # current vertex
                             for(i in 1:len){ # scan the arguments
                                 if(Ords[i]>preOrd){ # new vertex
                                     res<-c(res,list(Args[i]))
                                     preOrd<-Ords[i]
                                     curVer<-curVer+1
                                 } else # same vertex
                                     res[[curVer]]<-c(res[[curVer]],Args[i])
                             } # end of for
                         } # end of if
                     # fill object
                     .Object@.Data<-res
                     } # end of if-else
                     return(.Object)
                 } # end of function
                ) # end of setMethod
# returns a 'directedEdge' object whose data part is a 'list' of 'integer' vectors

# constructor method for class 'edgeList'
setMethod("initialize","edgeList",
                 function(.Object,...){
                     Args <- list(...) # get the arguments
                    if ( length(Args) >1) { # more than one argument
                        # just keep those of class 'edge'
                        Args<-Args[unlist(lapply(Args,"is",class2="edge"))]
                    }else if((length(Args) == 1)&&(!is(Args[[1]],"edge"))&&(is(Args[[1]],"list"))){
                        # single 'list' argument (not an 'edge')
                        # keep its elements of class 'edge', if any
                        if (length(Args[[1]])>0){ # at least one element
                            Args<-Args[[1]][unlist(lapply(Args[[1]],"is",class2="edge"))]
                        } else Args<-Args[[1]] # empty list
                    } # end of if-else
                     callNextMethod(.Object,Args) # create new object
                 } # end of function
                ) # end of setMethod
# returns an 'edgeList' object whose data part is a 'list'
# containing (only) the arguments of class 'edge'

# wrapper for undirected edge construction
u <- function(...) new("undirectedEdge",...)
# just call the right initialize method

# wrapper for directed edge construction
d <- function(...) new("directedEdge",...)
# just call the right initialize method

# wrapper for reverse construction
r <- function(...) {
           Args<-list(...)
           if(length(Args)==1) Args<-Args[[1]]
           new("directedEdge",rev(Args))
       } # end of function
# just call the right initialize method

# show method for class 'undirectedEdge'
setMethod("show","undirectedEdge",
                 function(object){
                     if(length(object)>1){
                         cat(object@.Data,sep="--",fill=TRUE)
                     } else if (length(object)==1){
                         cat(object@.Data,"<>",object@.Data,sep="",fill=TRUE)
                     } else{ # 'length(object)==0'
                         cat("--",fill=TRUE)
                    } # end of if-else if-else
                 } # end of function
                ) # end of setMethod
# for example '--' denotes the empty 'undirectedEdge'
# '1<>1' denotes a loop and '1--2' an ordinary 'undirectedEdge'
# while 1--2--3 denotes a three vertex 'undirectedEdge'

# 'showRel' method for class 'undirectedEdge' with reference to a 'vertexSet' object
setMethod("showRel",c(object="undirectedEdge",code="vertexSet"),
                 function(object,code){
                     if(length(object)>1){
                         cat(code@.Data[object@.Data],sep="--",fill=TRUE)
                     } else if (length(object)==1){
                         cat(code@.Data[object@.Data],"<>",code@.Data[object@.Data],sep="",fill=TRUE)
                     } else{ # 'length(object)==0'
                         cat("--",fill=TRUE)
                    } # end of if-else if-else
                 } # end of function
                ) # end of setMethod
# for example '--' denotes the empty 'undirectedEdge'
# 'a<>a' denotes a loop and 'a--b' an ordinary 'undirectedEdge'
# while a--b--c denotes a three vertex 'undirectedEdge'

# show method for class 'directedEdge'
setMethod("show","directedEdge",
                 function(object){
                     len<-length(object)
                     if(len>1) # proper directed edge
                         buf<-paste(unlist(lapply(object@.Data,"paste",collapse="--")),collapse="->")
                     else{ # loop or empty
                         buf<-"->"
                         if (len==1) buf<-paste(buf,paste(object[[1]],collapse="--"),sep="")
                     } # end of if-else
                     cat(buf,fill=TRUE)
                 } # end of function
                ) # end of setMethod
# for example '->' denotes the empty 'directedEdge'
# '->1' denotes a loop and '1->2' an ordinary 'directedEdge'
# while 1->2--3->4 denotes a general 'directedEdge'

# 'showRel' method for class 'directedEdge' with reference to a 'vertexSet' object
setMethod("showRel",c(object="directedEdge",code="vertexSet"),
                 function(object,code){
                     len<-length(object)
                     if(len>1) # proper directed edge
                         buf<-paste(unlist(lapply(object@.Data,function(x){
                                                         paste(code@.Data[x],collapse="--")
                                                         } # end of function
                                                 )),collapse="->")
                     else{ # loop or empty
                         buf<-"->"
                         if (len==1) buf<-paste(buf,paste(code@.Data[object[[1]]],collapse="--"),sep="")
                     } # end of if-else
                     cat(buf,fill=TRUE)
                 } # end of function
                ) # end of setMethod
# for example '->' denotes the empty 'directedEdge'
# '->a' denotes a loop and 'a->b' an ordinary 'directedEdge'
# while a->b--c->d denotes a general 'directedEdge'

# show method for class 'edgeList'
setMethod("show","edgeList",
                 function(object){
                     cat("{",fill=T)
                     lapply(object,"show")
                     cat("}",fill=T)
                 } # end of function
                ) # end of setMethod
# just show all edges in the list

# 'showRel' method for class 'edgeList'
setMethod("showRel",c(object="edgeList",code="vertexSet"),
                 function(object,code){
                     cat("{",fill=T)
                     lapply(object,"showRel",code=code)
                     cat("}",fill=T)
                 } # end of function
                ) # end of setMethod
# just show all edges in the list

## getting and setting information

# 'maxId' method for 'undirectedEdge'
setMethod("maxId","undirectedEdge", function(x) if(length(x)==0) return(0) else return(max(x)))
# gets the maximum numeric identifier of the edge

# 'maxId' operator for 'directedEdge'
setMethod("maxId","directedEdge", function(x) if(length(x)==0) return(0) else return(max(unlist(x))))
# gets the maximum numeric identifier of the edge

# 'maxId' operator for 'edgeList'
setMethod("maxId","edgeList", function(x) if(length(x)==0) return(0) else return(max(unlist(lapply(x,"maxId")))))
# gets the maximum numeric identifier of the multi-set

# recode method for class 'undirectedEdge' (from a vertex set to a new one)
setMethod("recode",c(object="undirectedEdge",src="vertexSet",dst="vertexSet"),
          function(object,src,dst) new("undirectedEdge",match(src[object],dst)))
# dropping all vertices that are not present in both of them

# recode method for class 'directedEdge' (from a vertex set to a new one)
setMethod("recode",c(object="directedEdge",src="vertexSet",dst="vertexSet"),
          function(object,src,dst){
            res<-lapply(object@.Data,function(x){
                match(src[x],dst)
            }) # end of lapply
            nas<-unlist(lapply(res,"is.na"))
            if(is(nas,"NULL"))
                return(new("directedEdge"))
            else
                return(new("directedEdge",res[!nas]))
          }) # end of setMethod
# dropping all vertices that are not present in both of them

# recode method for class 'edgeList' (from a vertex set to a new one)
setMethod("recode",c(object="edgeList",src="vertexSet",dst="vertexSet"),
          function(object,src,dst) new("edgeList",lapply(object,"recode",src=src,dst=dst))
         ) #Êend of setMethod
# every edge in the multi-set is dealt with by the appropriate recode method

# 'card' method for class 'directedEdge'
setMethod("card","directedEdge",
                 function(object,...){
                     length(unlist(object))
                 } # end of function
                ) # end of setMethod
# just unlist before taking length

## property checking

# 'isEmpty' method is inherited from class 'vector'

# 'isPresent' method for 'edge' in 'edgeList'
setMethod("isPresent",c("edge","edgeList"),
          function(el,ou) any(unlist(lapply(ou,"areTheSame",y=el)))
         ) # end of setMethod
# a 'logical' value answering the question is returned

# comparison method for class 'undirectedEdge'
setMethod("areTheSame", c("undirectedEdge", "undirectedEdge"), function(x,y) setequal(x@.Data,y@.Data))
# a 'logical' value answering the question is returned

# comparison method for class 'directedEdge'
setMethod("areTheSame", c("directedEdge", "directedEdge"),
                 function(x,y){
                     len<-length(x)
                     if(len==length(y)){ # maybe
                         res<-TRUE
                         i<-1
                         while(i<=len){
                             res<-res&&areTheSame(x[[i]],y[[i]])
                             i<-i+1
                         } # end of while
                     }else # no for sure
                         res<-FALSE
				     return(res)
                 } # end of function
               ) # end of setMethod
# a 'logical' value answering the question is returned

# comparison method for class 'edge' versus 'edge'
setMethod("areTheSame", c("edge", "edge"), function(x,y) return(FALSE))
# this will be used when comparing edges of different type

# comparison method for class 'edgeList'
setMethod("areTheSame", c("edgeList", "edgeList"),
          function(x,y){
            res<-(length(x)==length(y))
            if(res&&!isEmpty(x)){ # maybe and non-trivial
                m<-length(x)
                unmatched<-rep(TRUE,m) # refers to 'y'
                for(i in 1:m){ # match 'x[[i]]' in 'y'
                    found<-FALSE
                    for(j in which(unmatched)){ # try all unmatched 'y[[i]]'
                        if(areTheSame(x[[i]],y[[j]])){ # matched
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

## extraction

# multi extractor method for class 'undirectedEdge'
setMethod("[","undirectedEdge",function(x,i,j=NA,drop=NA) new("undirectedEdge",x@.Data[i]))
# an 'undirectedEdge' is extracted

# single extractor method for class 'undirectedEdge' gives an 'integer' object

# multi extractor method for class 'directedEdge'
setMethod("[","directedEdge",function(x,i,j=NA,drop=NA) new("directedEdge",x@.Data[i]))
# a 'directedEdge' is extracted

# single extractor method for class 'directedEdge'
setMethod("[[","directedEdge",function(x,i,j=NA,drop=NA) new("undirectedEdge",x@.Data[[i]]))
# an 'undirectedEdge' is extracted

# multiple extractor method for class 'edgeList'
setMethod("[","edgeList",function(x,i,j=NA,drop=NA)new("edgeList",x@.Data[i]))
# a sublist of edges is extracted

## typecasting

# typecasting method from 'vector' to 'undirectedEdge'
setAs("vector","undirectedEdge",function(from,to) new(to,from))
# just call the constructor with 'vector' input

# on the other hand, typecasting from "undirectedEdge" to "vector"
# is automatic and gives the '.Data' slot of the 'undirectedEdge' object

# typecasting method from 'vector' to 'directedEdge'
setAs("vector","directedEdge",function(from,to) new(to,from))
# just call the constructor with 'vector' input

# on the other hand, typecasting from 'directedEdge' to 'vector'
# is automatic and gives the '.Data' slot of the 'directedEdge' object

# typecasting from 'undirectedEdge' to 'directedEdge'
setAs("undirectedEdge","directedEdge",function(from,to) new(to,from))
# just call the constructor and exploit inheritance

# typecasting from 'directedEdge' to 'undirectedEdge'
setAs("directedEdge","undirectedEdge",function(from,to) new(to,from))
# just call the constructor and exploit inheritance

# typecasting from 'list' to 'edgeList'
setAs("list","edgeList",function(from,to) new(to,from))
# just call the constructor with 'list' input

# on the other hand, typecasting from "edgeList" to "list"
# is automatic and gives the '.Data' slot of the 'edgeList' object

## operators

# add operator for 'edgeList' and 'edge'
setMethod("+",c("edgeList","edge"),
                 function(e1,e2){
                     e1@.Data[[length(e1@.Data)+1]]<-e2
                     return(e1)
                 } # end of function
                ) # end of setMethod
# just add the edge to the list
setMethod("+",c("edge","edgeList"),function(e1,e2){e2+e1})
# and make the operation symmetric

# drop operator for 'edgeList' and 'edge' (not symmetric)
setMethod("-",c("edgeList","edge"),
                 function(e1,e2){
                     w<-which(as(unlist(lapply(e1,"areTheSame",y=e2)),"logical"))
                     if(length(w)==0){ # edge is not present in the list
                         return(e1)
                     } else{ # remove the first occurrence
                         return(e1[-w[1]])
                     } # end of if-else
                 } # end of function
                ) # end of setMethod
# look for the first copy of the edge in the list and remove it
