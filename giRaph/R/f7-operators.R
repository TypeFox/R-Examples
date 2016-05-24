## f7-operators.R --- 
## Author          : Jens Henrik Badsberg, Claus Dethlefsen, Luca la Rocca
## Created On      : Wed Oct 27 10:59:27 2004
## Last Modified By: Luca La Rocca
## Last Modified On: Sun Feb 26 18:36:00 2006
## Update Count    : 159
## Status          : Unknown, Use with caution!
######################################################

## ----------------------------------------------
## representation and vertex set
## ----------------------------------------------

# operation 'incidenceList' + 'vertexSet'
setMethod("+",signature=c("incidenceList","vertexSet"),
          function(e1,e2) {
            e1@V<-e1@V+e2
            return(e1)
          } # end of function
         ) # end of setMethod
# note that it is not necessary to recode the edges

# operation 'incidenceList' - 'vertexSet'
setMethod("-",signature=c("incidenceList","vertexSet"),
          function(e1,e2) {
            Vnames<-names(e1)
            keepV<-match(setdiff(Vnames,names(e2)),Vnames)
            return(e1[keepV])
          } #Êend of function
         ) # end of setMethod
# note that in general some edges are dropped

# operation 'incidenceList' * 'vertexSet'
setMethod("*",signature=c("incidenceList","vertexSet"),
          function(e1,e2) {
            Vnames<-names(e1)
            keepV<-match(intersect(Vnames,names(e2)),Vnames)
            return(e1[keepV])
          } #Êend of function
         ) # end of setMethod
# note that in general some edges are dropped

# operation 'incidenceMatrix' + 'vertexSet'
setMethod("+",signature=c("incidenceMatrix","vertexSet"),
          function(e1,e2) {
            newids<-setdiff(names(e2),names(e1))
            addN<-length(newids)
            if (addN>0){
              oldC<-ncol(e1@.Data)
              e1@.Data <- cbind(e1@.Data,matrix(0,nrow(e1@.Data),addN))
              colnames(e1@.Data) <- c(colnames(e1@.Data)[1:oldC],newids)
            } # end of if
            return(e1)
          } # end of function
         ) #Êend of setMethod
# note that the new vertices will be isolated vertices

# operation 'incidenceMatrix' - 'vertexSet'
setMethod("-",signature=c("incidenceMatrix","vertexSet"),
          function(e1,e2) {
            Vnames<-names(e1)
            Vkeep<-match(setdiff(Vnames,names(e2)),Vnames)
            if(isEmpty(Vkeep))
                Ekeep<-rep(FALSE,nrow(e1))
            else
                Ekeep<-apply(e1@.Data,1,function(x){sum(x[-Vkeep])==0})
            colnames(e1@.Data)<-NULL
            e1@.Data<-matrix(e1@.Data[Ekeep,Vkeep],sum(Ekeep),length(Vkeep))
            colnames(e1@.Data)<-Vnames[Vkeep]
            return(e1)
          } # end of function
         ) # end of setMethod
# note that in general some edges are dropped

# operation 'incidenceMatrix' * 'vertexSet'
setMethod("*",signature=c("incidenceMatrix","vertexSet"),
          function(e1,e2) {
            Vnames<-names(e1)
            Vkeep<-match(intersect(Vnames,names(e2)),Vnames)
            if(isEmpty(Vkeep))
                Ekeep<-rep(FALSE,nrow(e1))
            else
                Ekeep<-apply(e1@.Data,1,function(x){sum(x[-Vkeep])==0})
            colnames(e1@.Data)<-NULL
            e1@.Data<-matrix(e1@.Data[Ekeep,Vkeep],sum(Ekeep),length(Vkeep))
            colnames(e1@.Data)<-Vnames[Vkeep]
            return(e1)
          } # end of function
         ) # end of setMethod
# note that in general some edges are dropped

# operation 'adjacencyList' + 'vertexSet'
setMethod("+",signature=c("adjacencyList","vertexSet"),
          function(e1,e2) {
            Vnames<-names(e1)
            Vadd<-setdiff(names(e2),Vnames)
            e1@.Data<-c(e1@.Data,rep(list(list()),length(Vadd)))
            names(e1)<-c(Vnames,Vadd)
            return(e1)
          } # end of function
         ) # end of setMethod
# note that the new vertices will be isolated vertices

# operation 'adjacencyList' - 'vertexSet'
setMethod("-",signature=c("adjacencyList","vertexSet"),
          function(e1,e2) {
            Vnames<-names(e1)
            return(e1[match(setdiff(Vnames,names(e2)),Vnames)])
          } # end of function
        ) # end of setMethod
# note that in general some edges are dropped

# operation 'adjacencyList' * 'vertexSet'
setMethod("*",signature=c("adjacencyList","vertexSet"),
          function(e1,e2) {
            Vnames<-names(e1)
            return(e1[match(intersect(Vnames,names(e2)),Vnames)])
          } # end of function
        ) # end of setMethod
# note that in general some edges are dropped

# operation 'adjacencyMatrix' + 'vertexSet'
setMethod("+",signature=c("adjacencyMatrix","vertexSet"),
          function(e1,e2) {
            newids<-setdiff(names(e2),names(e1))
            addN<-length(newids)
            if (addN>0){
              oldN<-nrow(e1@.Data)
              e1@.Data <- cbind(e1@.Data,matrix(0,oldN,addN))
              e1@.Data <- rbind(e1@.Data,matrix(0,addN,oldN+addN))
              rownames(e1@.Data) <- c(rownames(e1@.Data)[1:oldN],newids)
              colnames(e1@.Data) <- rownames(e1@.Data)
            } # end if
            return(e1)
          } # end of function
         ) # end of setMethod
# note that the new vertices will be isolated vertices

# operation 'adjacencyMatrix' - 'vertexSet'
setMethod("-",signature=c("adjacencyMatrix","vertexSet"),
          function(e1,e2) {
            Vnames<-names(e1)
            keepV<-match(setdiff(Vnames,names(e2)),Vnames)
            newN<-length(keepV)
            rownames(e1@.Data)<-NULL
            colnames(e1@.Data)<-NULL
            e1@.Data<-matrix(e1@.Data[keepV,keepV],newN,newN)
            rownames(e1@.Data)<-Vnames[keepV]
            colnames(e1@.Data)<-rownames(e1@.Data)
            return(e1)
          } # end of function
         ) # end of setMethod
# note that in general some edges are dropped

# operation 'adjacencyMatrix' * 'vertexSet'
setMethod("*",signature=c("adjacencyMatrix","vertexSet"),
          function(e1,e2) {
            Vnames<-names(e1)
            keepV<-match(intersect(Vnames,names(e2)),Vnames)
            newN<-length(keepV)
            rownames(e1@.Data)<-NULL
            colnames(e1@.Data)<-NULL
            e1@.Data<-matrix(e1@.Data[keepV,keepV],newN,newN)
            rownames(e1@.Data)<-Vnames[keepV]
            colnames(e1@.Data)<-rownames(e1@.Data)
            return(e1)
          } # end of function
         ) # end of setMethod
# note that in general some edges are dropped

## ----------------------------------------------
## representation and edge
## ----------------------------------------------

# operation 'incidenceList' + 'edge'
setMethod("+",signature=c("incidenceList","edge"),
          function(e1,e2) {
          if(maxId(e2)<=length(e1@V)) e1@E<-e1@E+e2
          return(e1)
          } # end of function
         ) # end of setMethod
# the edge is added to the multi-set of edges (if meaningful for the vertex set)

# operation 'incidenceList' - 'edge'
setMethod("-",signature=c("incidenceList","edge"),
              function(e1,e2){
                e1@E<-e1@E-e2
                return(e1)
              } # end of function
            ) # end of setMethod
# the edge is removed from the multi-set of edges (if present)

# operation 'incidenceMatrix' + 'undirectedEdge'
setMethod("+",signature=c("incidenceMatrix","undirectedEdge"),
          function(e1,e2) {
            if((length(e2)>0)&&max(e2)<=ncol(e1)){ # add
                e1@.Data<-rbind(e1@.Data,rep(0,ncol(e1)))
                e1@.Data[nrow(e1@.Data),e2@.Data]<-1
            } # end of if
            return(e1)
          } # end of function
         ) # end of setMethod
# the edge is added to the multi-set of edges
# (if non-empty and meaningful for the vertex set)

# operation 'incidenceMatrix' - 'undirectedEdge'
setMethod("-",signature=c("incidenceMatrix","undirectedEdge"),
          function(e1,e2) {
            if((length(e2)>0)&&(max(e2)<=ncol(e1))&&(nrow(e1)>0)){
                eline<-rep(0,ncol(e1))
                eline[e2@.Data]<-1
                where<-match(T,apply(e1@.Data,1,function(x){all(x==eline)}))
                if(!is.na(where)) e1@.Data<-matrix(e1@.Data[-where,],nrow(e1)-1,ncol(e1))
            } # end of if
            return(e1)
          } # end of function
         ) # end of setMethod
# the edge is removed from the multi-set of edges
# (if non-empty, meaningful for the vertex set and present)

# operation 'incidenceMatrix' + 'directedEdge'
setMethod("+",signature=c("incidenceMatrix","directedEdge"),
          function(e1,e2) {
            if((length(e2)>1)&&(maxId(e2)<=ncol(e1))){ # add
                e1@.Data<-rbind(e1@.Data,rep(0,ncol(e1)))
                e1@.Data[nrow(e1@.Data),unlist(e2)]<-rep(1:length(e2),unlist(lapply(e2@.Data,length)))
            } # end of if
            return(e1)
          } # end of function
         ) # end of setMethod
# the edge is added to the multi-set of edges
# (if proper and meaningful for the vertex set)

# operation 'incidenceMatrix' - 'directedEdge'
setMethod("-",signature=c("incidenceMatrix","directedEdge"),
          function(e1,e2) {
            if((length(e2)>1)&&(maxId(e2)<=ncol(e1))&&(nrow(e1)>0)){
                eline<-rep(0,ncol(e1))
                eline[unlist(e2)]<-rep(1:length(e2),unlist(lapply(e2@.Data,length)))
                where<-match(T,apply(e1@.Data,1,function(x){all(x==eline)}))
                if(!is.na(where)) e1@.Data<-matrix(e1@.Data[-where,],nrow(e1)-1,ncol(e1))
            } # end of if
            return(e1)
          } # end of function
         ) # end of setMethod
# the edge is removed from the multi-set of edges
# (if non-empty, meaningful for the vertex set and present)

# operation 'adjacencyList' + 'undirectedEdge'
setMethod("+",signature=c("adjacencyList","undirectedEdge"),
          function(e1,e2) {
          if((length(e2)==2)&&(max(e2)<=length(e1))){ # add
            one<-e2@.Data[1]
            two<-e2@.Data[2]
            e1@.Data[[one]]$ne<-c(e1@.Data[[one]]$ne,two)
            e1@.Data[[two]]$ne<-c(e1@.Data[[two]]$ne,one)
          }else if((length(e2)==1)&&(e2<=length(e1))){ # loop
            e1@.Data[[e2]]$ne<-c(e1@.Data[[e2]]$ne,e2)
          } #Êend of if-else if
          return(e1)
          } # end of function
         ) # end of setMethod
# the edge is added to the multi-set of edges
# (if ordinary and meaningful for the vertex set)

# operation 'adjacencyList' - 'undirectedEdge'
setMethod("-",signature=c("adjacencyList","undirectedEdge"),
          function(e1,e2) {
            if((length(e2))==1&&(e2<=length(e1))){ # loop
                try<-match(e2,e1@.Data[[e2]]$ne) # first match
                if(!is.na(try)) e1@.Data[[e2]]$ne<-e1@.Data[[e2]]$ne[-try]
            }else if((length(e2)==2)&&(max(e2)<=length(e1))){
                one<-e2@.Data[1]
                two<-e2@.Data[2]
                try<-match(two,e1@.Data[[one]]$ne) # first match
                if(!is.na(try)){
                    e1@.Data[[one]]$ne<-e1@.Data[[one]]$ne[-try]
                    e1@.Data[[two]]$ne<-e1@.Data[[two]]$ne[-match(one,e1@.Data[[two]]$ne)]
                } #Êend of if
            } # end of if-else if
            return(e1)
          } # end of function
         ) # end of setMethod
# the edge is removed from the multi-set of edges (if present)

# operation 'adjacencyList' + 'directedEdge'
setMethod("+",signature=c("adjacencyList","directedEdge"),
          function(e1,e2) {
          if((length(e2)==2)&&(card(e2)==2)&&(maxId(e2)<=length(e1))){ # add
            one<-e2@.Data[[1]]
            two<-e2@.Data[[2]]
            e1@.Data[[one]]$ch<-c(e1@.Data[[one]]$ch,two)
            e1@.Data[[two]]$pa<-c(e1@.Data[[two]]$pa,one)
          } #Êend of if
          return(e1)
          } # end of function
         ) # end of setMethod
# the edge is added to the multi-set of edges
# (if ordinary and meaningful for the vertex set)

# operation 'adjacencyList' - 'directedEdge'
setMethod("-",signature=c("adjacencyList","directedEdge"),
          function(e1,e2) {
            if((length(e2)==2)&&(card(e2)==2)&&(maxId(e2)<=length(e1))){
                one<-e2@.Data[[1]]
                two<-e2@.Data[[2]]
                try<-match(two,e1@.Data[[one]]$ch) # first match
                if(!is.na(try)){
                    e1@.Data[[one]]$ch<-e1@.Data[[one]]$ch[-try]
                    e1@.Data[[two]]$pa<-e1@.Data[[two]]$pa[-match(one,e1@.Data[[two]]$pa)]
                } #Êend of if
            } # end of if
            return(e1)
          } # end of function
         ) # end of setMethod
# the edge is removed from the multi-set of edges (if present)

# operation 'adjacencyMatrix' + 'undirectedEdge'
setMethod("+",signature=c("adjacencyMatrix","undirectedEdge"),
          function(e1,e2) {
          if((length(e2)==2)&&(max(e2)<=nrow(e1))){ # add
            one<-e2@.Data[1]
            two<-e2@.Data[2]
            e1@.Data[one,two]<-1
            e1@.Data[two,one]<-1
          } #Êend of if
          return(e1)
          } # end of function
         ) # end of setMethod
# the edge is added to the set of edges (if it is ordinary,
# it is not a loop, and it is meaningful for the vertex set)

# operation 'adjacencyMatrix' - 'undirectedEdge'
setMethod("-",signature=c("adjacencyMatrix","undirectedEdge"),
          function(e1,e2) {
          if((length(e2)==2)&&(max(e2)<=nrow(e1))){ # maybe remove
            one<-e2@.Data[1]
            two<-e2@.Data[2]
            if(e1@.Data[one,two]&&e1@.Data[two,one]){ # edge is there
              e1@.Data[one,two]<-0
              e1@.Data[two,one]<-0
            } # end of if "edge is there"
          } #Êend of if "maybe remove"
          return(e1)
          } # end of function
         ) # end of setMethod
# the edge is removed from the set of edges

# operation 'adjacencyMatrix' + 'directedEdge'
setMethod("+",signature=c("adjacencyMatrix","directedEdge"),
          function(e1,e2) {
          if((length(e2)==2)&&(card(e2)==2)&&(maxId(e2)<=nrow(e1))){ # add
            one<-e2@.Data[[1]]
            two<-e2@.Data[[2]]
            e1@.Data[one,two]<-1
          } #Êend of if
          return(e1)
          } # end of function
         ) # end of setMethod
# the edge is added to the set of edges (if it is ordinary,
# and it is meaningful for the vertex set)

# operation 'adjacencyMatrix' - 'directedEdge'
setMethod("-",signature=c("adjacencyMatrix","directedEdge"),
          function(e1,e2) {
          if((length(e2)==2)&&(card(e2)==2)&&(maxId(e2)<=nrow(e1))){ # maybe remove
            one<-e2@.Data[[1]]
            two<-e2@.Data[[2]]
            if(e1@.Data[one,two]&&!e1@.Data[two,one]) # edge is ther
              e1@.Data[one,two]<-0
          } #Êend of if ("maybe remove")
          return(e1)
          } # end of function
         ) # end of setMethod
# the edge is removed from the set of edges

## ----------------------------------------------
## graph and vertex set
## ----------------------------------------------

# operation 'anyGraph' + 'vertexSet'
setMethod("+",signature=c("anyGraph","vertexSet"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList+e2
            return(e1)
          } # end of function
         ) # end of setMethod
# isolated vertices are added to the graph

# operation 'generalGraph' + 'vertexSet'
setMethod("+",signature=c("generalGraph","vertexSet"),
          function(e1,e2){
            if(!isEmpty(e1@incidenceList)){ # incidence list available
              e1@incidenceList<-e1@incidenceList+e2
              if(!isEmpty(e1@incidenceMatrix)) # incidence matrix too
                e1@incidenceMatrix<-e1@incidenceMatrix+e2
            }else{ # incidence matrix only
              e1@incidenceMatrix<-e1@incidenceMatrix+e2
            } # end of if else
            return(e1)
          } # end of function
         ) # end of setMethod
# isolated vertices are added to the graph

# operation 'multiGraph' + 'vertexSet'
setMethod("+",signature=c("multiGraph","vertexSet"),
          function(e1,e2){
            if(!isEmpty(e1@incidenceList)){ # incidence list available
              e1@incidenceList<-e1@incidenceList+e2
              if(!isEmpty(e1@incidenceMatrix)) # incidence matrix too
                e1@incidenceMatrix<-e1@incidenceMatrix+e2
              if(!isEmpty(e1@adjacencyList)) # adjacency list too
                e1@adjacencyList<-e1@adjacencyList+e2
            }else if(!isEmpty(e1@incidenceMatrix)){ # no incidence list, but incidence matrix available
              e1@incidenceMatrix<-e1@incidenceMatrix+e2
              if(!isEmpty(e1@adjacencyList)) # adjacency list too
                e1@adjacencyList<-e1@adjacencyList+e2
            }else{ # adjacency list only
              e1@adjacencyList<-e1@adjacencyList+e2
            } # end of if-else if-else
            return(e1)
          } # end of function
         ) # end of setMethod
# isolated vertices are added to the graph

# operation 'simpleGraph' + 'vertexSet'
setMethod("+",signature=c("simpleGraph","vertexSet"),
          function(e1,e2){
            if(!isEmpty(e1@incidenceList)){ # incidence list available
              e1@incidenceList<-e1@incidenceList+e2
              if(!isEmpty(e1@incidenceMatrix)) # incidence matrix too
                e1@incidenceMatrix<-e1@incidenceMatrix+e2
              if(!isEmpty(e1@adjacencyList)) # adjacency list too
                e1@adjacencyList<-e1@adjacencyList+e2
              if(!isEmpty(e1@adjacencyMatrix)) # adjacency matrix too
                e1@adjacencyMatrix<-e1@adjacencyMatrix+e2
            }else if(!isEmpty(e1@incidenceMatrix)){ # no incidence list, but incidence matrix available
              e1@incidenceMatrix<-e1@incidenceMatrix+e2
              if(!isEmpty(e1@adjacencyList)) # adjacencyList too
                e1@adjacencyList<-e1@adjacencyList+e2
              if(!isEmpty(e1@adjacencyMatrix)) # adjacency matrix too
                e1@adjacencyMatrix<-e1@adjacencyMatrix+e2
            }else if(!isEmpty(e1@adjacencyList)){ # no incidence list, nor incidence matrix, but adjacency list
              e1@adjacencyList<-e1@adjacencyList+e2
              if(!isEmpty(e1@adjacencyMatrix)) # adjacency matrix too
                e1@adjacencyMatrix<-e1@adjacencyMatrix+e2
            }else{ # adjacency matrix only
              e1@adjacencyMatrix<-e1@adjacencyMatrix+e2
            }# end of if-else if-else if-else
            return(e1)
          } # end of function
         ) # end of setMethod
# isolated vertices are added to the graph

# operation 'anyGraph' - 'vertexSet'
setMethod("-",signature=c("anyGraph","vertexSet"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList-e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph

# operation 'generalGraph' - 'vertexSet'
setMethod("-",signature=c("generalGraph","vertexSet"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList-e2
            e1@incidenceMatrix<-e1@incidenceMatrix-e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph

# operation 'multiGraph' - 'vertexSet'
setMethod("-",signature=c("multiGraph","vertexSet"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList-e2
            e1@incidenceMatrix<-e1@incidenceMatrix-e2
            e1@adjacencyList<-e1@adjacencyList-e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph

# operation 'simpleGraph' - 'vertexSet'
setMethod("-",signature=c("simpleGraph","vertexSet"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList-e2
            e1@incidenceMatrix<-e1@incidenceMatrix-e2
            e1@adjacencyList<-e1@adjacencyList-e2
            e1@adjacencyMatrix<-e1@adjacencyMatrix-e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph

# operation 'anyGraph' * 'vertexSet'
setMethod("*",signature=c("anyGraph","vertexSet"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList*e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph

# operation 'generalGraph' * 'vertexSet'
setMethod("*",signature=c("generalGraph","vertexSet"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList*e2
            e1@incidenceMatrix<-e1@incidenceMatrix*e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph

# operation 'multiGraph' * 'vertexSet'
setMethod("*",signature=c("multiGraph","vertexSet"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList*e2
            e1@incidenceMatrix<-e1@incidenceMatrix*e2
            e1@adjacencyList<-e1@adjacencyList*e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph

# operation 'simpleGraph' * 'vertexSet'
setMethod("*",signature=c("simpleGraph","vertexSet"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList*e2
            e1@incidenceMatrix<-e1@incidenceMatrix*e2
            e1@adjacencyList<-e1@adjacencyList*e2
            e1@adjacencyMatrix<-e1@adjacencyMatrix*e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph

## ----------------------------------------------
## graph and edge
## ----------------------------------------------

# operation 'anyGraph' + 'edge'
setMethod("+",signature=c("anyGraph","edge"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList+e2
            return(e1)
          } # end of function
         ) # end of setMethod
# an edge is added to any graph

# operation 'generalGraph' + 'edge'
setMethod("+",signature=c("generalGraph","edge"),
          function(e1,e2){
            if(is(e2,"undirectedEdge")&&card(e2)>0||is(e2,"directedEdge")&&length(e2)>1){ # add edge
              e1@incidenceList<-e1@incidenceList+e2
              e1@incidenceMatrix<-e1@incidenceMatrix+e2
            } # otherwise do nothing
            return(e1)
          } # end of function
         ) # end of setMethod
# a proper directed or undirected edge is added to a general graph

# operation 'multiGraph' + 'edge'
setMethod("+",signature=c("multiGraph","edge"),
          function(e1,e2){
            if(is(e2,"undirectedEdge")&&card(e2)>0&&card(e2)<3||
               is(e2,"directedEdge")&&length(e2)>1&&card(e2)<3){ # add edge
              e1@incidenceList<-e1@incidenceList+e2
              e1@incidenceMatrix<-e1@incidenceMatrix+e2
              e1@adjacencyList<-e1@adjacencyList+e2
            } # otherwise do nothing
            return(e1)
          } # end of function
         ) # end of setMethod
# a proper non-hyper directed or undirected edge is added to a multi graph

# operation 'simpleGraph' + 'undirectedEdge'
setMethod("+",signature=c("simpleGraph","undirectedEdge"),
          function(e1,e2){
            if(card(e2)==2&&!isPresent(e2,e1)){ # add edge
              arrow<-new("directedEdge",e2@.Data[1],e2@.Data[2])
              worra<-new("directedEdge",e2@.Data[2],e2@.Data[1])
              # update incidence list
              e1@incidenceList<-e1@incidenceList-arrow
              e1@incidenceList<-e1@incidenceList-worra
              e1@incidenceList<-e1@incidenceList+e2
              # update incidence matrix
              e1@incidenceMatrix<-e1@incidenceMatrix-arrow
              e1@incidenceMatrix<-e1@incidenceMatrix-worra
              e1@incidenceMatrix<-e1@incidenceMatrix+e2
              # update adjacency list
              e1@adjacencyList<-e1@adjacencyList-arrow
              e1@adjacencyList<-e1@adjacencyList-worra
              e1@adjacencyList<-e1@adjacencyList+e2
              # update adjacency matrix
              e1@adjacencyMatrix<-e1@adjacencyMatrix+e2
            } # otherwise do nothing
            return(e1)
          } # end of function
         ) # end of setMethod
# a proper non-hyper undirected edge is added to a simple graph (if not already present)

# operation 'simpleGraph' + 'directedEdge'
setMethod("+",signature=c("simpleGraph","directedEdge"),
          function(e1,e2){
            if(length(e2)==2&&card(e2)==2&&!isPresent(e2,e1)){ # maybe add edge
              line<-new("undirectedEdge",e2@.Data[[1]],e2@.Data[[2]])
              worra<-new("directedEdge",e2@.Data[[2]],e2@.Data[[1]])
              if(!isPresent(line,e1)){ #Êindeed add edge
                if(isPresent(worra,e1)){ # gives an undirected edge
                  # update incidence list
                  e1@incidenceList<-e1@incidenceList-worra
                  e1@incidenceList<-e1@incidenceList+line
                  # update incidence matrix
                  e1@incidenceMatrix<-e1@incidenceMatrix-worra
                  e1@incidenceMatrix<-e1@incidenceMatrix+line
                  # update adjacency list
                  e1@adjacencyList<-e1@adjacencyList-worra
                  e1@adjacencyList<-e1@adjacencyList+line
                }else{
                  # update incidence list, incidence matrix and adjacency list
                  e1@incidenceList<-e1@incidenceList+e2
                  e1@incidenceMatrix<-e1@incidenceMatrix+e2
                  e1@adjacencyList<-e1@adjacencyList+e2
                } # end of if-else
                # update adjacency matrix
                e1@adjacencyMatrix<-e1@adjacencyMatrix+e2
              } # otherwise not really to be added
            } # otherwise do nothing
            return(e1)
          } # end of function
         ) # end of setMethod
# a proper non-hyper directed edge is added to a simple graph (if not already present)
# possibly originating an undirected edge

# operation 'anyGraph' - 'edge'
setMethod("-",signature=c("anyGraph","edge"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList-e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph

# operation 'generalGraph' - 'edge'
setMethod("-",signature=c("generalGraph","edge"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList-e2
            e1@incidenceMatrix<-e1@incidenceMatrix-e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph

# operation 'multiGraph' - 'edge'
setMethod("-",signature=c("multiGraph","edge"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList-e2
            e1@incidenceMatrix<-e1@incidenceMatrix-e2
            e1@adjacencyList<-e1@adjacencyList-e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph

# operation 'simpleGraph' - 'edge'
setMethod("-",signature=c("simpleGraph","edge"),
          function(e1,e2){
            e1@incidenceList<-e1@incidenceList-e2
            e1@incidenceMatrix<-e1@incidenceMatrix-e2
            e1@adjacencyList<-e1@adjacencyList-e2
            e1@adjacencyMatrix<-e1@adjacencyMatrix-e2
            return(e1)
          } # end of function
         ) # end of setMethod
# vertices and related edges are removed from the graph
