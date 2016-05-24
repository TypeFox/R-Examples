## f8-interfaces.R --- 
## Author          : Jens Henrik Badsberg, Claus Dethlefsen, Luca La Rocca
## Created On      : Wed Nov 03 19:22:57 2004
## Last Modified By: Claus Dethlefsen
## Last Modified On: Thu Apr 19 13:49:17 2007
## Update Count    : 99
## Status          : Unknown, Use with caution!
######################################################

# Typecasting from 'simpleGraph' to 'dg.simple.graph'
#setAs("simpleGraph", "dg.simple.graph",
#      function(from, to) {
#        iList <- incidenceList(from)
#        new("dg.simple.graph", vertex.names = names(iList),
#            edge.list = lapply(iList@E, function(i) unlist(i)),
#            oriented  = lapply(iList@E, function(i) is(i, "directedEdge"))
#        ) #Êend of new
#      } # end of funcion
#     ) # end of setAs
# a new 'dg.simple.graph' object is created

# Typecasting from 'simpleGraph' to 'dg.graph'
#setAs("simpleGraph", "dg.graph", function(from, to) as(as(from,"dg.simple.graph"), to))
# a new 'dg.graph' object is created

# Typecasting from 'multiGraph' to 'dg.simple.graph'
# via 'simpleGraph' i.e. dropping loops and parallel edges
#setAs("multiGraph", "dg.simple.graph", function(from, to) as(as(from,"simpleGraph"), to))
# a new 'dg.simple.graph' object is created

# Typecasting from 'multiGraph' to 'dg.graph'
# via 'simpleGraph' i.e. dropping loops and parallel edges
#setAs("multiGraph", "dg.graph", function(from, to) as(as(from,"dg.simple.graph"), to))
# a new 'dg.graph' object is created

# Typecasting from 'generalGraph' to 'dg.graph'
#setAs("generalGraph", "dg.graph",
#      function(from, to) {
#        iList <- incidenceList(from)
#        hyper<-unlist(lapply(iList@E,function(e) card(e)>2))
#        hyperList<-iList@E[hyper] # put hyper edges aside
#        iList@E<-iList@E[!hyper] # consider ordinary edges
#        iList<-incidenceList(new("simpleGraph",incidenceList=iList)) # drop loops and parallel edges
#        edges<-lapply(iList@E,function(i) unlist(i))
#        orien<-lapply(iList@E,function(i) is(i,"directedEdge"))
#        nOrd<-length(edges) # remember the number of ordinary edges
#        if(length(hyperList)>0){ # consider hyper edges
#          iList@E<-new("edgeList") # expand all of them
#          for(it in seq(1,length(hyperList))){
#            e<-hyperList[[it]]
#            if(is(e,"undirectedEdge")) # undirected edge
#              for (i in seq(1,length(e)-1))
#                for (j in seq(i+1,length(e))){
#                  iList@E<-iList@E+new("undirectedEdge",e[[i]],e[[j]])
#            }else{ # directed edge
#              for (i in seq(1,length(e))){ # for all components
#                if(length(e[[i]])>1) # not a singleton component
#                  for(j in seq(1,length(e[[i]])-1)) # add undirected edges
#                    for (k in seq(j+1,length(e[[i]]))){
#                      iList@E<-iList@E+new("undirectedEdge",e[[i]][[j]],e[[i]][[k]])
#                } # end of if-for-for
#                if(i<length(e)) # not the last component
#                  for (j in e[[i]]) # add directed edges
#                    for (k in e[[i+1]]){
#                      iList@E<-iList@E+new("directedEdge",j,k)
#                } # end of if-for-for
#              } #Êend of for all components                                      
#            } # end of if-else
#          } #Êend of for
#        iList<-incidenceList(new("simpleGraph",incidenceList=iList)) # drop parallel edges
#        for(it in 1:length(iList@E)){ # add edges resulting from expansion when possible
#            x<-unlist(iList@E[[it]])
#            if(nOrd==0||!any(unlist(lapply(edges[1:nOrd],function(y) all(x==y)||all(x==y[2:1]))))){
#              edges[[length(edges)+1]]<-x
#              orien[[length(orien)+1]]<-is(iList@E[[it]],"directedEdge")
#            } # end of if
#          } # end of for
#        } # end of if
#                  
#        res<-simpleGraphToGraph(new("dg.simple.graph", vertex.names = names(iList),
#                                    edge.list = edges, oriented  = orien))
#
#        if(nOrd<length(edges)) # dashed expanded hyperedges
#          for(i in seq(nOrd+1,length(edges)))
#            res@edgeList[[i]]@dash<-"-"
#           
#        return(res)
#      } # end of function
#     ) # end of setAs
# a new 'dg.graph' object is created with dashed expanded hyper edges

# Typecasting from 'anyGraph' to 'dg.graph'
# via 'generalGraph' i.e. keeping non-empty undirected edges and proper directed edges only
#setAs("anyGraph", "dg.graph",function(from, to) as(as(from,"generalGraph"), to))
# a new 'dg.graph' object is created with dashed expanded hyper edges

# Dynamic graphical representation of simple graphs via package 'dynamicGraph'
#setMethod("dynamic.Graph", signature(object = "simpleGraph"),
#          function(object, ...){
#            require(dynamicGraph)
#            dg(simpleGraphToGraph(as(object,"dg.simple.graph"), ...),modelObject = new("dg.Model"),...)
#                               } #Êend of function
#         ) # end of setMethod
# based on typecasting to class 'dg.simple.graph'

# Dynamic graphical representation of multi graphs via package 'dynamicGraph'
#setMethod("dynamic.Graph", signature(object = "multiGraph"),
#          function(object, ...){
#            require(dynamicGraph)
#            dg(simpleGraphToGraph(as(object,"dg.simple.graph"), ...),modelObject = new("dg.Model"),...)
#                               } #Êend of function
#         ) # end of setMethod
# based on typecasting to class 'dg.simple.graph'

# Dynamic graphical representation of general graphs via package 'dynamicGraph'
#setMethod("dynamic.Graph", signature(object = "generalGraph"),
#          function(object, ...){
#            require(dynamicGraph)
#            dg(as(object,"dg.graph"),modelObject = new("dg.Model"),...)
#                               } #Êend of function
#         ) # end of setMethod
# based on typecasting to class 'dg.graph'

# Dynamic graphical representation of any graph via package 'dynamicGraph'
#setMethod("dynamic.Graph", signature(object = "anyGraph"),
#          function(object, ...){
#            require(dynamicGraph)
#            dg(as(object,"dg.graph"),modelObject = new("dg.Model"),...)
#                               } #Êend of function
#         ) # end of setMethod
# based on typecasting to class 'dg.graph'

## Bioconductors 'graph' (or, rather graphNEL)
##
## setAs("simpleGraph","graphNEL",function(from,to){
##   require(graph)
##   Nodes <- names(from)
##   G <- incidenceList(from)
##   Edges <- G@edgeList
##   
##   n <- length(Edges)
##   ndirected <- sum(unlist(lapply(Edges,is,"directedEdge")))
##   nundirected<- n - ndirected
## 
##   if (ndirected>=nundirected) edgemode <- "directed"
##   else edgemode <- "undirected"
## 
##   if (edgemode=="directed") Edges <- lapply(Edges,d)
##   else Edges <- lapply(Edges,u)
## 
##   G@edgeList <- new("edgeList",Edges)
##   A <- as(G,"adjacencyList")
##   A <- lapply(A,function(x) {
##     if (!is(x,"reverseEdge")) return(unlist(x))
##   })
## 
##   A <- lapply(A,function(i) list(edges=Nodes[i]))
##   
##   new(to,nodes=Nodes,edgeL=A,edgemode=edgemode)
##   })
## 
## setAs("graphNEL","simpleGraph",function(from,to){
##   require(graph)
## 
##   ## apparently not working
##   
##   N <- nodes(from)
##   A <- edges(from)
##   mode <- edgemode(p)
## 
##   A <- lapply(A, function(i) match(i, N))
##   
##   if (mode=="directed")
##     A <- lapply(A,function(x) lapply(x,d))
##   else
##     A <- lapply(A,function(x) lapply(x,u))
## 
##   I <- new("adjacencyList",A)
##   X <- as(I,"adjacencyMatrix")
## 
##   new(to,adjacencyMatrix=X)
##   })

# Make S3 informal 'mathgraph' class an S4 formal class
setOldClass("mathgraph")
# note that loading package 'mathgraph' is not necessary
# to work with this class, but it gives a better show method

# Typecasting from 'simpleGraph' to 'mathgraph'
setAs("simpleGraph","mathgraph",
      function(from,to){
        
        I <- as(adjacencyMatrix(from), "incidenceList")
        M <- matrix( unlist(lapply(I@E,unlist)), ncol=2, byrow=TRUE)
        colnames(M) <- c("e1","e2")
        directed <- unlist(lapply(I@E,is,"directedEdge"))
        
        X <- adjacencyMatrix(from)
        rsum <- apply(X@.Data,2,sum)
        csum <- apply(X@.Data,1,sum)
        idx.iso <- (1:ncol(X))[rsum==0&csum==0]
        if (length(idx.iso)>0) {
          M <- rbind(M,matrix(rep(idx.iso,2),nrow=length(idx.iso)))
        } # end of if
        
        attr(M,"directed") <- directed
        class(M) <- "mathgraph"
        return(M)
      } # end of funcion
    ) # end of setAs
# returns a 'mathgraph' object

# Typecasting from 'mathgraph' to 'simpleGraph'
setAs("mathgraph","simpleGraph",
      function(from,to){
        if(nrow(from)==0){ # no vertices => empty simple graph
          return(new("simpleGraph"))
        }else{ # there are vertices
          iList<-new("incidenceList",V=new("vertexSet",1:max(from)))
          from<-unclass(from) # so that 'from[i,j]' works
          for(i in 1:nrow(from)){ # add all edges
            if(attributes(from)$directed[i])
              iList@E<-iList@E+new("directedEdge",from[i,1],from[i,2])
            else
              iList@E<-iList@E+new("undirectedEdge",from[i,1],from[i,2])
          } # end of for
          return(new("simpleGraph",incidenceList=iList))
        } # end of if-else
      } # end of function
     ) # end of setAs
# returns a 'simpleGraph' object

# Static graphical representation of simple graphs via package 'mathgraph'
setMethod("display","simpleGraph",
          function(x,...){
            require(mathgraph)
            plot(as(x,"mathgraph"),...)
          } # end of function
         ) # end of setAs
# based on typecasting to class 'mathgraph'

# Static graphical representation of multi graphs via package 'mathgraph'
setMethod("display","multiGraph",
          function(x,...){
            display(as(x,"simpleGraph"),...)
          } # end of function
         ) # end of setAs
# implemented via typecasting to class 'simpleGraph'

# Static graphical representation of general graphs via package 'mathgraph'
setMethod("display","generalGraph",
          function(x,...){
            display(as(x,"simpleGraph"),...)
          } # end of function
         ) # end of setAs
# implemented via typecasting to class 'simpleGraph'

# Static graphical representation of any graph via package 'mathgraph'
setMethod("display","anyGraph",
          function(x,...){
            display(as(x,"simpleGraph"),...)
          } # end of function
         ) # end of setAs
# implemented via typecasting to class 'simpleGraph'
