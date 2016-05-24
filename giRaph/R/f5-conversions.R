## f5-conversions.R --- 
## Author          : Jens Henrik Badsberg, Claus Dethlefsen, Luca La Rocca
## Created On      : Fri Jun 24 10:55:00 2005
## Last Modified By: Luca La Rocca
## Last Modified On: Sat Feb 25 16:23:00 2006
## Update Count    : 17
## Status          : Unknown, Use with caution!
######################################################

## -----------------------------------------------------------
## CONVERSIONS between representations
## -----------------------------------------------------------
#f\t|  X  |  A  |  I  |  G
#X  |  *  |  /  |  /  |  +
#A  |  /  |  *  |  /  | (+) 
#I  |  +  |  +  |  *  | (+)
#G  |  +  | (+) |  +  |  *

## Note that conversions from a more general representation
## to a simpler representation silently (i.e. no warning)
## drop all edges not "available" in the new representation.

## ---------------------------------------------------
## 'incidenceList' -> ...
## ---------------------------------------------------

### typecasting from 'incidenceList' to 'incidenceMatrix'
setAs("incidenceList", "incidenceMatrix", function(from,to) {
  
  n <- card(from)$v # number of vertices
  m <- card(from)$e # number of edge occurrences
  
  I <- matrix(0, ncol = n, nrow = m)

  counter<-1 # counts edges kept (plus one)
  if (m>0) { # there are edges
    for (e in 1:m) { # for all edges
      edge <- from@E[[e]]
      if (is(edge,"undirectedEdge")&&!isEmpty(edge)){ # undirected and non-empty edge
          I[counter, edge@.Data ] <- 1
          counter <- counter + 1
      } else if (is(edge,"directedEdge")&&length(edge)>1) { # directed and proper edge
          I[counter, unlist(edge@.Data) ] <- rep(1:length(edge),unlist(lapply(edge@.Data,length)))
          counter <- counter + 1
      } ## else do nothing (empty undirected edges, improper directed edges and other types of edge are ignored)
    } # end of for
  } # end of if
  
  if(counter==1){ # no edges (kept)
    colnames(I) <- names(from)
  }else{
    I<-matrix(I[seq(1,counter-1),],nrow=counter-1,ncol=n,dimnames=list(NULL,names(from)))
  } # end of if-else
  
  return(new(to,I))
})
# an 'incidenceMatrix' object is returned

### typecasting from 'incidenceList' to 'adjacencyList'
setAs("incidenceList","adjacencyList", function(from,to) {

  n <- card(from)$v # number of vertices
  m <- card(from)$e # number of edge occurrences

  A <- new(to,id=names(from))

  if (m>0) { # there are edges
    for (h in 1:m) { # for all edges
      edge <- from@E[[h]]
      if (is(edge,"undirectedEdge")){
          q <- card(edge)
          if (q==2) { # (non-empty) ordinary (non-loop) edge
              A@.Data[[edge@.Data[1]]]$ne[length(A@.Data[[edge@.Data[1]]]$ne)+1] <- edge@.Data[2]
              A@.Data[[edge@.Data[2]]]$ne[length(A@.Data[[edge@.Data[2]]]$ne)+1] <- edge@.Data[1]
          } else if (q==1){ # loop
              A@.Data[[edge@.Data[1]]]$ne[length(A@.Data[[edge@.Data[1]]]$ne)+1] <- edge@.Data[1]
          } # end of "undirected edge"
      } else if (is(edge,"directedEdge")&&length(edge)>1) {
          if ( length(edge)>=2 && card(edge)<=2 ) { # proper ordinary edge
            A@.Data[[edge@.Data[[1]]]]$ch[length(A@.Data[[edge@.Data[[1]]]]$ch)+1] <- edge@.Data[[2]]
            A@.Data[[edge@.Data[[2]]]]$pa[length(A@.Data[[edge@.Data[[2]]]]$pa)+1] <- edge@.Data[[1]]
          } # end of "directed edge"
      } ## else do nothing
    } # end of 'for (h in 1:m)'
  } # end of 'if(m>0)'
  
  return(A)
})
# an 'adjacencyList' object is returned

### typecasting from 'incidenceList' to 'adjacencyMatrix'
setAs("incidenceList", "adjacencyMatrix", function(from, to) {

  n <- card(from)$v # number of vertices
  m <- card(from)$e # number of edge occurrences

  X <- matrix(0, ncol = n, nrow = n)
  dimnames(X) <- list(names(from),names(from))
  
  if (m>0){ # there are edges
    for (i in 1:m) { # for all edges
      edge <- from@E[[i]]
      if (is(edge,"undirectedEdge")){
        if (card(edge)==2){
          X[edge@.Data[[1]],edge@.Data[[2]]] <- 1
          X[edge@.Data[[2]],edge@.Data[[1]]] <- 1
        } # end of "undirected edge"
      } else if (is(edge,"directedEdge")&&length(edge)>1){
        if (card(edge)==2){
          X[edge@.Data[[1]],edge@.Data[[2]]] <- 1
        } # end of "directed edge"
      } ## else do nothing (hyperedges are ignored)
        ## note that multiple edges are reduced to a single edge
        ## and that 1->2, 2<-1 becomes 1-2 and so on.
    } # end of 'for (h in 1:m)'
  } # end of 'if(m>0)'

  return(new(to,X))
})
# an 'adjacencyMatrix' object is returned

## ---------------------------------------------------
## 'incidenceMatrix' -> ...
## ---------------------------------------------------

### typecasting from 'incidenceMatrix' to 'incidenceList'
setAs("incidenceMatrix", "incidenceList", function(from,to) {

  n <- card(from)$v # number of vertices
  m <- card(from)$e # number of edge occurrences

  if(n==0) return(new(to)) # empty representation

  E <- as(rep(NA,m),"list")

  if (m>0) {
      for (i in 1:m) {
          edge <- from@.Data[i,]
          edgeorder <- edge[!edge==0]
          if (all(edgeorder==1))
              E[[i]] <- new("undirectedEdge",(1:n)[!edge==0])
          else if (!any(duplicated(edgeorder)))
              E[[i]] <- new("directedEdge",(1:n)[!edge==0][sort.list(edgeorder)])
          else {
              res <- rep( list( list()), max(edge) )
              for (j in unique(edgeorder)) {res[[j]] <- (1:n)[edge==j]}
              E[[i]] <- new("directedEdge",res)
          } # end of if-else if-else
      } # end of 'for (i in 1:m)'
  } # end of 'if (m>0)'
  
  return(new(to,V=names(from),E=E))

})
# an 'incidenceList' object is returned

### typecasting from 'incidenceMatrix' to 'adjacencyList'
setAs("incidenceMatrix", "adjacencyList", function(from,to){

  n <- card(from)$v # number of vertices
  m <- card(from)$e # number of edge occurrences
  
  A <- new(to,id=names(from))

  if(m>0){ # there are edges
      for (i in 1:m) { # for all edges
        edge <- from@.Data[i,]
        edgeorder <- edge[!edge==0] 
        q <- length(edgeorder)
        if (all(edgeorder==1)) { # "undirected edge"
              edge <- (1:n)[edge!=0]
              if (q==2) { # (non-empty) ordinary (non-loop) edge
                A@.Data[[ edge[1] ]]$ne <- c(A@.Data[[ edge[1] ]]$ne, edge[2])
                A@.Data[[ edge[2] ]]$ne <- c(A@.Data[[ edge[2] ]]$ne, edge[1])
              } else if (q==1) { # loop
                A@.Data[[ edge[1] ]]$ne <- c(A@.Data[[ edge[1] ]]$ne, edge[1])
              } # end of if-else
        } else { # "directed edge"
              if (q==2) { # proper ordinary edge
                start <- (1:n)[edge==1]
                end   <- (1:n)[edge==2]
                A@.Data[[ start ]]$ch <- c(A@.Data[[ start ]]$ch, end)
                A@.Data[[ end ]]$pa   <- c(A@.Data[[ end ]]$pa, start)
              } # end of if
        } # end of if-else
      } # end of 'for (i in 1:m)'
  } # end of 'if(m>0)'
  
  return((A))
})
# an 'adjacencyList' object is returned

### typecasting from 'incidenceMatrix' to 'adjacencyMatrix'
setAs("incidenceMatrix", "adjacencyMatrix", function(from,to) {

  n <- card(from)$v # number of vertices
  m <- card(from)$e # number of edge occurrences
  
  X <- matrix(0, nrow = n, ncol = n)
  dimnames(X) <- list(names(from),names(from))
  
  if (m>0) { # there are edges
    for (i in 1:m) { # for all edges
      edge <- from@.Data[i,]
      if (sum(edge)==3) # "directed edge"
          X[edge==1,edge==2] <- 1
      else if (sum(edge)==2) { # "undirected edge"
        idx <- (1:n)[edge==1]
        X[idx[1],idx[2]] <- 1
        X[idx[2],idx[1]] <- 1
      } # end of if-else if
    } # end of 'for (i in 1:m)'
  } # end of 'if(m>0)'
  
  return(new(to,X))
})
# an 'adjacencyMatrix' object is returned

## ---------------------------------------------------
## 'adjacencyList' -> ...
## ---------------------------------------------------

### typecasting from 'adjacencyList' to 'incidenceList'
setAs("adjacencyList","incidenceList", function(from,to) {

  n <- card(from)$v # number of vertices
  
  if(n==0) # empty
      return(new(to))
  else{ # not empty
      E <- list()
      for(i in 1:n){ # for all vertices
          a <- from@.Data[[i]]
          if (length(a$ne[a$ne>=i])>0) for (j in a$ne[a$ne>=i])
            E[[length(E)+1]] <- new("undirectedEdge",i,j)
          if (length(a$ch)>0) for (j in a$ch)
            E[[length(E)+1]] <- new("directedEdge",i,j)
      } # end of for
  return(new(to,V=names(from),E=E))
  } # end of if-else

})
# an 'incidenceList' object is returned

### typecasting from 'adjacencyList' to 'incidenceMatrix'
setAs("adjacencyList", "incidenceMatrix", function(from,to) {

  n <- card(from)$v # number of vertices

  aux <- list()
  
  if(n>0){ # not empty
      for(i in 1:n){
          a <- from@.Data[[i]]
          if (length(a$ne[a$ne>=i])>0) for (j in a$ne[a$ne>=i])
            aux[[length(aux)+1]] <- c(i,j) # undirected edge
          if (length(a$ch)>0) for (j in a$ch)
            aux[[length(aux)+1]] <- list(i,j) # directed edge
      } # end of for
  } # end of if

  m<-length(aux)
  
  I <- matrix(0, ncol = n, nrow = m)
  colnames(I) <- names(from)

  if(m>0){ # there are edges
      for(e in 1:m){ # for all edges
          edge<-aux[[e]]
          if(is(edge,"list")){ # directed edge
              I[e,edge[[1]]]<-1
              I[e,edge[[2]]]<-2
          }else{ # undirected edge
              I[e,edge]<-1
          } # end of if-else
      } # end of for
  } # end of if
  
  return(new(to,I))
})
# an 'incidenceMatrix' object is returned

### typecasting from 'adjacencyList' to 'adjacencyMatrix'
setAs("adjacencyList","adjacencyMatrix", function(from, to) {

  n <- card(from)$v # number of vertices

  X <- matrix(0, ncol = n, nrow = n)
  dimnames(X) <- list( names(from), names(from))

  if(n>0){ # not empty
      for(i in 1:n){
          a <- from@.Data[[i]]
          if (length(a$ne[a$ne>i])>0)
              for (j in a$ne[a$ne>i]){ # undirected edges
                  X[i,j]<-1
                  X[j,i]<-1
              } # end of for
          if (length(a$ch)>0)
              for (j in a$ch) # directed edges
                  X[i,j]<-1
      } # end of 'for (i in 1:n)'
  } # end of 'if(n>0)'
  
  return(new(to,X))
  })
# an 'adjacencyMatrix' object is returned

## ---------------------------------------------------
## 'adjacencyMatrix' -> ...
## ---------------------------------------------------

### typecasting from 'adjacencyMatrix' to 'incidenceList'
setAs("adjacencyMatrix","incidenceList", function(from,to) {

  n <- card(from)$v # number of vertices

  E<-list()

  if(n>1){ # maybe there are edges
  for(i in seq(1,n-1)){
          for(j in seq(i+1,n)){
              if(from@.Data[i,j]){ # edge
                  if(from@.Data[j,i]){ # undirected
                      E[[length(E)+1]]<-new("undirectedEdge",i,j)
                  }else{ # directed edge
                      E[[length(E)+1]]<-new("directedEdge",i,j)
                  } # end of if-else
              } else if (from@.Data[j,i]){ # directed
                  E[[length(E)+1]]<-new("directedEdge",j,i)
              } # end of if-else
          } # end of for (j)
      } # end of for (i)
  } # end of if
  
  new(to,V=names(from),E=E)
})
# an 'incidenceList' object is returned

### typecasting from 'adjacencyMatrix' to 'incidenceMatrix'
setAs("adjacencyMatrix", "incidenceMatrix", function(from,to) {
  
  n <- card(from)$v # number of vertices
  
  aux <- list()
  
  if(n>1){ # maybe there are edges
  for(i in seq(1,n-1)){
          for(j in seq(i+1,n)){
              if(from@.Data[i,j]){ # edge
                  if(from@.Data[j,i]){ # undirected
                      aux[[length(aux)+1]]<-c(i,j)
                  }else{ # directed edge
                      aux[[length(aux)+1]]<-list(i,j)
                  } # end of if-else
              } else if (from@.Data[j,i]){ # directed
                  aux[[length(aux)+1]]<-list(j,i)
              } # end of if-else
          } # end of for (j)
      } # end of for (i)
  } # end of if

  m<-length(aux)
  
  I <- matrix(0, ncol = n, nrow = m)
  colnames(I) <- names(from)

  if(m>0){ # there are edges
      for(e in 1:m){ # for all edges
          edge<-aux[[e]]
          if(is(edge,"list")){ # directed edge
              I[e,edge[[1]]]<-1
              I[e,edge[[2]]]<-2
          }else{ # undirected edge
              I[e,edge]<-1
          } # end of if-else
      } # end of for
  } # end of if
  
  return(new(to,I))
})
# an 'incidenceMatrix' object is returned

### typecasting from 'adjacencyMatrix' to 'adjacencyList'
setAs("adjacencyMatrix","adjacencyList", function(from, to) {

  n <- card(from)$v # number of vertices
  
  A <- new(to,id=names(from))
  
  if(n>1){ # maybe there are edges
  for(i in seq(1,n-1)){
          for(j in seq(i+1,n)){
              if(from@.Data[i,j]){ # edge
                  if(from@.Data[j,i]){ # undirected
                      A@.Data[[i]]$ne[length(A@.Data[[i]]$ne)+1]<-j
					  A@.Data[[j]]$ne[length(A@.Data[[j]]$ne)+1]<-i
                  }else{ # directed edge
				      A@.Data[[i]]$ch[length(A@.Data[[i]]$ch)+1]<-j
				      A@.Data[[j]]$pa[length(A@.Data[[j]]$pa)+1]<-i
                  } # end of if-else
              } else if (from@.Data[j,i]){ # directed
				      A@.Data[[i]]$pa[length(A@.Data[[i]]$pa)+1]<-j
				      A@.Data[[j]]$ch[length(A@.Data[[j]]$ch)+1]<-i
              } # end of if-else
          } # end of for (j)
      } # end of for (i)
  } # end of if

  return((A))
  })
# an 'adjacencyList' object is returned
