DAGGER <-
function(Maps, Linearize = TRUE, method = "QP", rescale = TRUE, MapFilename = NULL, GraphFilename = NULL, GraphDistances = FALSE) {

#############################################################
Merge <- function(M,G = NULL) {

#sort M
SortIdx <- sort(M[,2],index.return=TRUE)$ix
M <- M[SortIdx,]

if (dim(M)[1]==1) {
print("Linkage map must have at least two markers.")
} else if (is.null(G)) {
# no G was passed
#convert M to diG 
G <- list(Markers = list(), ZeroEdges = list(), ReverseEdges = list(), ForwardEdges = list(), Weights = list(), MarkerNames = M[,1], Vertices = rep(0,dim(M)[1]), Nvert = 1, Nmark = dim(M)[1])
G$Markers[[1]] <- 1
G$Vertices[1] <- 1
G$ReverseEdges[[1]] <- integer(0)
G$ForwardEdges[[1]] <- integer(0)
G$ZeroEdges[[1]] <- integer(0)
G$Weights[[1]] <- numeric(0)

for (i in 2:G$Nmark) {
    dist <- M[i,2] - M[i-1,2]
    if (dist==0) {
    #same linkage group
    G$Markers[[G$Nvert]] <- c(G$Markers[[G$Nvert]],i)
    G$Vertices[i] <- G$Nvert
    } else {
    #new linkage group
    G$ForwardEdges[[G$Nvert]] <- G$Nvert+1
    G$Weights[[G$Nvert]] <- dist
    G$Nvert <- G$Nvert + 1
    G$Markers[[G$Nvert]] <- i
    G$Vertices[i] <- G$Nvert
    G$ForwardEdges[[G$Nvert]] <- integer(0) #initialize to NULL
    G$ZeroEdges[[G$Nvert]] <- integer(0)
    G$Weights[[G$Nvert]] <- numeric(0)
    G$ReverseEdges[[G$Nvert]] <- G$Nvert-1
    } #end else
  } #end for

} else {
# G was passed, need to do merge sort on marker lists for M and G

MarkerCoding <- rep(0,dim(M)[1]) #this array gives the marker number in G for the markers in M
for (i in 1:length(MarkerCoding)) {
  pos <- which(G$MarkerNames==M[i,1])
  if (length(pos) == 0) {
    #this is a new marker
    G$Nmark <- G$Nmark+1
    MarkerCoding[i] <- G$Nmark
    G$MarkerNames[G$Nmark] <- M[i,1]
    G$Vertices[G$Nmark] <- 0
  } else {
    MarkerCoding[i] <- pos
  }
} #end for

#construct G2 from M
G2 <- list(Markers = list(), ZeroEdges = list(), ReverseEdges = list(), ForwardEdges = list(), Weights = list(), Nvert = 1, Nmark = length(MarkerCoding))
G2$Markers[[1]] <- MarkerCoding[1]
G2$Vertices[1] <- 1
G2$ReverseEdges[[1]] <- integer(0)
G2$ForwardEdges[[1]] <- integer(0)
G2$ZeroEdges[[1]] <- integer(0)
G2$Weights[[1]] <- numeric(0)

for (i in 2:G2$Nmark) {
    dist <- M[i,2] - M[i-1,2]
    if (dist==0) {
    #same linkage group
    G2$Markers[[G2$Nvert]] <- c(G2$Markers[[G2$Nvert]],MarkerCoding[i])
    } else {
    #new linkage group
    G2$ForwardEdges[[G2$Nvert]] <- G2$Nvert+1
    G2$Weights[[G2$Nvert]] <- dist
    G2$Nvert <- G2$Nvert + 1
    G2$Markers[[G2$Nvert]] <- MarkerCoding[i]
    G2$ReverseEdges[[G2$Nvert]] <- G2$Nvert-1
    G2$ForwardEdges[[G2$Nvert]] <- integer(0) #initialize to NULL
    G2$ZeroEdges[[G2$Nvert]] <- integer(0)
    G2$Weights[[G2$Nvert]] <- numeric(0)
    } #end else
  } #end for

#Now merge DAGS
#Add first eq class from G2
LinkFrom <- integer(0)
NumEquiv <- G$Nvert
Added <- integer(0)
W <- G2$Markers[[1]]  #markers in first vertex of G2

#which vertices in G contain markers in W
Vlist <- integer(0)
for (j in W) {Vlist <- union(Vlist,G$Vertices[j])}
Vlist <- setdiff(Vlist,0)

if (length(Vlist) > 0) {
for (j in Vlist) {
#partition bin if needed
  Vj <- G$Markers[[j]]
  Q <- intersect(Vj,W)
  if (setequal(Q,Vj)) {
     # Entire set Vj is contained in W
     LinkFrom <- c(LinkFrom,j)
  } else if (length(Q) > 0) {
  # need to partition Vj 
  NumEquiv <- NumEquiv + 1
  G$Markers[[j]] <- setdiff(Vj,Q)
  for (k in Q) {G$Vertices[k] <- NumEquiv}
  G$Markers[[NumEquiv]] <- Q
  G$ForwardEdges[[NumEquiv]] <- G$ForwardEdges[[j]]
  G$Weights[[NumEquiv]] <- G$Weights[[j]]
  G$ReverseEdges[[NumEquiv]] <- G$ReverseEdges[[j]]
  G$ZeroEdges[[NumEquiv]] <- c(G$ZeroEdges[[j]],j)
  for (k in G$ReverseEdges[[j]]) {
     u <- which(G$ForwardEdges[[k]]==j)
     G$ForwardEdges[[k]] <- c(G$ForwardEdges[[k]],rep(NumEquiv,length(u)))
     G$Weights[[k]] <- c(G$Weights[[k]],G$Weights[[k]][u])
  }
  for (k in G$ForwardEdges[[j]]) {
     G$ReverseEdges[[k]] <- c(G$ReverseEdges[[k]],NumEquiv)
  }  
  LinkFrom <- c(LinkFrom,NumEquiv)
  } #end if
  Added <- c(Added,Q)
}
} #end if length(Vlist) > 0

#now see if there are any markers left
Q <- setdiff(W,Added)
if (length(Q) > 0) {
  NumEquiv <- NumEquiv + 1
  G$Markers[[NumEquiv]] <- Q
  for (k in Q) {G$Vertices[k] <- NumEquiv}
  G$ForwardEdges[[NumEquiv]] <- integer(0)  #this initializes to NULL
  G$ReverseEdges[[NumEquiv]] <- integer(0)
  G$ZeroEdges[[NumEquiv]] <- integer(0)
  G$Weights[[NumEquiv]] <- numeric(0)
  LinkFrom <- c(LinkFrom,NumEquiv)
}    

#LinkFrom is set of vertices which contain markers in W.  They need zero edges to each other.
  n.LinkFrom <- length(LinkFrom)
  if (n.LinkFrom > 1) {
    for (j in 1:(n.LinkFrom-1)) {
       G$ZeroEdges[[LinkFrom[j]]] <- c(G$ZeroEdges[[LinkFrom[j]]],LinkFrom[(j+1):n.LinkFrom])
    }
  }


for (i in 2:G2$Nvert) {  
  W <- G2$Markers[[i]]
  LinkTo <- integer(0)
  Added <- integer(0)

  #which vertices in G contain markers in W
  Vlist <- integer(0)
  for (j in W) {Vlist <- union(Vlist,G$Vertices[j])}
  Vlist <- setdiff(Vlist,0)

  if (length(Vlist) > 0) {
  for (j in Vlist)  {
    #partition bin if needed
    Vj <- G$Markers[[j]]
    Q <- intersect(Vj,W)
    if (setequal(Q,Vj)) {
       LinkTo <- c(LinkTo,j)
    } else if (length(Q)>0) {
    # need to partition Vj
    NumEquiv <- NumEquiv + 1
    G$Markers[[j]] <- setdiff(Vj,Q)
    G$Markers[[NumEquiv]] <- Q
    for (k in Q) {G$Vertices[k] <- NumEquiv}
    G$ForwardEdges[[NumEquiv]] <- G$ForwardEdges[[j]]
    G$Weights[[NumEquiv]] <- G$Weights[[j]]
    G$ReverseEdges[[NumEquiv]] <- G$ReverseEdges[[j]]
    G$ZeroEdges[[NumEquiv]] <- c(G$ZeroEdges[[j]],j)
    for (k in G$ReverseEdges[[j]]) {
      u <- which(G$ForwardEdges[[k]]==j)
      G$ForwardEdges[[k]] <- c(G$ForwardEdges[[k]],rep(NumEquiv,length(u)))
      G$Weights[[k]] <- c(G$Weights[[k]],G$Weights[[k]][u])
    }  
    for (k in G$ForwardEdges[[j]]) {
      G$ReverseEdges[[k]] <- c(G$ReverseEdges[[k]],NumEquiv)
    }   
    LinkTo <- c(LinkTo,NumEquiv)
    } #end if
  Added <- c(Added,Q)
  }
  } #end if length(Vlist) > 0

  #now see if there are any markers left
  Q <- setdiff(W,Added)
  if (length(Q) > 0) {
    NumEquiv <- NumEquiv + 1
    G$Markers[[NumEquiv]] <- Q
    for (k in Q) {G$Vertices[k] <- NumEquiv}
    G$ForwardEdges[[NumEquiv]] <- integer(0)  #this initializes to NULL
    G$ReverseEdges[[NumEquiv]] <- integer(0)
    G$ZeroEdges[[NumEquiv]] <- integer(0)
    G$Weights[[NumEquiv]] <- numeric(0)
    LinkTo <- c(LinkTo,NumEquiv)
  }    
  
  for (j in LinkFrom) {
      G$ForwardEdges[[j]] <- c(G$ForwardEdges[[j]],LinkTo)
      G$Weights[[j]] <- c(G$Weights[[j]],rep(G2$Weights[[i-1]],length(LinkTo)))
  }
  for (j in LinkTo) {
      G$ReverseEdges[[j]] <- c(G$ReverseEdges[[j]],LinkFrom)
  }  

  #LinkTo is set of vertices which contain markers in W.  They need zero edges to each other.
  n.LinkTo <- length(LinkTo)
  if (n.LinkTo > 1) {
    for (j in 1:(n.LinkTo-1)) {
       G$ZeroEdges[[LinkTo[j]]] <- c(G$ZeroEdges[[LinkTo[j]]],LinkTo[(j+1):n.LinkTo])
    }
  }

  LinkFrom <- LinkTo
} #end for i
G$Nvert <- NumEquiv

} #end if length(G)==0

G  #return DAG
} #end function Merge

##############################################################

NumMaps <- length(Maps)
linkage.map.lengths <- rep(0,NumMaps)

if (NumMaps < 2) {
print("Must have at least two maps.")
} else {
   G <- Merge(Maps[[1]])
   linkage.map.lengths[1] <- Maps[[1]][nrow(Maps[[1]]),2]
   for (i in 2:NumMaps) {
      G <- Merge(Maps[[i]],G)
      linkage.map.lengths[i] <- Maps[[i]][nrow(Maps[[i]]),2]
   }

#Check for inconsistencies using SCC

#perform dfs on G_reverse
pre <- rep(0,G$Nvert)
post <- rep(0,G$Nvert)
clock <- 1

explore<-function(v) {
pre[v] <<- clock
clock <<- clock+1
for (w in G$ReverseEdges[[v]]) {
  if (pre[w] == 0) explore(w)
}
post[v] <<- clock
clock <<- clock+1
} #end explore

for (v in 1:G$Nvert) {
  if (pre[v] == 0) explore(v)
}

postsort <- sort(post,decreasing=TRUE,index.return=TRUE)$ix

#perform dfs on G
pre <- rep(0,G$Nvert)
post <- rep(0,G$Nvert)
scc <- rep(0,G$Nvert)
clock <- 1
sccnum <- 0

explore<-function(v) {
pre[v] <<- clock
scc[v] <<- sccnum
clock <<- clock+1
for (w in G$ForwardEdges[[v]]) {
  if (pre[w] == 0) explore(w)
}
post[v] <<- clock
clock <<- clock+1
} #end explore

for (v in postsort) {
  if (pre[v] == 0) {
    sccnum <- sccnum + 1
    explore(v)
  }
} # finished SCC identification

#construct list of names for the bins
BinNames <- character(0)

for (i in 1:G$Nvert) {
   MarkerList <- G$Markers[[i]]
   Name <- paste("\"",G$MarkerNames[MarkerList[1]],sep="")
   M <- length(MarkerList)
   if (M > 1) {
   for (j in 2:M) {
     Name <- paste(Name,G$MarkerNames[MarkerList[j]],sep="\\n")
   }
   }
   Name <- paste(Name,"\"",sep="")   
   BinNames[i] <- Name
} #end for i

if (sccnum < G$Nvert) {
# consensus graph has non-singleton SCCs

print("Consensus graph has inconsistencies.") 

if (!is.null(GraphFilename)) {

#write SCC to file
con <- file(GraphFilename,open="w")
writeLines("digraph SCC {",con)

for (k in 1:sccnum) {
  Vlist <- which(scc==k)
  if (length(Vlist) > 1) {
    #this scc is non-singleton
    for (i in Vlist) {
    count <- 0
    for (j in intersect(Vlist,G$ForwardEdges[[i]])) {
     count <- count + 1
     writeLines(paste(BinNames[i],"->",BinNames[j],"[label=",paste("\"",as.character(round(G$Weights[[i]][count],digits=2)),"\"",sep=""),"];"),con)
    } #end for j
    } #end for i
  } #end if
} #end for k
writeLines("}",con)
close(con)
} #end if GraphFilename

if (Linearize==FALSE) {FALSE} else {NULL}  #if Linearize=TRUE, return NULL

} else { 
#graph is acyclic (no inconsistencies)

if (Linearize==TRUE) {

#####################################
LinearizeMap<-function(G,method="LP") {
N <- G$Nvert
A <- Matrix(nrow=0,ncol=N,sparse=TRUE)  #adjacency matrix for nonzero edges
B <- Matrix(nrow=0,ncol=N,sparse=TRUE)  #adjacency matrix for zero edges
d <- numeric(0)
NormWeights <- numeric(0)
ZeroWeights <- numeric(0)

for (i in 1:N) {
  v <- rep(0,N)
  v[i] <- -1
  Edges <- G$ForwardEdges[[i]]
  Ne <- length(Edges)
  if (Ne > 0) {
  for (j in 1:Ne) {
    NormWeights <- c(NormWeights,length(G$Markers[[i]])*length(G$Markers[[Edges[j]]]))
    w <- v
    w[Edges[j]] <- 1
    A <- rBind(A,w)
    d <- c(d,G$Weights[[i]][j])
  } #end for j
  } #end if Ne > 0

  Zero.Edges <- G$ZeroEdges[[i]]
  Nz <- length(Zero.Edges)
  if (Nz > 0) {
   for (j in 1:Nz) {
     w <- v
     w[Zero.Edges[j]] <- 1
     B <- rBind(B,w)
     ZeroWeights <- c(ZeroWeights,length(G$Markers[[i]])*length(G$Markers[[Zero.Edges[j]]]))
   } #end for j
   } #end if Nz

} #end for i

n.zero <- nrow(B)
n.edge <- nrow(A)

if (method=="LP") {
 A.B <- rBind(A,B)
 d.0 <- c(d,rep(0,n.zero))
 LHS <- rBind(cBind(A,Matrix(0,nrow=n.edge,ncol=n.edge+n.zero,sparse=TRUE)),cBind(A.B,Diagonal(n.edge+n.zero)),cBind(-A.B,Diagonal(n.edge+n.zero)))
 RHS <- c(rep(0,n.edge),d.0,-d.0)
 f <- c(rep(0,N),NormWeights,ZeroWeights)
 dir <- rep(">=",3*n.edge+2*n.zero)
 H <- Rglpk_solve_LP(f,LHS,dir,RHS)

 if (H$status > 0) {
  print("Error in LPsolver.")
  } else {
  c.map <- H$solution[1:N]
  c.map - min(c.map)  #map starts at zero
 } #end if/else H
} else {
# method QP
 A <- rBind(c(1,rep(0,N-1)),A)
 d <- c(0,d)
 A.B <- rBind(A,B)
 d.0 <- c(d,rep(0,n.zero))
 Q <- Diagonal(n.edge+n.zero+1,c(1,NormWeights,ZeroWeights))
 QP.ans <- solve.QP(t(A.B)%*%Q%*%A.B,t(A.B)%*%Q%*%d.0,t(A),meq=1)
 c.map <- QP.ans$solution
 c.map - min(c.map)
} #end if method
} #end function
###################################################3

D <- LinearizeMap(G,method)  #map positions, guaranteed to be nonnegative

if (rescale == TRUE) {
  D <- D/max(D)*mean(linkage.map.lengths)
}

#sort distances in order of increasing map distance
SortIdx = sort(D,index.return=TRUE)$ix 

if (!is.null(MapFilename)) {
  #write linear map to file
  con <- file(MapFilename,open="w")

  count <- 0
  for (i in SortIdx) {
    count <- count + 1
    z <- chartr(old="\\n",new="  ",BinNames[i])
    writeLines(paste(count,as.character(round(D[i],digits=2)),substr(z,2,nchar(z)-1)),con)
  } #end for i
  close(con)
} #end if MapFilename

#construct data frame with linear map
NameArray <- character(0)
MapArray <- numeric(0)
count <- 0
for (i in SortIdx) {
  MapPos <- round(D[i],digits=2)
  for (j in G$Markers[[i]]) {
    count <- count + 1
    MapArray[count] <- MapPos
    NameArray[count] <- G$MarkerNames[j]
  } #end for j
} #end for i

Lmap <- data.frame(Name=NameArray,Position=MapArray)

}  #end if Linearize==TRUE

if (!is.null(GraphFilename)) {
  #need to write graph to file

  if (GraphDistances==FALSE) {
  #Remove extra edges and condense vertices so that only ordering is preserved

  #construct reachability list and prune edges
  ReachList <- list()
  for (v in postsort) {
    Vlist <- union(integer(0),G$ForwardEdges[[v]])  #removes multiple edges to same vertex
    G$ReverseEdges[[v]] <- union(integer(0),G$ReverseEdges[[v]])  #remove multiples from reverse edges
    NewEdges <- Vlist
    temp <- Vlist
    for (w in Vlist) {
       temp <- union(temp,ReachList[[w]])
       NewEdges <- setdiff(NewEdges,ReachList[[w]])
    } #end for w
    ReachList[[v]] <- temp
    Pruned <- setdiff(Vlist,NewEdges)  #these are forward edges that were pruned
    for (w in Pruned) {G$ReverseEdges[[w]] <- setdiff(G$ReverseEdges[[w]],v)}
    G$ForwardEdges[[v]] <- NewEdges 
  } #end for v

  #Note: Weights are no longer correct

  #Condense vertices
  for (i in 2:sccnum) {
    v <- postsort[i]
    Forward1 <- G$ForwardEdges[[v]]
    Reverse1 <- G$ReverseEdges[[v]]
    CheckList <- setdiff(postsort[(i-1):1],ReachList[[v]])

    for (w in CheckList) {
      Forward2 <- G$ForwardEdges[[w]]
      Reverse2 <- G$ReverseEdges[[w]]
  
      if (setequal(Forward1,Forward2) & setequal(Reverse1,Reverse2)) {
        #Eliminate vertex w, adding its markers to v 
        G$Markers[[v]] <- union(G$Markers[[v]],G$Markers[[w]])
        G$Vertices[G$Markers[[w]]] <- v
        G$Markers[[w]] <- integer(0)
        for (u in Forward2) {G$ReverseEdges[[u]]=setdiff(G$ReverseEdges[[u]],w)}
        for (u in Reverse2) {G$ForwardEdges[[u]]=setdiff(G$ForwardEdges[[u]],w)}
      } #end if setequal
   } #end for w 
  } #end for i

  #Note that G$Nvert is no longer correct for condensed graph
  }  #end if GraphDistances == FALSE

  #write graph to file
  con <- file(GraphFilename,open="w")
  writeLines("digraph DAG {",con)

  #construct labels for bins
  Bin2Bin <- rep(0,G$Nvert)  # this is needed to get right number of bins in condensed graph
  NumBins <- 0
  for (i in 1:G$Nvert) {
   MarkerList <- G$Markers[[i]]
   M <- length(MarkerList)
   if (M > 0) {  #this check is needed because some vertices have been condensed
     NumBins <- NumBins + 1
     Bin2Bin[i] <- NumBins
     label <- paste(NumBins," [label=\"",G$MarkerNames[MarkerList[1]],sep="")
     if (M > 1) {
       for (j in 2:M) {label <- paste(label,G$MarkerNames[MarkerList[j]],sep="\\n")}
     } #end if M > 1
   label <- paste(label,"\"];",sep="")   
   writeLines(label,con)
   } #end if M > 0
  } #end for i
 
  for (i in which(Bin2Bin > 0)) {
   count <- 0
   for (j in G$ForwardEdges[[i]]) {
     count <- count + 1
     if (GraphDistances==TRUE) {
      writeLines(paste(Bin2Bin[i],"->",Bin2Bin[j],"[label=",paste("\"",as.character(round(G$Weights[[i]][count],digits=2)),"\"",sep=""),"];"),con)
     } else {
      writeLines(paste(Bin2Bin[i],"->",Bin2Bin[j],";"),con)
     } #end if/else
   } #end for j
  if (GraphDistances==TRUE) {
   for (j in G$ZeroEdges[[i]]) {
     writeLines(paste(Bin2Bin[i],"->",Bin2Bin[j],"[label=0,arrowhead=none];"),con)
   } #end for j
  } #end if GraphDistances
  } #end for i
  writeLines("}",con)
  close(con)

}# end if GraphFilename 

if (Linearize==TRUE) {
  Lmap  #return linearized map as data frame
} else {
  TRUE #return TRUE indicating graph was acyclic
}
} #end else sscnum < N
} #end if NumMaps < 2
}  #end dagger  

