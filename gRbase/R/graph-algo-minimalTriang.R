# R function returning a minimal triangulation of an undirected graph by the recursive thinning method
# Author: Clive Bowsher

# Inputs:
# uG: graphNEL representation of the undirected graph

# Output: returns a minimal triangulation TuG of the undirected graph uG; the mcwh method is used to
#         obtain the initial triangulation prior to applying the Recursive Thinning algorithm of Olesen & Madsen 2002

# date originated: 23.03.09
# checked:  using the DAGs of Fig.1,2 of Olesen & Madsen 2002, Fig1a and Fig2 of Leimer93 (CGB 04.03.09).
#			using TestMinTriang2.R (CGB 08.03.09), see comments on that file
# issues: clearly the triangulation is not guranteed to be a minimUM one, but this is not necessary for our purposes here


minimalTriang <- function(object, tobject=triangulate(object), result=NULL, details=0){
    UseMethod("minimalTriang")
}

minimalTriang.default <- function(object, tobject=triangulate(object), result=NULL, details=0){
    cls <- match.arg(class( object ),
                     c("graphNEL","matrix","dgCMatrix"))
    if (is.null( result ))
        result <- cls

    switch(cls,
           "graphNEL" ={tt<-.minimalTriang(object, TuG=tobject, details=details) },
           "dgCMatrix"=,
           "matrix"   ={ #FIXME: minimalTriang: Not sure if this is correct...
               object2 <- as(object,  "graphNEL")
               tobject2<- as(tobject, "graphNEL")
               tt<-.minimalTriang(object2, TuG=tobject2, details=details)
           })
    as(tt, result)
}


minimalTriangMAT <- function(amat, tamat=triangulateMAT(amat), details=0){
  as.adjMAT(.minimalTriang(as(amat, "graphNEL"), TuG=as(tamat, "graphNEL"), details=details))
}


.minimalTriang <- function(uG, TuG=triangulate(uG, method="mcwh"), details=0) {

  removed <- 0

  if (is.triangulated(uG)) # no triangulation of G is needed
    {
      if (details>0)
        cat(sprintf("Graph is already triangulated\nNumber of edges removed in recursive thinngs: %i\n", removed))
      return(uG)
    }
  else
    {
      ##TuG <- triangulate(uG, method="mcwh")	# calls C implementation of mcwh method

      ## Soren modification
      uGmat  <- as(uG,"matrix")
      di     <- dimnames(uGmat)
      TuGmat <- as(TuG,"matrix")
      TuGmat <- TuGmat[di[[1]],di[[1]]]
      TT <- TuGmat - uGmat # edges added in triangulatn in adjacency representatn
      ## ends here...

      ## TT <- as(TuG,"matrix") - as(uG,"matrix")	# edges added in triangulatn in adjacency representatn

      TT <- as(TT,"graphNEL")
      Tn <- edgeList(TT)

      if (details>0)
        cat(sprintf("Number of edges added in triangulation step: %i \n", length(Tn)))

      Rn <- Tn
      exclT <- new("graphNEL", nodes = graph::nodes(uG), edgeL = vector("list", length = 0))

      repeat{
        if (length(Rn) == 0) {	# if R' is empty so is TT' and so cannot execute the for loop below
          break
        }
        for(i in 1:length(Rn)){
          neX <- graph::adj(TuG,Rn[[i]][1])
          neY <- graph::adj(TuG,Rn[[i]][2])
          neXY <- neX[[1]][neX[[1]] %in% neY[[1]]] 	#select elements of neX in neY to obtain intersection

##           sg <<- subGraph(neXY,TuG)
##           print(sg)
          if(is.complete(graph::subGraph(neXY,TuG)))
            {
              TuG   <- graph::removeEdge(Rn[[i]][1],Rn[[i]][2],TuG)	# directly updates TuG
              exclT <- graph::addEdge(Rn[[i]][1],Rn[[i]][2],exclT)	# keep track of excluded edges as a graph
              removed <- removed + 1
            }
        }
        if (length(edgeList(exclT)) == 0){
          break
        }
        ## Soren modification
        uGmat  <- as(uG,"matrix")
        di <- dimnames(uGmat)
        TuGmat <- as(TuG,"matrix")
        TuGmat <- TuGmat[di[[1]],di[[1]]]
        TT <- TuGmat - uGmat # edges added in triangulatn in adjacency representatn
        ## ends here...
        ##TT <- as(TuG,"matrix") - as(uG,"matrix")		# recompute the triangulation edges retained

        TT <- as(TT,"graphNEL")
        Tn <- edgeList(TT)

        ## form vector of nodes corresponding to the edge list (possibly with repeated elements)
        exclT <- c(edgeList(exclT),recursive=TRUE)

        Rn <- vector("list", length = 0)
        j <- 0
        for(k in 1:length(Tn)) {
          if((Tn[[k]][1] %in% exclT) | (Tn[[k]][2] %in% exclT)) {
            j <- j+1
            Rn[[j]] <- Tn[[k]]
          }
        }

        ## NOTE: TuG already updated at line 43
        exclT <- new("graphNEL", nodes = graph::nodes(uG), edgeL = vector("list", length = 0))
      }
      if (details>0)
        cat(sprintf("Number of edges removed in recursive thinngs: %i\n", removed))
      if (!(is.triangulated(TuG))){
        cat("WARNING: minimal triangulation step failed\n")
      }

      return(TuG)
    }
}
